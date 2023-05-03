#if !defined(CLING) || defined(ROOTCLING)

#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ITSMFTSimulation/Hit.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITStracking/IOUtils.h"
#include "DataFormatsParameters/GRPObject.h"

#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"

#include <gsl/gsl>
#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TSystemDirectory.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "CommonDataFormat/RangeReference.h"
#include "StrangenessTracking/StrangenessTracker.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "StrangenessTracking/StrangenessTracker.h"
#include "GPUCommonArray.h"
#include "DetectorsBase/Propagator.h"

#include "strTrackingUtils.h"

#endif

using namespace o2;
using namespace o2::framework;

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using Cascade = o2::dataformats::Cascade;

using MCTrack = o2::MCTrack;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using CompClusterExt = o2::itsmft::CompClusterExt;
using ITSCluster = o2::BaseCluster<float>;
using Vec3 = ROOT::Math::SVector<double, 3>;
using StrangeTrack = o2::dataformats::StrangeTrack;

const int motherPDG = 3334; // Omega
const int firstDaughterPDG = 3122; // pdg of the V0  = Lambda
const int secondDaughterPDG = 321; // pdg of the bachelor = K
const int firstV0dauPDG = 2212;    // pdg of the V0 daughter = p
const int secondV0dauPDG = 211;    // pdg of the V0 daughter = Pi

std::array<int, 2> matchCascToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec, Cascade &casc);
std::array<int,2> getBachMCtrack(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map,std::vector<V0> *v0vec, Cascade &casc);
std::array<int,2> getV0MCtrack(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map,std::vector<V0> *v0vec, Cascade &casc);
double calcCascAlpha(const Cascade &cascade);
double calcMass(const Cascade &casc, int firstV0dauPDG, int secondV0dauPDG);



template <typename TrackT>
std::vector<TrackT>* fetchTracks(TFile* file, const char* treename, const char* branchname)
{
  auto tree = (TTree*)file->Get(treename);
  auto br = tree->GetBranch(branchname);
  std::vector<TrackT>* tracks = nullptr;
  br->SetAddress(&tracks);
  br->GetEntry(0);
  return tracks;
}



void CascTrackTreeBuilder(std::string path, std::string outSuffix = "")
{
    auto strUtils = strTrackingUtils(motherPDG, firstDaughterPDG, secondDaughterPDG);

    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    std::vector<TString> kine_files;

    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (file.substr(0, 2) == "tf" &&  std::stoi(file.substr(2, 3))<10)//201)
        {
            dirs.push_back(path + "/" + file);
            auto innerdir = (TSystemDirectory *)fileObj;
            auto innerfiles = innerdir->GetListOfFiles();
            // if (! innerfiles->Contains("o2simdigitizerworkflow_configuration.ini")){
            //     std::cout<<"Skipping "<<file<<": Not complete.\n";
            //     continue;
            // }
            // dirs.push_back(path + "/" + file);

            for (auto innerfileObj : *innerfiles)
            {
                TString innerfile = ((TSystemFile *)innerfileObj)->GetName();
                if (innerfile.EndsWith("Kine.root") && innerfile.Contains("sgn"))
                {
                    kine_files.push_back(innerfile);
                }
            }
        }
    }

    TFile outFile = TFile(Form("TrackedCascTree%s.root", outSuffix.data()), "recreate");
    TTree *outTree = new TTree("CascTree", "CascTree");
    float gPt, gCt, gR2, recoPt, recoR2, recoMass, ITStrackPt, trackedPt, trackedR2, trackedMass;
    float recoDCAXY, recoDCAZ, trackedDCAXY, trackedDCAZ;
    float bachDCAXY, bachDCAZ;
    float recoPvX, recoPvY, recoPvZ;
    bool isTrueVertex, isTrueCasc, isTrackedCasc, isTrueCascTrack, isCascTrackable, isDuplicated;
    float matchchi2;

    outTree->Branch("gPt", &gPt);
    outTree->Branch("gCt", &gCt);
    outTree->Branch("gR2", &gR2);
    outTree->Branch("isTrueCasc", &isTrueCasc);
    outTree->Branch("isTrueVertex", &isTrueVertex);
    outTree->Branch("recoPvX", &recoPvX);
    outTree->Branch("recoPvY", &recoPvY);
    outTree->Branch("recoPvZ", &recoPvZ);
    outTree->Branch("recoPt", &recoPt);
    outTree->Branch("recoR2", &recoR2);
    outTree->Branch("recoMass", &recoMass);
    outTree->Branch("ITStrackPt", &ITStrackPt);
    outTree->Branch("isCascTrackable", &isCascTrackable);
    outTree->Branch("isTrackedCasc", &isTrackedCasc);
    outTree->Branch("isTrueCascTrack", &isTrueCascTrack);
    outTree->Branch("isDuplicated", &isDuplicated);
    outTree->Branch("trackedPt", &trackedPt);
    outTree->Branch("trackedR2", &trackedR2);
    outTree->Branch("trackedMass", &trackedMass);
    outTree->Branch("recoDCAXY", &recoDCAXY);
    outTree->Branch("recoDCAZ", &recoDCAZ);
    outTree->Branch("trackedDCAXY", &trackedDCAXY);
    outTree->Branch("trackedDCAZ", &trackedDCAZ);

    //create daughters tree
    TTree *dauTree=new TTree("DauTree", "DauTree");
    float bach_gPt, v0_gPt, bachelorPt, v0_recoPt, v0_trackedPt;
    float v0_recoVX, v0_recoVY, v0_recoVZ;
    float v0decay_recoVX, v0decay_recoVY, v0decay_recoVZ;
    dauTree->Branch("bach_gPt", &bach_gPt);
    dauTree->Branch("v0_gPt", &v0_gPt);
    dauTree->Branch("bachelorPt", &bachelorPt);
    dauTree->Branch("v0_recoPt", &v0_recoPt);
    dauTree->Branch("v0_trackedPt", &v0_trackedPt);
    dauTree->Branch("v0_recoVX", &v0_recoVX);
    dauTree->Branch("v0_recoVY", &v0_recoVY);
    dauTree->Branch("v0_recoVZ", &v0_recoVZ);
    dauTree->Branch("v0decay_recoVX", &v0decay_recoVX);
    dauTree->Branch("v0decay_recoVY", &v0decay_recoVY);
    dauTree->Branch("v0decay_recoVZ", &v0decay_recoVZ);
    dauTree->Branch("bachDCAZ", &bachDCAZ);
    dauTree->Branch("bachDCAXY", &bachDCAXY);
    dauTree->Branch("matchchi2", &matchchi2);

    // create MC tree for efficiency calculation
    TTree *mcTree = new TTree("MCTree", "MCTree");
    float mcPt, mcCt, mcR2;
    mcTree->Branch("mcPt", &mcPt);
    mcTree->Branch("mcCt", &mcCt);
    mcTree->Branch("mcR2", &mcR2);

    // Geometry
    o2::base::GeometryManager::loadGeometry(dirs[0] + "/o2sim_geometry.root");
    auto gman = o2::its::GeometryTGeo::Instance();
    gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));

    const auto grp = o2::parameters::GRPObject::loadFrom(dirs[0] + "/" + "o2sim_grp.root");
    o2::base::Propagator::initFieldFromGRP(grp);
    auto propagator = o2::base::Propagator::Instance();

    int counter = 0;
    for (unsigned int i = 0; i < dirs.size(); i++)
    {
        auto &dir = dirs[i];
        auto &kine_file = kine_files[i];
        LOG(info) << "Processing " << dir;
        // Files
        auto fMCTracks = TFile::Open((TString(dir + "/") + kine_file));
        auto fStrangeTracks = TFile::Open((dir + "/o2_strange_tracks.root").data());
        auto fSecondaries = TFile::Open((dir + "/o2_secondary_vertex.root").data());
        auto fPrimaryVertex = TFile::Open((dir + "/o2_primary_vertex.root").data());

        auto fITSTPC = TFile::Open((dir + "/o2match_itstpc.root").data());
        auto fTPCTOF = TFile::Open((dir + "/o2match_tof_tpc.root").data());
        auto fTPCTRD = TFile::Open((dir + "/trdmatches_tpc.root").data());
        auto fITSTPCTOF = TFile::Open((dir + "/o2match_tof_itstpc.root").data());
        auto fITS = TFile::Open((dir + "/o2trac_its.root").data());
        auto fClusITS = TFile::Open((dir + "/o2clus_its.root").data());
        auto fTPC = TFile::Open((dir + "/tpctracks.root").data());
        auto fITSTPCTRD = TFile::Open((dir + "/trdmatches_itstpc.root").data());
        auto fTPCTRDTOF = TFile::Open((dir + "/o2match_tof_tpctrd.root").data());
        auto fITSTPCTRDTOF = TFile::Open((dir + "/o2match_tof_itstpctrd.root").data());

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treePrimaryVertex = (TTree *)fPrimaryVertex->Get("o2sim");
        auto treeStrangeTracks = (TTree *)fStrangeTracks->Get("o2sim");
        auto treeSecondaries = (TTree *)fSecondaries->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
        auto treeITSTPCTOF = (TTree *)fITSTPCTOF->Get("matchTOF");
        auto treeTPCTOF = (TTree *)fTPCTOF->Get("matchTOF");
        auto treeTPCTRD = (TTree *)fTPCTRD->Get("tracksTRD");
        auto treeITSTPCTRD = (TTree *)fITSTPCTRD->Get("tracksTRD");
        auto treeTPCTRDTOF = (TTree *)fTPCTRDTOF->Get("matchTOF");
        auto treeITSTPCTRDTOF = (TTree *)fITSTPCTRDTOF->Get("matchTOF");

        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSclus = (TTree *)fClusITS->Get("o2sim");
        auto treeTPC = (TTree *)fTPC->Get("tpcrec");

        // MC Tracks
        std::vector<o2::MCTrack> *MCtracks = nullptr;
        std::vector<o2::itsmft::Hit> *ITSHits = nullptr;

        // Primary Vertex
        std::vector<o2::dataformats::PrimaryVertex> *primVertices = nullptr;

        // Hypertracks
        std::vector<StrangeTrack> *strangeTrackVec = nullptr;
        std::vector<o2::track::TrackParametrizationWithError<float>> *hypertrackVec = nullptr;
        std::vector<o2::strangeness_tracking::ClusAttachments> *nAttachments = nullptr;

        // Secondary Vertices
        std::vector<Cascade> *cascVec = nullptr;
        std::vector<V0> *v0Vec = nullptr;

        // ITS tracks
        std::vector<o2::its::TrackITS> *ITStracks = nullptr;

        // other tracks
        std::vector<o2::tpc::TrackTPC> *TPCtracks = nullptr;
        std::vector<o2::dataformats::TrackTPCITS> *TPCITStracks = nullptr;

        // Labels
        std::vector<o2::MCEventLabel> *pvMcArr = nullptr;
        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTRDvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTRDvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTRDTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTRDTOFvec = nullptr;

        // Clusters
        std::vector<CompClusterExt> *ITSclus = nullptr;
        o2::dataformats::MCTruthContainer<o2::MCCompLabel> *clusLabArr = nullptr;
        std::vector<int> *ITSTrackClusIdx = nullptr;
        std::vector<unsigned char> *ITSpatt = nullptr;

        // Setting branches
        treeStrangeTracks->SetBranchAddress("StrangeTracks", &strangeTrackVec);
        treeStrangeTracks->SetBranchAddress("ClusUpdates", &nAttachments);
        treeSecondaries->SetBranchAddress("Cascades", &cascVec);
        treeSecondaries->SetBranchAddress("V0s", &v0Vec);
        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
        treePrimaryVertex->SetBranchAddress("PrimaryVertex", &primVertices);
        treePrimaryVertex->SetBranchAddress("PVMCTruth", &pvMcArr);

        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeITS->SetBranchAddress("ITSTrack", &ITStracks);
        treeTPC->SetBranchAddress("TPCTracks", &TPCtracks);
        treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
        treeITSTPC->SetBranchAddress("TPCITS", &TPCITStracks);
        treeTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTOFvec);
        treeTPCTRD->SetBranchAddress("labels", &labTPCTRDvec);
        treeITSTPCTRD->SetBranchAddress("labelsTRD", &labITSTPCTRDvec);
        treeTPCTRDTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTRDTOFvec);
        treeITSTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTOFvec);
        treeITSTPCTRDTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTRDTOFvec);
        treeITS->SetBranchAddress("ITSTrackClusIdx", &ITSTrackClusIdx);
        treeITSclus->SetBranchAddress("ITSClusterComp", &ITSclus);
        treeITSclus->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);

        // define detector map
        std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}, {"TPC", labTPCvec}, {"ITS-TPC", labITSTPCvec}, {"TPC-TOF", labTPCTOFvec}, {"TPC-TRD", labTPCTRDvec}, {"ITS-TPC-TOF", labITSTPCTOFvec}, {"ITS-TPC-TRD", labITSTPCTRDvec}, {"TPC-TRD-TOF", labTPCTRDTOFvec}, {"ITS-TPC-TRD-TOF", labITSTPCTRDTOFvec}};
        std::map<std::string, std::vector<o2::MCCompLabel> *> trackMap{{"ITS", labITSvec}, {"TPC", labTPCvec}, {"ITS-TPC", labITSTPCvec}, {"TPC-TOF", labTPCTOFvec}, {"TPC-TRD", labTPCTRDvec}, {"ITS-TPC-TOF", labITSTPCTOFvec}, {"ITS-TPC-TRD", labITSTPCTRDvec}, {"TPC-TRD-TOF", labTPCTRDTOFvec}, {"ITS-TPC-TRD-TOF", labITSTPCTRDTOFvec}};

        // fill MC matrix
        std::vector<std::vector<o2::MCTrack>> mcTracksMatrix;

        auto nev = treeMCTracks->GetEntriesFast();
        mcTracksMatrix.resize(nev);
        for (int n = 0; n < nev; n++)
        { // loop over MC events
            treeMCTracks->GetEvent(n);
            mcTracksMatrix[n].resize(MCtracks->size());
            for (unsigned int mcI{0}; mcI < MCtracks->size(); ++mcI)
            {
                mcTracksMatrix[n][mcI] = MCtracks->at(mcI);
                if (abs(mcTracksMatrix[n][mcI].GetPdgCode()) == motherPDG)
                {
                    mcPt = mcTracksMatrix[n][mcI].GetPt();
                    mcR2 = strUtils.calcR2(MCtracks, mcTracksMatrix[n][mcI]);
                    mcCt = strUtils.calcLifetime(MCtracks, mcTracksMatrix[n][mcI]);
                    mcTree->Fill();
                }
            }
        }

        // Starting matching Cascades and ITS tracks
        int counterV0 = 0;
        for (int frame = 0; frame < treeITS->GetEntriesFast(); frame++)
        {

            if (!treeITS->GetEvent(frame) || !treeSecondaries->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeITSclus->GetEvent(frame) || !treeTPCTRD->GetEvent(frame) || !treeStrangeTracks->GetEvent(frame) ||
                !treeITSTPCTRD->GetEvent(frame) || !treeTPCTRDTOF->GetEvent(frame) || !treeITSTPCTRDTOF->GetEvent(frame) || !treePrimaryVertex->GetEvent(frame))
                {

                continue;
                }

            for (unsigned int iCascVec = 0; iCascVec < cascVec->size(); iCascVec++)
            {
                std::cout << "Cascades: " << iCascVec << "/" << cascVec->size() <<  std::endl;
                bool isTreeFilled = false;

                // default tree values
                gPt = -1, gCt = -1, gR2 = -1, recoPt = -1, recoR2 = -1, recoMass = -1, ITStrackPt = -1, trackedPt = -1, trackedR2 = -1, trackedMass = -1;
                recoPvX = -1, recoPvY = -1, recoPvZ = -1;
                recoDCAXY = -1, recoDCAZ = -1, trackedDCAXY = -1, trackedDCAZ = -1;
                bachDCAXY = -1; bachDCAZ = -1;
                isTrueVertex = false, isTrueCasc = false, isTrackedCasc = false, isTrueCascTrack = false, isCascTrackable = false, isDuplicated = false;
                bach_gPt = -1, v0_gPt = -1, bachelorPt = -1, v0_recoPt = -1, v0_trackedPt = -1;
                v0_recoVX = -1, v0_recoVY = -1, v0_recoVZ = -1;
                v0decay_recoVX = -1, v0decay_recoVY = -1, v0decay_recoVZ = -1;
                matchchi2 = -1;
                

                auto &casc = cascVec->at(iCascVec);
                bool isMatter = calcCascAlpha(casc) > 0;
                auto cascMCref = matchCascToMC(mcTracksMatrix, map, v0Vec, casc);
                std::array<int,2> bachMCref = getBachMCtrack(mcTracksMatrix, map,v0Vec, casc);
                std::array<int,2> v0MCref = getV0MCtrack(mcTracksMatrix, map,v0Vec, casc);

                recoPt = casc.getPt();
                recoR2 = sqrt(casc.calcR2());
                recoMass = isMatter ? calcMass(casc, firstDaughterPDG, secondDaughterPDG) : calcMass(casc, secondDaughterPDG, firstDaughterPDG);

                auto &primaryVertex = primVertices->at(casc.getVertexID());
                auto &pvLabel = pvMcArr->at(casc.getVertexID());
                int itsIdx = -1; // loop over ITS track for signal V0s

                auto &secondaryVertex = v0Vec->at(casc.getV0ID());
                //omega decay vertex
                v0_recoVX =casc.getX();
                v0_recoVY = casc.getY();
                v0_recoVZ = casc.getZ();

                //lambda decay vertex
                v0decay_recoVX = secondaryVertex.getX();
                v0decay_recoVY = secondaryVertex.getY();
                v0decay_recoVZ = secondaryVertex.getZ();

                // LOG(info)<<"************************** index"<<casc.getBachelorID().getIndex();
                
                // o2::its:t

                o2::track::TrackParF bachelorTrack;

                switch (casc.getBachelorID().getSource()) {
                    // case o2::dataformats::V0::GIndex::ITS:
                    case o2::dataformats::V0::GIndex::ITS:
                    // case o2::dataformats::VtxTrackIndex::Source::ITS:
                        bachelorPt = ITStracks->at(casc.getBachelorID().getIndex()).getPt();
                        LOG(info)<<"********* ITS " << ITStracks->size() << " " << bachelorPt;
                        // auto ITSbachelor =  ITStracks->at(casc.getBachelorID().getIndex());
                        bachelorTrack = ITStracks->at(casc.getBachelorID().getIndex());
                        break;
                    case o2::dataformats::V0::GIndex::TPC:
                        bachelorPt = TPCtracks->at(casc.getBachelorID().getIndex()).getPt();
                        LOG(info)<<"********* TPC "<< TPCtracks->size()<<"   " <<bachelorPt;
                        // auto TPCbachelor = TPCtracks->at(casc.getBachelorID().getIndex());
                        bachelorTrack= TPCtracks->at(casc.getBachelorID().getIndex());
                        break;
                    case o2::dataformats::V0::GIndex::ITSTPC:
                        auto trackTPCITS = TPCITStracks->at(casc.getBachelorID().getIndex());
                        matchchi2= trackTPCITS.getChi2Match();
                        auto source = trackTPCITS.getRefITS().getSource();
                        if (source !=24 && matchchi2>1 ){
                            bachelorPt = trackTPCITS.getPt();
                            LOG(info)<<"***** ITS TPC not 24 "<<TPCITStracks->size()<<"   "<< bachelorPt;
                            // auto bachelor = trackTPCITS;
                            bachelorTrack = trackTPCITS;
                            break;
                        }
                        else
                            LOG(info)<<"********* TPC ITS24: skipping ";
                            break;
                }  


                v0_recoPt = secondaryVertex.getPt();

                if (v0MCref[0]!= -1 && v0MCref[1] != -1)
                {
                    auto &mcV0 = mcTracksMatrix[v0MCref[0]][v0MCref[1]];
                    v0_gPt= mcV0.GetPt();
                }

                if (bachMCref[0]!= -1 && bachMCref[1] != -1)
                {
                    auto &mcBach = mcTracksMatrix[bachMCref[0]][bachMCref[1]];
                    bach_gPt= mcBach.GetPt();
                }


                if (cascMCref[0] != -1 && cascMCref[1] != -1)
                {
                    counterV0++;

                    isTrueCasc = true;
                    auto &mcCasc = mcTracksMatrix[cascMCref[0]][cascMCref[1]];
                    if (pvLabel.getEventID() == cascMCref[0])
                        isTrueVertex = true;

                    recoPvX = primaryVertex.getX();
                    recoPvY = primaryVertex.getY();
                    recoPvZ = primaryVertex.getZ();

                    gPt = mcCasc.GetPt();
                    gR2 = strUtils.calcR2(&mcTracksMatrix[cascMCref[0]], mcCasc);
                    gCt = strUtils.calcLifetime(&mcTracksMatrix[cascMCref[0]], mcCasc);

                    gpu::gpustd::array<float, 2> dca{-999.f, -999.f};
                    if (propagator->propagateToDCA(primaryVertex.getXYZ(), casc, propagator->getNominalBz(), 2.f, o2::base::PropagatorF::MatCorrType::USEMatCorrLUT, &dca))
                    {
                        recoDCAXY = dca[0];
                        recoDCAZ = dca[1];
                    }


                    // gpu::gpustd::array<float, 2> dcaBach{-999.f, -999.f};
                    // if (propagator->propagateToDCA(primaryVertex.getXYZ(), bachelorTrack, propagator->getNominalBz(), 2.f, o2::base::PropagatorF::MatCorrType::USEMatCorrLUT, &dcaBach))
                    // {
                    //     bachDCAXY = dcaBach[0];
                    //     bachDCAZ = dcaBach[1];
                    // }

                    // Matching ITS tracks to MC tracks and V0
                    std::array<int, 2> ITSref = {-1, 1};
                    o2::its::TrackITS ITStrack;

                    for (unsigned int iITStrack = 0; iITStrack < ITStracks->size(); iITStrack++)
                    {
                        auto &labITS = (*labITSvec)[iITStrack];
                        auto &trackIdx = (*ITSTrackClusIdx)[iITStrack];

                        ITSref = strUtils.matchITStracktoMC(mcTracksMatrix, labITS);
                        ITStrack = (*ITStracks)[iITStrack];

                        if (ITSref[0] == cascMCref[0] && ITSref[1] == cascMCref[1])
                        {
                            // LOG(info) << "++++++++++++++++";
                            // LOG(info) << "Casc + ITS track found! ";
                            // LOG(info) << "Casc Momentum: " << casc.getPt() << ", ITS track momentum: " << ITStrack.getPt();
                            isCascTrackable = true;
                            ITStrackPt = ITStrack.getPt();
                            itsIdx = iITStrack;
                            break;
                        }
                    }
                }

                // check whether the corresponding strange track is found: both signal and background checked

                int cloneCounter = 0;
                for (unsigned int iStTr = 0; iStTr < strangeTrackVec->size(); iStTr++)
                {
                    auto &strangeTrack = strangeTrackVec->at(iStTr);
                    auto &itsRef = strangeTrack.mITSRef;
                    auto &cascRef = strangeTrack.mDecayRef;
                    if (!(strangeTrack.mPartType == o2::dataformats::kStrkCascade) || cascRef != iCascVec)
                    {
                        continue;
                    }

                    auto &sTrack = strangeTrack.mMother;
                    isTrackedCasc = true;
                    isDuplicated = cloneCounter > 0;
                    trackedPt = sqrt(strangeTrack.mDecayMom[0] * strangeTrack.mDecayMom[0] + strangeTrack.mDecayMom[1] * strangeTrack.mDecayMom[1]);
                    trackedR2 = strangeTrack.mDecayVtx[0] * strangeTrack.mDecayVtx[0] + strangeTrack.mDecayVtx[1] * strangeTrack.mDecayVtx[1];

                    ///////////?????
                    v0_trackedPt = sqrt(strangeTrack.mDecayVtx[0] * strangeTrack.mDecayVtx[0] + strangeTrack.mDecayVtx[1] * strangeTrack.mDecayVtx[1]);
                    if (int(itsRef) == itsIdx && cascRef == iCascVec)
                    {
                        // LOG(info) << "Strange track found! ";
                        isTrueCascTrack = true;
                        gpu::gpustd::array<float, 2> dcaT{-999.f, -999.f};
                        if (propagator->propagateToDCA(primaryVertex.getXYZ(), sTrack, propagator->getNominalBz(), 2.f, o2::base::PropagatorF::MatCorrType::USEMatCorrLUT, &dcaT))
                        {
                            trackedDCAXY = dcaT[0];
                            trackedDCAZ = dcaT[1];
                        }
                    }
                    else
                    {
                        isTrueCascTrack = false;
                    }
                    outTree->Fill();
                    dauTree->Fill();
                    isTreeFilled = true; // double filling if two strange tracks are associated to the same V0
                    cloneCounter++;
                }

                if (!isTreeFilled)
                {
                    outTree->Fill();
                    dauTree->Fill();
                }
            }
        }
    }

    outFile.cd();
    outTree->Write();
    mcTree->Write();
    dauTree->Write();
}

std::array<int, 2> matchCascToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec, Cascade &casc)
{
    std::array<int, 2> motherVec{-1, -1};
    std::array<std::array<int, 2>, 2> v0DauRefs;
    std::array<int, 2> bachRef;

    auto v0Idx = casc.getV0ID();
    auto &v0 = v0vec->at(v0Idx);

    auto bachID = casc.getBachelorID();
    if (!map[bachID.getSourceName()])
        return motherVec;
    auto &bachLab = map[bachID.getSourceName()]->at(bachID.getIndex());
    if (bachLab.isValid())
        bachRef = {bachLab.getEventID(), bachLab.getTrackID()};
    else
        return motherVec;

    for (unsigned int iV0 = 0; iV0 < 2; iV0++)
    {
        v0DauRefs[iV0] = {-1, -1};
        if (map[v0.getProngID(iV0).getSourceName()])
        {
            auto labTrackType = map[v0.getProngID(iV0).getSourceName()];
            auto lab = labTrackType->at(v0.getProngID(iV0).getIndex());

            int trackID, evID, srcID;
            bool fake;
            lab.get(trackID, evID, srcID, fake);
            if (lab.isValid())
            {
                v0DauRefs[iV0] = {lab.getEventID(), lab.getTrackID()};
            }
        }
        else
        {
            LOG(info) << "No MC label for " << v0.getProngID(iV0).getSourceName();
        }
    }

    if (v0DauRefs[0][1] == -1 || v0DauRefs[1][1] == -1)
        return motherVec;

    auto &dau1MC = mcTracksMatrix[v0DauRefs[0][0]][v0DauRefs[0][1]];
    auto &dau2MC = mcTracksMatrix[v0DauRefs[1][0]][v0DauRefs[1][1]];

    if (!(std::abs(dau1MC.GetPdgCode()) == firstV0dauPDG && std::abs(dau2MC.GetPdgCode()) == secondV0dauPDG) && !(std::abs(dau1MC.GetPdgCode()) == secondV0dauPDG && std::abs(dau2MC.GetPdgCode()) == firstV0dauPDG))
        return motherVec;

    if (!dau1MC.isSecondary() || !dau2MC.isSecondary() || dau1MC.getMotherTrackId() != dau2MC.getMotherTrackId())
    {
        // LOG(info) << dau1MC.GetPdgCode() << " " << dau2MC.GetPdgCode() << " " << dau1MC.isSecondary() << " " << dau2MC.isSecondary() << " " << dau1MC.getMotherTrackId() << " " << dau2MC.getMotherTrackId();
        return motherVec;
    }

    auto v0MC = mcTracksMatrix[v0DauRefs[0][0]][dau1MC.getMotherTrackId()];
    auto &bachMC = mcTracksMatrix[bachRef[0]][bachRef[1]];
    // LOG(info) << "++++++++++ bach"<<bachMC.getMotherTrackId();


    if (std::abs(v0MC.GetPdgCode()) != firstDaughterPDG || !v0MC.isSecondary() || !bachMC.isSecondary())
        return motherVec;

    auto cascMC = mcTracksMatrix[v0DauRefs[0][0]][v0MC.getMotherTrackId()];
    if (v0MC.getMotherTrackId() != bachMC.getMotherTrackId() || std::abs(cascMC.GetPdgCode()) != motherPDG)
        return motherVec;

    motherVec = {v0DauRefs[0][0], v0MC.getMotherTrackId()};
    return motherVec;
}

double calcMass(const Cascade &casc, int firstDauPdg, int secondDauPdg)
{
    double dauMass[2];
    int dauCharges[2];
    auto firstV0dau = TDatabasePDG::Instance()->GetParticle(firstDauPdg);
    auto secondV0dau = TDatabasePDG::Instance()->GetParticle(secondDauPdg);
    if (firstDauPdg == 1000020030 || secondDauPdg == 1000020030)
    {
        int indexHe3 = firstDauPdg == 1000020030 ? 0 : 1;
        dauMass[indexHe3] = 2.808391;
        dauCharges[indexHe3] = 2;
        dauMass[1 - indexHe3] = firstDauPdg == 1000020030 ? secondV0dau->Mass() : firstV0dau->Mass();
        dauCharges[1 - indexHe3] = firstDauPdg == 1000020030 ? secondV0dau->Charge() / 3 : firstV0dau->Charge() / 3;
    }

    else
    {
        dauMass[0] = firstV0dau->Mass();
        dauCharges[0] = firstV0dau->Charge();
        dauMass[1] = secondV0dau->Mass();
        dauCharges[1] = secondV0dau->Charge();
    }

    std::vector<o2::dataformats::Cascade::Track> dauTracks = {casc.getProng(0), casc.getProng(1)};
    TLorentzVector moth, prong;
    std::array<float, 3> p;
    for (int i = 0; i < 2; i++)
    {
        auto &track = dauTracks[i];
        auto &mass = dauMass[i];
        track.getPxPyPzGlo(p);
        int charge = dauCharges[i];
        prong.SetVectM({charge * p[0], charge * p[1], charge * p[2]}, mass);
        moth += prong;
    }
    return moth.M();
}

double calcCascAlpha(const Cascade &casc)
{
    std::array<float, 3> fV0mom, fPmom, fNmom = {0, 0, 0};
    casc.getProng(0).getPxPyPzGlo(fPmom);
    casc.getProng(1).getPxPyPzGlo(fNmom);
    casc.getPxPyPzGlo(fV0mom);

    TVector3 momNeg(fNmom[0], fNmom[1], fNmom[2]);
    TVector3 momPos(fPmom[0], fPmom[1], fPmom[2]);
    TVector3 momTot(fV0mom[0], fV0mom[1], fV0mom[2]);

    Double_t lQlNeg = momNeg.Dot(momTot) / momTot.Mag();
    Double_t lQlPos = momPos.Dot(momTot) / momTot.Mag();

    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

std::array<int,2> getBachMCtrack(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map,std::vector<V0> *v0vec, Cascade &casc)
{
    std::array<int, 2> null_vec{-1, -1};
    std::array<int, 2> bachRef;

    auto v0Idx = casc.getV0ID();
    auto &v0 = v0vec->at(v0Idx);

    auto bachID = casc.getBachelorID();
    if (!map[bachID.getSourceName()])
        return null_vec;
    auto &bachLab = map[bachID.getSourceName()]->at(bachID.getIndex());
    if (bachLab.isValid())
        bachRef = {bachLab.getEventID(), bachLab.getTrackID()};
    else
        return null_vec;

    auto &bachMC = mcTracksMatrix[bachRef[0]][bachRef[1]];

    if (std::abs(bachMC.GetPdgCode())!= secondDaughterPDG)
        return null_vec;

    // LOG(info) << "++++++++++ bach"<<std::abs(bachMC.GetPdgCode());

    return bachRef;}



std::array<int, 2> getV0MCtrack(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec, Cascade &casc)
{
    std::array<int, 2> null_vec{-1, -1};
    std::array<std::array<int, 2>, 2> v0DauRefs;
    std::array<int, 2> bachRef;

    auto v0Idx = casc.getV0ID();
    auto &v0 = v0vec->at(v0Idx);

    for (unsigned int iV0 = 0; iV0 < 2; iV0++)
    {
        v0DauRefs[iV0] = {-1, -1};
        if (map[v0.getProngID(iV0).getSourceName()])
        {
            auto labTrackType = map[v0.getProngID(iV0).getSourceName()];
            auto lab = labTrackType->at(v0.getProngID(iV0).getIndex());

            int trackID, evID, srcID;
            bool fake;
            lab.get(trackID, evID, srcID, fake);
            if (lab.isValid())
            {
                v0DauRefs[iV0] = {lab.getEventID(), lab.getTrackID()};
            }
        }
        else
        {
            LOG(info) << "No MC label for " << v0.getProngID(iV0).getSourceName();
        }
    }

    if (v0DauRefs[0][1] == -1 || v0DauRefs[1][1] == -1)
        return null_vec;

    auto &dau1MC = mcTracksMatrix[v0DauRefs[0][0]][v0DauRefs[0][1]];
    auto &dau2MC = mcTracksMatrix[v0DauRefs[1][0]][v0DauRefs[1][1]];

    if (!(std::abs(dau1MC.GetPdgCode()) == firstV0dauPDG && std::abs(dau2MC.GetPdgCode()) == secondV0dauPDG) && !(std::abs(dau1MC.GetPdgCode()) == secondV0dauPDG && std::abs(dau2MC.GetPdgCode()) == firstV0dauPDG))
        return null_vec;

    if (!dau1MC.isSecondary() || !dau2MC.isSecondary() || dau1MC.getMotherTrackId() != dau2MC.getMotherTrackId())
    {
        // LOG(info) << dau1MC.GetPdgCode() << " " << dau2MC.GetPdgCode() << " " << dau1MC.isSecondary() << " " << dau2MC.isSecondary() << " " << dau1MC.getMotherTrackId() << " " << dau2MC.getMotherTrackId();
        return null_vec;
    }

    auto &v0MC = mcTracksMatrix[v0DauRefs[0][0]][dau1MC.getMotherTrackId()];
    std::array<int,2> v0MCRef = {v0DauRefs[0][0], dau1MC.getMotherTrackId()};

    if (std::abs(v0MC.GetPdgCode()) != firstDaughterPDG)
        return null_vec;

    // LOG(info) << "++++++++++ VO"<<std::abs(v0MC.GetPdgCode());

    return v0MCRef;
}
