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
#include "DetectorsVertexing/DCAFitterN.h"
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
using StrangeTrack = o2::strangeness_tracking::StrangeTrack;

const int motherPDG = 1010010030;
const int firstDaughterPDG = 1000020030;
const int secondDaughterPDG = 211;

std::array<int, 2> matchV0ToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, V0 &v0);
double calcV0alpha(const V0 &v0);
double calcMass(const V0 &v0, int firstV0dauPDG, int secondV0dauPDG);

void V0TrackTreeBuilder(std::string path, std::string outSuffix = "")
{
    auto strUtils = strTrackingUtils(motherPDG, firstDaughterPDG, secondDaughterPDG);

    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    std::vector<TString> kine_files;

    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (file.substr(0, 2) == "tf")
        {
            if (stoi(file.substr(2)) > 1)
                continue;
            dirs.push_back(path + "/" + file);
            auto innerdir = (TSystemDirectory *)fileObj;
            auto innerfiles = innerdir->GetListOfFiles();
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

    TFile outFile = TFile(Form("TrackedV0Tree%s.root", outSuffix.data()), "recreate");
    TTree *outTree = new TTree("V0Tree", "V0Tree");
    float gPt, gCt, gR2, gX, gY, gZ, recoPt, recoR2, recoMass, ITStrackPt, trackedPt, trackedR2, trackedMass;
    float recoDCAXY, recoDCAZ, trackedDCAXY, trackedDCAZ;
    float recoPvX, recoPvY, recoPvZ;
    bool isTrueVertex, isTrueV0, isTrackedV0, isTrueV0Track, isV0Trackable, isDuplicated;

    outTree->Branch("gX", &gX);
    outTree->Branch("gY", &gY);
    outTree->Branch("gZ", &gZ);
    outTree->Branch("gPt", &gPt);
    outTree->Branch("gCt", &gCt);
    outTree->Branch("gR2", &gR2);
    outTree->Branch("isTrueV0", &isTrueV0);
    outTree->Branch("isTrueVertex", &isTrueVertex);
    outTree->Branch("recoPvX", &recoPvX);
    outTree->Branch("recoPvY", &recoPvY);
    outTree->Branch("recoPvZ", &recoPvZ);
    outTree->Branch("recoPt", &recoPt);
    outTree->Branch("recoR2", &recoR2);
    outTree->Branch("recoMass", &recoMass);
    outTree->Branch("ITStrackPt", &ITStrackPt);
    outTree->Branch("isV0Trackable", &isV0Trackable);
    outTree->Branch("isTrackedV0", &isTrackedV0);
    outTree->Branch("isTrueV0Track", &isTrueV0Track);
    outTree->Branch("isDuplicated", &isDuplicated);
    outTree->Branch("trackedPt", &trackedPt);
    outTree->Branch("trackedR2", &trackedR2);
    outTree->Branch("trackedMass", &trackedMass);
    outTree->Branch("recoDCAXY", &recoDCAXY);
    outTree->Branch("recoDCAZ", &recoDCAZ);
    outTree->Branch("trackedDCAXY", &trackedDCAXY);
    outTree->Branch("trackedDCAZ", &trackedDCAZ);

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
        treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
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


        // Starting matching V0s and ITS tracks
        int counterV0 = 0;
        for (int frame = 0; frame < treeITS->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeSecondaries->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeITSclus->GetEvent(frame) || !treeTPCTRD->GetEvent(frame) || !treeStrangeTracks->GetEvent(frame) ||
                !treeITSTPCTRD->GetEvent(frame) || !treeTPCTRDTOF->GetEvent(frame) || !treeITSTPCTRDTOF->GetEvent(frame) || !treePrimaryVertex->GetEvent(frame))
                continue;

            for (unsigned int iV0Vec = 0; iV0Vec < v0Vec->size(); iV0Vec++)
            {
                bool isTreeFilled = false;

                // default tree values
                gPt = -1, gCt = -1, gR2 = -1, recoPt = -1, recoR2 = -1, recoMass = -1, ITStrackPt = -1, trackedPt = -1, trackedR2 = -1, trackedMass = -1;
                recoPvX = -1, recoPvY = -1, recoPvZ = -1;
                recoDCAXY = -1, recoDCAZ = -1, trackedDCAXY = -1, trackedDCAZ = -1;
                isTrueVertex = false, isTrueV0 = false, isTrackedV0 = false, isTrueV0Track = false, isV0Trackable = false, isDuplicated = false;

                auto &v0 = v0Vec->at(iV0Vec);
                bool isMatter = calcV0alpha(v0) > 0;
                auto v0MCref = matchV0ToMC(mcTracksMatrix, map, v0);

                recoPt = v0.getPt();
                recoR2 = sqrt(v0.calcR2());
                recoMass = isMatter ? calcMass(v0, firstDaughterPDG, secondDaughterPDG) : calcMass(v0, secondDaughterPDG, firstDaughterPDG);

                auto &primaryVertex = primVertices->at(v0.getVertexID());
                auto &pvLabel = pvMcArr->at(v0.getVertexID());

                int itsIdx = -1; // loop over ITS track for signal V0s

                if (v0MCref[0] != -1 && v0MCref[1] != -1)
                {
                    counterV0++;

                    isTrueV0 = true;
                    auto &mcV0 = mcTracksMatrix[v0MCref[0]][v0MCref[1]];
                    if (pvLabel.getEventID() == v0MCref[0])
                        isTrueVertex = true;

                    recoPvX = primaryVertex.getX();
                    recoPvY = primaryVertex.getY();
                    recoPvZ = primaryVertex.getZ();

                    gPt = mcV0.GetPt();
                    gR2 = strUtils.calcR2(&mcTracksMatrix[v0MCref[0]], mcV0);
                    gCt = strUtils.calcLifetime(&mcTracksMatrix[v0MCref[0]], mcV0);

                    gpu::gpustd::array<float, 2> dca{-999.f, -999.f};
                    if (propagator->propagateToDCA(primaryVertex.getXYZ(), v0, propagator->getNominalBz(), 2.f, o2::base::PropagatorF::MatCorrType::USEMatCorrLUT, &dca))
                    {
                        recoDCAXY = dca[0];
                        recoDCAZ = dca[1];
                    }

                    // Matching ITS tracks to MC tracks and V0
                    std::array<int, 2> ITSref = {-1, 1};
                    o2::its::TrackITS ITStrack;
                    std::array<std::array<int, 2>, 7> clsRef;

                    for (unsigned int iITStrack = 0; iITStrack < ITStracks->size(); iITStrack++)
                    {
                        auto &labITS = (*labITSvec)[iITStrack];
                        auto &trackIdx = (*ITSTrackClusIdx)[iITStrack];

                        ITSref = strUtils.matchITStracktoMC(mcTracksMatrix, labITS);
                        ITStrack = (*ITStracks)[iITStrack];

                        if (ITSref[0] == v0MCref[0] && ITSref[1] == v0MCref[1])
                        {
                            LOG(info) << "++++++++++++++++";
                            LOG(info) << "V0 + ITS track found! ";
                            LOG(info) << "V0 Momentum: " << v0.getPt() << ", ITS track momentum: " << ITStrack.getPt();
                            isV0Trackable = true;
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
                    auto &v0Ref = strangeTrack.mDecayRef;
                    if (!(strangeTrack.mPartType == o2::strangeness_tracking::kV0) || v0Ref != int(iV0Vec))
                    {
                        continue;
                    }

                    auto &sTrack = strangeTrack.mMother;
                    isTrackedV0 = true;
                    isDuplicated = cloneCounter > 0;
                    trackedPt = sqrt(strangeTrack.decayMom[0] * strangeTrack.decayMom[0] + strangeTrack.decayMom[1] * strangeTrack.decayMom[1]);
                    trackedR2 = strangeTrack.decayVtx[0] * strangeTrack.decayVtx[0] + strangeTrack.decayVtx[1] * strangeTrack.decayVtx[1];

                    if (itsRef == int(itsIdx) && v0Ref == int(iV0Vec))
                    {
                        LOG(info) << "Strange track found! ";
                        isTrueV0Track = true;
                        gpu::gpustd::array<float, 2> dcaT{-999.f, -999.f};
                        if (propagator->propagateToDCA(primaryVertex.getXYZ(), sTrack, propagator->getNominalBz(), 2.f, o2::base::PropagatorF::MatCorrType::USEMatCorrLUT, &dcaT))
                        {
                            trackedDCAXY = dcaT[0];
                            trackedDCAZ = dcaT[1];
                        }
                    }
                    else
                    {
                        isTrueV0Track = false;
                    }
                    outTree->Fill();
                    isTreeFilled = true; // double filling if two strange tracks are associated to the same V0
                    cloneCounter++;
                }

                if (!isTreeFilled)
                {
                    outTree->Fill();
                }
            }
        }
    }

    outFile.cd();
    outTree->Write();
    mcTree->Write();
}

std::array<int, 2> matchV0ToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, V0 &v0)
{
    std::array<int, 2> motherVec{-1, -1};
    std::array<std::array<int, 2>, 2> v0DauRefs;

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
    }

    if (v0DauRefs[0][1] == -1 || v0DauRefs[1][1] == -1)
        return motherVec;

    auto &dau1MC = mcTracksMatrix[v0DauRefs[0][0]][v0DauRefs[0][1]];
    auto &dau2MC = mcTracksMatrix[v0DauRefs[1][0]][v0DauRefs[1][1]];

    if (!(std::abs(dau1MC.GetPdgCode()) == firstDaughterPDG && std::abs(dau2MC.GetPdgCode()) == secondDaughterPDG) && !(std::abs(dau1MC.GetPdgCode()) == secondDaughterPDG && std::abs(dau2MC.GetPdgCode()) == firstDaughterPDG))
        return motherVec;

    if (!dau1MC.isSecondary() || !dau2MC.isSecondary() || dau1MC.getMotherTrackId() != dau2MC.getMotherTrackId())
        return motherVec;
    auto v0MC = mcTracksMatrix[v0DauRefs[0][0]][dau1MC.getMotherTrackId()];

    if (std::abs(v0MC.GetPdgCode()) != motherPDG)
        return motherVec;

    motherVec = {v0DauRefs[0][0], dau1MC.getMotherTrackId()};
    return motherVec;
}

double calcMass(const V0 &v0, int firstV0dauPDG, int secondV0dauPDG)
{
    double dauMass[2];
    int dauCharges[2];
    auto firstV0dau = TDatabasePDG::Instance()->GetParticle(firstV0dauPDG);
    auto secondV0dau = TDatabasePDG::Instance()->GetParticle(secondV0dauPDG);
    if (firstV0dauPDG == 1000020030 || secondV0dauPDG == 1000020030)
    {
        int indexHe3 = firstV0dauPDG == 1000020030 ? 0 : 1;
        dauMass[indexHe3] = 2.808391;
        dauCharges[indexHe3] = 2;
        dauMass[1 - indexHe3] = firstV0dauPDG == 1000020030 ? secondV0dau->Mass() : firstV0dau->Mass();
        dauCharges[1 - indexHe3] = firstV0dauPDG == 1000020030 ? secondV0dau->Charge() / 3 : firstV0dau->Charge() / 3;
    }

    else
    {
        dauMass[0] = firstV0dau->Mass();
        dauCharges[0] = firstV0dau->Charge();
        dauMass[1] = secondV0dau->Mass();
        dauCharges[1] = secondV0dau->Charge();
    }

    std::vector<o2::dataformats::V0::Track> dauTracks = {v0.getProng(0), v0.getProng(1)};
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

double calcV0alpha(const V0 &v0)
{
    std::array<float, 3> fV0mom, fPmom, fNmom = {0, 0, 0};
    v0.getProng(0).getPxPyPzGlo(fPmom);
    v0.getProng(1).getPxPyPzGlo(fNmom);
    v0.getPxPyPzGlo(fV0mom);

    TVector3 momNeg(fNmom[0], fNmom[1], fNmom[2]);
    TVector3 momPos(fPmom[0], fPmom[1], fPmom[2]);
    TVector3 momTot(fV0mom[0], fV0mom[1], fV0mom[2]);

    Double_t lQlNeg = momNeg.Dot(momTot) / momTot.Mag();
    Double_t lQlPos = momPos.Dot(momTot) / momTot.Mag();

    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}