#if !defined(CLING) || defined(ROOTCLING)

#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventLabel.h"
#include "SimulationDataFormat/MCTrack.h"

#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITStracking/IOUtils.h"
#include "DataFormatsParameters/GRPObject.h"

#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include <TLorentzVector.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TSystemDirectory.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "StrangenessTracking/StrangenessTracker.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "ReconstructionDataFormats/StrangeTrack.h"
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

const int motherPDG = 3312;
const int firstDaughterPDG = 3122; // pdg of the V0
const int secondDaughterPDG = 211; // pdg of the bachelor
const int firstV0dauPDG = 2212;    // pdg of the V0 daughter
const int secondV0dauPDG = 211;    // pdg of the V0 daughter

struct GPart
{
    float gPt, gCt, gR2, recoPt, recoR2, recoMass, ITStrackPt, trackedPt, trackedR2, trackedMass;
    float recoDCAXY, recoDCAZ, trackedDCAXY, trackedDCAZ;
    float recoPvX, recoPvY, recoPvZ;
    bool isTrueVertex = false, isReconstructed = false, isStTracked = false, isFakeMatching = false, hasITSTrack = false, isDuplicated = false;
    bool isTrueCasc = false;
    int itsRef = -1;
};

bool fillStrangeTrackInfo(GPart &part, std::vector<StrangeTrack> *strangeTrackVec, int cascID, o2::dataformats::PrimaryVertex &pv);
std::array<int, 2> matchCascToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec, Cascade &casc);
double calcCascAlpha(const Cascade &cascade);
double calcMass(const Cascade &casc, int firstV0dauPDG, int secondV0dauPDG);

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
        if (file.substr(0, 2) == "tf")
        {
            // if (stoi(file.substr(2)) > 100)
            //     continue;
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

    // dump the gPart structures to a tree
    TFile outFile = TFile(Form("TrackedCascTree%s.root", outSuffix.data()), "recreate");
    TTree *outTree = new TTree("CascTree", "CascTree");
    // declare branch
    GPart outputPart;
    outTree->Branch("Casc", &outputPart);

    // Geometry
    o2::base::GeometryManager::loadGeometry(dirs[0] + "/o2sim_geometry.root");
    auto gman = o2::its::GeometryTGeo::Instance();
    gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));

    const auto grp = o2::parameters::GRPObject::loadFrom(dirs[0] + "/" + "o2sim_grp.root");
    o2::base::Propagator::initFieldFromGRP(grp);
    auto propagator = o2::base::Propagator::Instance();
    LOG(info) << "Propagator mag field: " << int(propagator->getNominalBz()) << " kG";

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

        // bkg vector
        std::vector<GPart> bkgVec;

        // fill MC matrix
        std::vector<std::vector<MCTrack>> mcTracksMatrix;
        std::vector<std::vector<GPart>> gPartMatrix;

        auto nev = treeMCTracks->GetEntriesFast();
        mcTracksMatrix.resize(nev);
        gPartMatrix.resize(nev);
        for (int n = 0; n < nev; n++)
        { // loop over MC events
            treeMCTracks->GetEvent(n);
            mcTracksMatrix[n].resize(MCtracks->size());
            gPartMatrix[n].resize(MCtracks->size());
            for (unsigned int mcI{0}; mcI < MCtracks->size(); ++mcI)
            {
                GPart gPart;
                auto &mcTrack = MCtracks->at(mcI);

                if (abs(mcTrack.GetPdgCode()) == motherPDG)
                {
                    gPart.gPt = mcTrack.GetPt();
                    gPart.gR2 = strUtils.calcR2(MCtracks, mcTrack);
                    gPart.gCt = strUtils.calcLifetime(MCtracks, mcTrack);
                }
                else
                {
                    gPart.gPt = -1;
                    gPart.gR2 = -1;
                    gPart.gCt = -1;
                }

                gPartMatrix[n][mcI] = gPart;
                mcTracksMatrix[n][mcI] = mcTrack;
            }
        }
        LOG(info) << "MC matrix filled";

        // Starting matching Cascades and ITS tracks
        int counterV0 = 0;
        for (int frame = 0; frame < treeITS->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeSecondaries->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeITSclus->GetEvent(frame) || !treeTPCTRD->GetEvent(frame) || !treeStrangeTracks->GetEvent(frame) ||
                !treeITSTPCTRD->GetEvent(frame) || !treeTPCTRDTOF->GetEvent(frame) || !treeITSTPCTRDTOF->GetEvent(frame) || !treePrimaryVertex->GetEvent(frame))
                continue;

            for (unsigned int iCascVec = 0; iCascVec < cascVec->size(); iCascVec++)
            {

                auto &casc = cascVec->at(iCascVec);
                bool isMatter = calcCascAlpha(casc) > 0;
                auto cascMCref = matchCascToMC(mcTracksMatrix, map, v0Vec, casc);

                auto &primaryVertex = primVertices->at(casc.getVertexID());
                auto &pvLabel = pvMcArr->at(casc.getVertexID());

                if (cascMCref[0] != -1 && cascMCref[1] != -1)
                {

                    auto &mcCasc = gPartMatrix[cascMCref[0]][cascMCref[1]];
                    if (pvLabel.getEventID() == cascMCref[0])
                        mcCasc.isTrueVertex = true;

                    mcCasc.isReconstructed = true;

                    mcCasc.recoPt = casc.getPt();
                    mcCasc.recoR2 = sqrt(casc.calcR2());
                    mcCasc.recoMass = calcMass(casc, firstDaughterPDG, secondDaughterPDG);

                    mcCasc.recoPvX = primaryVertex.getX();
                    mcCasc.recoPvY = primaryVertex.getY();
                    mcCasc.recoPvZ = primaryVertex.getZ();

                    gpu::gpustd::array<float, 2> dca{-999.f, -999.f};
                    if (propagator->propagateToDCA(primaryVertex.getXYZ(), casc, propagator->getNominalBz(), 2.f, o2::base::PropagatorF::MatCorrType::USEMatCorrLUT, &dca))
                    {
                        mcCasc.recoDCAXY = dca[0];
                        mcCasc.recoDCAZ = dca[1];
                    }

                    // look for any strange track
                    auto isStrangeTracked = fillStrangeTrackInfo(mcCasc, strangeTrackVec, iCascVec, primaryVertex);
                    if (isStrangeTracked)
                    {

                        // match ITS track to MC
                        auto &labITS = labITSvec->at(mcCasc.itsRef);
                        auto ITSref = strUtils.matchITStracktoMC(mcTracksMatrix, labITS);

                        // check fake matching
                        if (ITSref[0] != cascMCref[0] || ITSref[1] != cascMCref[1])
                        {
                            mcCasc.isFakeMatching = true;
                        }
                    }
                }

                else // this will fill the background vector
                {
                    GPart bkgPart;
                    bkgPart.isReconstructed = true;
                    bkgPart.recoPt = casc.getPt();
                    bkgPart.recoR2 = sqrt(casc.calcR2());
                    bkgPart.recoMass = calcMass(casc, firstDaughterPDG, secondDaughterPDG);
                    bkgPart.gPt = -1;
                    bkgPart.gR2 = -1;
                    bkgPart.gCt = -1;
                    auto isStrangeTracked = fillStrangeTrackInfo(bkgPart, strangeTrackVec, iCascVec, primaryVertex);
                    bkgPart.isFakeMatching = true;
                    bkgVec.push_back(bkgPart);
                }
            }
            LOG(info) << "Cascades matched";

            // Finally, Check if the ITS track is available. Relevant only for signal

            for (unsigned int iITStrack = 0; iITStrack < ITStracks->size(); iITStrack++)
            {
                std::array<int, 2> ITSref = {-1, 1};
                auto &labITS = (*labITSvec)[iITStrack];
                auto &trackIdx = (*ITSTrackClusIdx)[iITStrack];

                ITSref = strUtils.matchITStracktoMC(mcTracksMatrix, labITS);
                auto &ITStrack = ITStracks->at(iITStrack);

                if (ITSref[0] != -1 && ITSref[1] != -1)
                {
                    auto &mcCasc = gPartMatrix[ITSref[0]][ITSref[1]];
                    mcCasc.hasITSTrack = true;
                    mcCasc.ITStrackPt = ITStrack.getPt();
                }
            }

            // first loop over the signal
            for (auto &part : gPartMatrix)
            {
                for (auto &mcPart : part)
                {
                    if (mcPart.gPt > 0)
                    {
                        outputPart = mcPart;
                        outTree->Fill();
                    }
                }
            }

            // now loop over the background
            for (auto &part : bkgVec)
            {
                outputPart = part;
                outTree->Fill();
            }
        }
    }
    outFile.cd();
    outTree->Write();
};

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
        LOG(info) << dau1MC.GetPdgCode() << " " << dau2MC.GetPdgCode() << " " << dau1MC.isSecondary() << " " << dau2MC.isSecondary() << " " << dau1MC.getMotherTrackId() << " " << dau2MC.getMotherTrackId();
        return motherVec;
    }

    auto v0MC = mcTracksMatrix[v0DauRefs[0][0]][dau1MC.getMotherTrackId()];
    auto &bachMC = mcTracksMatrix[bachRef[0]][bachRef[1]];

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

    auto firstV0dau = TDatabasePDG::Instance()->GetParticle(firstDauPdg);
    auto secondV0dau = TDatabasePDG::Instance()->GetParticle(secondDauPdg);

    double dauMass[2] = {firstV0dau->Mass(), secondV0dau->Mass()};

    std::vector<o2::dataformats::Cascade::Track> dauTracks = {casc.getV0Track(), casc.getBachelorTrack()};
    TLorentzVector moth, prong;
    std::array<float, 3> p;

    for (int i = 0; i < 2; i++)
    {
        auto &track = dauTracks[i];
        auto &mass = dauMass[i];
        track.getPxPyPzGlo(p);
        prong.SetVectM({p[0], p[1], p[2]}, mass);
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

bool fillStrangeTrackInfo(GPart &part, std::vector<StrangeTrack> *strangeTrackVec, int cascID, o2::dataformats::PrimaryVertex &pv)
{
    // check whether the corresponding strange track is found: both signal and background checked
    bool isFound = false;
    int index = -1;
    for (unsigned int iStTr = 0; iStTr < strangeTrackVec->size(); iStTr++)
    {
        auto &strangeTrack = strangeTrackVec->at(iStTr);
        if (int(strangeTrack.mDecayRef) == cascID)
        {
            isFound = true;
            index = iStTr;
            break;
        }
    }
    if (!isFound)
        return false;

    auto &StrangeTrack = strangeTrackVec->at(index);
    part.isStTracked = true;
    part.trackedPt = sqrt(StrangeTrack.mDecayMom[0] * StrangeTrack.mDecayMom[0] + StrangeTrack.mDecayMom[1] * StrangeTrack.mDecayMom[1]);
    part.trackedR2 = StrangeTrack.mDecayVtx[0] * StrangeTrack.mDecayVtx[0] + StrangeTrack.mDecayVtx[1] * StrangeTrack.mDecayVtx[1];

    gpu::gpustd::array<float, 2> dcaT{-999.f, -999.f};
    auto propagator = o2::base::Propagator::Instance();
    if (propagator->propagateToDCA(pv.getXYZ(), StrangeTrack.mMother, propagator->getNominalBz(), 2.f, o2::base::PropagatorF::MatCorrType::USEMatCorrLUT, &dcaT))
    {
        part.trackedDCAXY = dcaT[0];
        part.trackedDCAZ = dcaT[1];
    }

    part.itsRef = StrangeTrack.mITSRef;

    return true;
}