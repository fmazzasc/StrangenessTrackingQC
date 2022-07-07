#if !defined(CLING) || defined(ROOTCLING)
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ITSMFTSimulation/Hit.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITStracking/IOUtils.h"

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

#endif

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

const int motherPDG = 3312;
const int v0PDG = 3122;
const int bachPDG = 211;
const int firstV0dauPDG = 2212;
const int secondV0dauPDG = 211;

double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);

std::array<int, 2> matchCascToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec, Cascade &casc, bool &isV0reco);
std::array<int, 2> checkV0mother(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, V0 &v0);
std::array<int, 2> matchITStracktoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel);
std::array<int, 2> matchCompLabelToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel);
std::vector<ITSCluster> getTrackClusters(const o2::its::TrackITS &ITStrack, const std::vector<ITSCluster> &ITSClustersArray, std::vector<int> *ITSTrackClusIdx);

void cascadeStudy()
{
    std::string path = "/data/fmazzasc/its_data/sim/xi/";
    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    std::vector<TString> kine_files;

    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (file.substr(0, 3) == "tf3")
        {
            dirs.push_back(path + file);
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

    TH2D *hHyperhisto = new TH2D("histo hyperV0", ";#it{p}_{T} (GeV/#it{c}); Radius^2 (cm) ; Counts", 20, 1, 10, 30, 1, 900);
    TH1D *hResV0histo = new TH1D("pT resolution before hypertracking", ";(#it{p}_{T}^{gen} - #it{p}_{T}^{rec})/#it{p}_{T}^{gen}; Counts", 20, -0.2, 0.2);
    TH1D *hResHyperhisto = new TH1D("pT resolution after hypertracking", ";(#it{p}_{T}^{gen} - #it{p}_{T}^{rec})/#it{p}_{T}^{gen}; Counts", 20, -0.2, 0.2);
    TH1D *hResCascR2 = new TH1D("R2 resolution before tracking", ";(R2^{gen} - R2^{rec})/R2^{gen}; Counts", 20, -0.2, 0.2);
    TH1D *hResCascTrackedR2 = new TH1D("R2 resolution after tracking", ";(R2^{gen} - R2^{rec})/R2^{gen}; Counts", 20, -0.2, 0.2);
    TH1D *hCascCounter = new TH1D("Casc counter", ";Casc counter; Counts", 1, 0.5, 1.5);
    TH1D *hGenXiCounter = new TH1D("Gen Casc counter", ";Casc counter; Counts", 1, 0.5, 1.5);
    TH1D *hFakeAssocCounter = new TH1D("Fake assoc counter", ";Fake assoc counter; Counts", 1, 0.5, 1.5);

    TH1D *hRecXiCounter = new TH1D("Rec Xi counter", "; ; Counts", 1, 0.5, 1.5);
    TH1D *hCascMomsInV0s = new TH1D("Rec V0s from Xi counter", "; ; Counts", 1, 0.5, 1.5);
    TH1D *hFindableBachfromXiCounter = new TH1D("Rec bachelors from Xi counter", "; ; Counts", 1, 0.5, 1.5);

    TH1D *hXiStats = new TH1D("cascade_stats", "; ; Counts", 3, 0.5, 3.5);

    int counter = 0;
    for (unsigned int i = 0; i < dirs.size(); i++)
    {
        auto &dir = dirs[i];
        auto &kine_file = kine_files[i];
        LOG(info) << "Processing " << dir;
        LOG(info) << "kine file " << kine_file;
        // Files
        auto fMCTracks = TFile::Open((TString(dir + "/") + kine_file));
        auto fStrangeTracks = TFile::Open((dir + "/o2_strange_tracks.root").data());
        auto fSecondaries = TFile::Open((dir + "/o2_secondary_vertex.root").data());
        auto fITSTPC = TFile::Open((dir + "/o2match_itstpc.root").data());
        auto fTPCTOF = TFile::Open((dir + "/o2match_tof_tpc.root").data());
        auto fTPCTRD = TFile::Open((dir + "/trdmatches_tpc.root").data());
        auto fITSTPCTOF = TFile::Open((dir + "/o2match_tof_itstpc.root").data());
        auto fITS = TFile::Open((dir + "/o2trac_its.root").data());
        auto fClusITS = TFile::Open((dir + "/o2clus_its.root").data());
        auto fTPC = TFile::Open((dir + "/tpctracks.root").data());

        // Geometry
        o2::base::GeometryManager::loadGeometry(dir + "/o2sim_geometry.root");

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeStrangeTracks = (TTree *)fStrangeTracks->Get("o2sim");
        auto treeSecondaries = (TTree *)fSecondaries->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
        auto treeITSTPCTOF = (TTree *)fITSTPCTOF->Get("matchTOF");
        auto treeTPCTOF = (TTree *)fTPCTOF->Get("matchTOF");
        auto treeTPCTRD = (TTree *)fTPCTRD->Get("tracksTRD");

        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSclus = (TTree *)fClusITS->Get("o2sim");
        auto treeTPC = (TTree *)fTPC->Get("tpcrec");

        // MC Tracks
        std::vector<o2::MCTrack> *MCtracks = nullptr;
        std::vector<o2::itsmft::Hit> *ITSHits = nullptr;

        // Hypertracks
        std::vector<StrangeTrack> *strangeTrackVec = nullptr;
        std::vector<int> *strangeITSrefVec = nullptr;
        std::vector<int> *strangeDecRefVec = nullptr;
        std::vector<o2::track::TrackParametrizationWithError<float>> *hypertrackVec = nullptr;
        std::vector<o2::strangeness_tracking::ClusAttachments> *nAttachments = nullptr;

        // Secondary Vertices
        std::vector<Cascade> *cascVec = nullptr;
        std::vector<V0> *v0Vec = nullptr;

        // ITS tracks
        std::vector<o2::its::TrackITS> *ITStracks = nullptr;

        // Labels
        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTRDvec = nullptr;

        // Clusters
        std::vector<CompClusterExt> *ITSclus = nullptr;
        o2::dataformats::MCTruthContainer<o2::MCCompLabel> *clusLabArr = nullptr;
        std::vector<int> *ITSTrackClusIdx = nullptr;
        std::vector<unsigned char> *ITSpatt = nullptr;

        // Setting branches
        treeStrangeTracks->SetBranchAddress("StrangeTracks", &strangeTrackVec);
        treeStrangeTracks->SetBranchAddress("ITSTrackRefs", &strangeITSrefVec);
        treeStrangeTracks->SetBranchAddress("DecayTrackRefs", &strangeDecRefVec);
        treeStrangeTracks->SetBranchAddress("ClusUpdates", &nAttachments);

        treeSecondaries->SetBranchAddress("Cascades", &cascVec);
        treeSecondaries->SetBranchAddress("V0s", &v0Vec);

        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);

        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeITS->SetBranchAddress("ITSTrack", &ITStracks);
        treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
        treeTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTOFvec);
        treeTPCTRD->SetBranchAddress("labels", &labTPCTRDvec);

        treeITSTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTOFvec);

        treeITS->SetBranchAddress("ITSTrackClusIdx", &ITSTrackClusIdx);
        treeITSclus->SetBranchAddress("ITSClusterComp", &ITSclus);
        treeITSclus->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);

        // define detector map
        std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}, {"TPC", labTPCvec}, {"ITS-TPC", labITSTPCvec}, {"TPC-TOF", labTPCTOFvec}, {"TPC-TRD", labTPCTRDvec}, {"ITS-TPC-TOF", labITSTPCTOFvec}};

        // load geometry
        auto gman = o2::its::GeometryTGeo::Instance();
        gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));

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
                if (abs(MCtracks->at(mcI).GetPdgCode()) == motherPDG)
                {
                    hGenXiCounter->Fill(1);
                }
            }
        }

        // Starting matching Cascades and ITS tracks
        int counterV0 = 0;
        for (int frame = 0; frame < treeITS->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeSecondaries->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeITSclus->GetEvent(frame) || !treeTPCTRD->GetEvent(frame) || !treeStrangeTracks->GetEvent(frame))
                continue;

            for (unsigned int iCascVec = 0; iCascVec < cascVec->size(); iCascVec++)
            {
                hCascCounter->Fill(1);
                auto &casc = cascVec->at(iCascVec);

                bool isV0reco = false;
                auto cascMCref = matchCascToMC(mcTracksMatrix, map, v0Vec, casc, isV0reco);

                if (cascMCref[0] == -1 || cascMCref[1] == -1)
                    continue;

                LOG(info) << "Casc found!";

                hRecXiCounter->Fill(1);

                auto &mcCasc = mcTracksMatrix[cascMCref[0]][cascMCref[1]];
                // LOG(info) << "---------------------------------------------------";
                // LOG(info) << "Cascade found, PDG: " << mcCasc.GetPdgCode() << ", MC pT: " << mcCasc.GetPt() << ", reco Pt: " << casc.getPt();

                // Matching ITS tracks to MC tracks and V0
                std::array<int, 2> ITSref = {-1, 1};
                o2::its::TrackITS ITStrack;
                std::array<std::array<int, 2>, 7> clsRef;

                int iTrack = -1;
                bool isMatched = false;

                for (unsigned int iITStrack = 0; iITStrack < ITStracks->size(); iITStrack++)
                {
                    auto &labITS = (*labITSvec)[iITStrack];
                    auto &trackIdx = (*ITSTrackClusIdx)[iITStrack];

                    ITSref = matchITStracktoMC(mcTracksMatrix, labITS);
                    ITStrack = (*ITStracks)[iITStrack];

                    if (ITSref[0] == cascMCref[0] && ITSref[1] == cascMCref[1])
                    {
                        LOG(info) << "ITS track found! " << ITStrack.getPt();
                        hXiStats->Fill(1);

                        auto firstClus = ITStrack.getFirstClusterEntry();
                        auto ncl = ITStrack.getNumberOfClusters();
                        for (int icl = 0; icl < ncl; icl++)
                        {
                            auto &labCls = (clusLabArr->getLabels(ITSTrackClusIdx->at(firstClus + icl)))[0];
                            auto &clus = (*ITSclus)[(*ITSTrackClusIdx)[firstClus + icl]];
                            auto layer = gman->getLayer(clus.getSensorID());
                            clsRef[layer] = matchCompLabelToMC(mcTracksMatrix, labCls);
                            if (clsRef[layer][0] > -1 && clsRef[layer][1] > -1)
                                LOG(info) << "Layer: " << layer << "PDG: " << mcTracksMatrix[clsRef[layer][0]][clsRef[layer][1]].GetPdgCode();
                            else
                                LOG(info) << "Layer: " << layer << ", No valid cluster ref";
                        }
                    }
                }
            }

            for (unsigned int iV0 = 0; iV0 < v0Vec->size(); iV0++)
            {
                auto &v0 = v0Vec->at(iV0);
                auto v0MCref = checkV0mother(mcTracksMatrix, map, v0);
                if (v0MCref[0] == -1 || v0MCref[1] == -1)
                    continue;
                LOG(info) << "V0 from Xi found!";
                hCascMomsInV0s->Fill(1);
            }

            // LOG(info) << "+++++++++++++++++++++++++++++++++++++++++++++";

            // for (unsigned int iHyperVec = 0; iHyperVec < strangeTrackVec->size(); iHyperVec++)
            // {

            //     auto &strangeTrack = strangeTrackVec->at(iHyperVec);
            //     if (!strangeTrack.isCascade)
            //     {
            //         continue;
            //     }

            //     auto &hyperChi2 = strangeTrack.mMatchChi2;
            //     auto &clusAttachments = nAttachments->at(iHyperVec);
            //     auto &ITStrack = ITStracks->at(strangeITSrefVec->at(iHyperVec));
            //     auto &ITStrackLab = labITSvec->at(strangeITSrefVec->at(iHyperVec));

            //     auto clusAttArr = clusAttachments.arr;
            //     auto &casc = cascVec->at(strangeDecRefVec->at(iHyperVec));

            //     auto ITStrackRef = matchCompLabelToMC(mcTracksMatrix, ITStrackLab);
            //     bool isV0reco = false;
            //     auto cascMCref = matchCascToMC(mcTracksMatrix, map, v0Vec, casc, isV0reco);

            //     if (cascMCref[0] == -1 || cascMCref[1] == -1 || ITStrackRef[0] == -1 || ITStrackRef[1] == -1)
            //         continue;

            //     if (ITStrackRef[0] == cascMCref[0] && ITStrackRef[1] == cascMCref[1])
            //     {
            //         hXiStats->Fill(2);
            //     }
            //     else
            //     {
            //         hFakeAssocCounter->Fill(1);
            //     }
            // }
        }
    }

    auto outFile = TFile("cascade_study.root", "recreate");

    hResV0histo->Write();
    hHyperhisto->Write();
    hResHyperhisto->Write();

    hGenXiCounter->Write();
    hCascCounter->Write();
    hXiStats->Write();
    hRecXiCounter->Write();
    hFakeAssocCounter->Write();
    hCascMomsInV0s->Write();
    outFile.Close();
}
std::array<int, 2> matchCascToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec, Cascade &casc, bool &isV0reco)
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
    }

    if (v0DauRefs[0][1] == -1 || v0DauRefs[1][1] == -1)
        return motherVec;

    auto &dau1MC = mcTracksMatrix[v0DauRefs[0][0]][v0DauRefs[0][1]];
    auto &dau2MC = mcTracksMatrix[v0DauRefs[1][0]][v0DauRefs[1][1]];

    if (!(std::abs(dau1MC.GetPdgCode()) == firstV0dauPDG && std::abs(dau2MC.GetPdgCode()) == secondV0dauPDG) && !(std::abs(dau1MC.GetPdgCode()) == secondV0dauPDG && std::abs(dau2MC.GetPdgCode()) == firstV0dauPDG))
        return motherVec;

    if (!dau1MC.isSecondary() || !dau2MC.isSecondary() || dau1MC.getMotherTrackId() != dau2MC.getMotherTrackId())
        return motherVec;

    LOG(info) << "-----------------------------------------";
    LOG(info) << "V0 daughters: " << dau1MC.GetPdgCode() << " " << dau2MC.GetPdgCode();
    LOG(info) << " Dau refs: " << dau1MC.getMotherTrackId() << " " << dau2MC.getMotherTrackId();

    auto v0MC = mcTracksMatrix[v0DauRefs[0][0]][dau1MC.getMotherTrackId()];
    auto &bachMC = mcTracksMatrix[bachRef[0]][bachRef[1]];
    LOG(info) << "V0 mother: " << v0MC.GetPdgCode() << ", Ref: " << v0MC.getMotherTrackId();
    LOG(info) << "bach PDG: " << bachMC.GetPdgCode() << ", Ref: " << bachMC.getMotherTrackId();

    if (v0MC.getMotherTrackId() >= 0)
    {
        LOG(info) << "PDG of Casc from V0:" << std::abs(mcTracksMatrix[v0DauRefs[0][0]][v0MC.getMotherTrackId()].GetPdgCode());

        if (std::abs(mcTracksMatrix[v0DauRefs[0][0]][v0MC.getMotherTrackId()].GetPdgCode()) == motherPDG)
            isV0reco = true;
    }

    if (std::abs(v0MC.GetPdgCode()) != v0PDG || !v0MC.isSecondary() || !bachMC.isSecondary())
        return motherVec;

    auto cascMC = mcTracksMatrix[v0DauRefs[0][0]][v0MC.getMotherTrackId()];
    if (v0MC.getMotherTrackId() != bachMC.getMotherTrackId() || std::abs(cascMC.GetPdgCode()) != motherPDG)
        return motherVec;

    LOG(info) << "Casc PDG: " << mcTracksMatrix[v0DauRefs[0][0]][v0MC.getMotherTrackId()].GetPdgCode();

    motherVec = {v0DauRefs[0][0], v0MC.getMotherTrackId()};
    return motherVec;
}

std::array<int, 2> checkV0mother(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, V0 &v0)
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

    if (!(std::abs(dau1MC.GetPdgCode()) == firstV0dauPDG && std::abs(dau2MC.GetPdgCode()) == secondV0dauPDG) && !(std::abs(dau1MC.GetPdgCode()) == secondV0dauPDG && std::abs(dau2MC.GetPdgCode()) == firstV0dauPDG))
        return motherVec;

    if (!dau1MC.isSecondary() || !dau2MC.isSecondary() || dau1MC.getMotherTrackId() != dau2MC.getMotherTrackId())
        return motherVec;

    LOG(info) << "-----------------------------------------";
    LOG(info) << "V0 daughters: " << dau1MC.GetPdgCode() << " " << dau2MC.GetPdgCode();
    LOG(info) << " Dau refs: " << dau1MC.getMotherTrackId() << " " << dau2MC.getMotherTrackId();

    auto v0MC = mcTracksMatrix[v0DauRefs[0][0]][dau1MC.getMotherTrackId()];
    LOG(info) << "V0 mother: " << v0MC.GetPdgCode() << ", Ref: " << v0MC.getMotherTrackId();

    auto cascMC = mcTracksMatrix[v0DauRefs[0][0]][v0MC.getMotherTrackId()];
    if (std::abs(cascMC.GetPdgCode()) != motherPDG)
        return motherVec;

    LOG(info) << "Casc from V0 PDG: " << mcTracksMatrix[v0DauRefs[0][0]][v0MC.getMotherTrackId()].GetPdgCode();

    motherVec = {v0DauRefs[0][0], v0MC.getMotherTrackId()};
    return motherVec;
}

std::array<int, 2> matchITStracktoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel)

{
    std::array<int, 2> outArray = {-1, -1};
    int trackID, evID, srcID;
    bool fake;
    ITSlabel.get(trackID, evID, srcID, fake);
    if (ITSlabel.isValid() && std::abs(mcTracksMatrix[evID][trackID].GetPdgCode()) == motherPDG)
    {
        outArray = {evID, trackID};
    }

    return outArray;
}

std::array<int, 2> matchCompLabelToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel compLabel)
{
    std::array<int, 2> compRef = {-1, -1};
    int trackID, evID, srcID;
    bool fake;
    compLabel.get(trackID, evID, srcID, fake);
    if (compLabel.isValid())
    {
        compRef = {evID, trackID};
    }
    return compRef;
}

std::vector<ITSCluster> getTrackClusters(const o2::its::TrackITS &ITStrack, const std::vector<ITSCluster> &ITSClustersArray, std::vector<int> *ITSTrackClusIdx)
{

    std::vector<ITSCluster> outVec;
    auto firstClus = ITStrack.getFirstClusterEntry();
    auto ncl = ITStrack.getNumberOfClusters();
    for (int icl = 0; icl < ncl; icl++)
    {
        outVec.push_back(ITSClustersArray[(*ITSTrackClusIdx)[firstClus + icl]]);
    }
    return outVec;
}

double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
{
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    for (auto iD{idStart}; iD < idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        if (std::abs(dauTrack.GetPdgCode()) == dauPDG)
        {
            auto decLength = (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) *
                                 (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) +
                             (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) *
                                 (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY());
            return decLength;
        }
    }
    return -1;
}