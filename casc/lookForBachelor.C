#if !defined(CLING) || defined(ROOTCLING)
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Cascade.h"
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
// #include "StrangenessTracking/StrangenessTracker.h"
#include "cascadeStudyUtils.h"

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
// using StrangeTrack = o2::strangeness_tracking::StrangeTrack;

std::array<int, 2> matchCascToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec, Cascade &casc, bool &isV0reco);
std::array<int, 2> checkV0mother(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, V0 &v0);

void lookForBachelor()
{
    std::string path = "/data/fmazzasc/its_data/sim/xi_test/";
    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    std::vector<TString> kine_files;

    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (file.substr(0, 3) == "tf1")
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

    TH1D *hGenXiRadius = new TH1D("gen_casc_r", "; Radius (cm); Counts", 40, 0., 40.);
    TH1D *hGenXiLifetime = new TH1D("gen_casc_lifetime", "; Radius (cm); Counts", 40, 0., 40.);

    TH1D *hGenXiMom = new TH1D("gen_casc_pt", "; #it{p}_{T}^{gen} (GeV/#it{c}); Counts", 40, 1, 10);
    TH1D *hGenEta = new TH1D("gen_casc_eta", "; #eta; Counts", 40, -2.5, 2.5);

    TH1D *hRecXiRadius = new TH1D("rec_casc_r", "; Radius (cm); Counts", 40, 0., 40.);
    TH1D *hRecXiMom = new TH1D("rec_casc_pt", "; #it{p}_{T}^{gen} (GeV/#it{c}); Counts", 40, 1, 10);
    TH1D *hFakeAssocCounter = new TH1D("Fake assoc counter", ";Fake assoc counter; Counts", 1, 0.5, 1.5);

    TH1D *hRecXiCounter = new TH1D("Rec Xi counter", "; ; Counts", 1, 0.5, 1.5);
    TH1D *hCascMomsInV0s = new TH1D("Rec V0s from Xi counter", "; ; Counts", 1, 0.5, 1.5);
    TH1D *hFindableBachfromXiCounter = new TH1D("Rec bachelors from Xi counter", "; ; Counts", 1, 0.5, 1.5);

    TH1D *hRecCascInvMass = new TH1D("Rec Casc InvMass", "; M (GeV/c^{2}); Counts", 150, 1.2, 1.5);
    TH1D *hStrangeTrackInvMass = new TH1D("Strange track inv mass", ";M (GeV/c^{2}); Counts", 150, 1.2, 1.5);
    TH1D *hXiStats = new TH1D("cascade_stats", "; ; Counts", 3, 0.5, 3.5);

    TH1D *hRecXiRad = new TH1D("r_trackable", ";Radius_{trackable} (cm); Counts", 20, 4, 40);
    TH1D *hTrackedXiR2 = new TH1D("R2 resolution after tracking", ";R2^{rec} (cm); Counts", 20, 4, 40);

    int counter = 0;
    for (unsigned int i = 0; i < dirs.size(); i++)
    {
        auto &dir = dirs[i];
        auto &kine_file = kine_files[i];
        LOG(info) << "Processing " << dir;
        LOG(info) << "kine file " << kine_file;
        // Files
        auto fMCTracks = TFile::Open((TString(dir + "/") + kine_file));
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
                    auto &motherTrack = mcTracksMatrix[n][mcI];
                    hGenXiCounter->Fill(1);
                    hGenXiMom->Fill(motherTrack.GetPt());
                    hGenXiRadius->Fill(calcDecLength(MCtracks, motherTrack, bachPDG));
                    hGenEta->Fill(motherTrack.GetEta());
                    hGenXiLifetime->Fill(calcLifetime(MCtracks, motherTrack, bachPDG));
                }
            }
        }

        // Starting matching Cascades and ITS tracks
        int counterV0 = 0;
        for (int frame = 0; frame < treeITS->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeSecondaries->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeITSclus->GetEvent(frame) || !treeTPCTRD->GetEvent(frame))
                continue;

            for (unsigned int iV0 = 0; iV0 < v0Vec->size(); iV0++)
            {
                auto &v0 = v0Vec->at(iV0);
                auto v0MCref = checkV0mother(mcTracksMatrix, map, v0);
                if (v0MCref[0] == -1 || v0MCref[1] == -1)
                    continue;
                // start looking for bachelor tracks
                LOG(info) << "---------------------------";
                LOG(info) << "V0, R2: " << v0.calcR2() << ", Pt: " << v0.getPt();

                for (unsigned int iGlo{0}; iGlo < map["ITS-TPC-TOF"]->size(); iGlo++)
                {
                    auto &gloLab = map["ITS-TPC-TOF"]->at(iGlo);
                    auto evID = gloLab.getEventID();
                    if (evID != v0MCref[0])
                        continue;
                    auto trkID = gloLab.getTrackID();
                    auto motherID = mcTracksMatrix[evID][trkID].getMotherTrackId();
                    if (motherID != v0MCref[1])
                        continue;
                    LOG(info) << "Bachelor ITS-TPC-TOF found: " << iGlo;
                }

                // Start looking for bachelor ITS-TPC tracks
                for (unsigned int iGlo{0}; iGlo < map["ITS-TPC"]->size(); iGlo++)
                {
                    auto &gloLab = map["ITS-TPC"]->at(iGlo);
                    auto evID = gloLab.getEventID();
                    if (evID != v0MCref[0])
                        continue;
                    auto trkID = gloLab.getTrackID();
                    auto motherID = mcTracksMatrix[evID][trkID].getMotherTrackId();
                    if (motherID != v0MCref[1])
                        continue;
                    LOG(info) << "Bachelor ITS-TPC found: " << iGlo;
                }

                // Start looking for bachelor TPC only tracks
                for (unsigned int iGlo{0}; iGlo < map["TPC"]->size(); iGlo++)
                {
                    auto &gloLab = map["TPC"]->at(iGlo);
                    auto evID = gloLab.getEventID();
                    if (evID != v0MCref[0])
                        continue;
                    auto trkID = gloLab.getTrackID();
                    auto motherID = mcTracksMatrix[evID][trkID].getMotherTrackId();
                    if (motherID != v0MCref[1])
                        continue;
                    LOG(info) << "Bachelor TPC found: " << iGlo;
                }
            }
        }
    }
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

    auto v0MC = mcTracksMatrix[v0DauRefs[0][0]][dau1MC.getMotherTrackId()];
    auto cascMC = mcTracksMatrix[v0DauRefs[0][0]][v0MC.getMotherTrackId()];
    if (std::abs(cascMC.GetPdgCode()) != motherPDG)
        return motherVec;

    motherVec = {v0DauRefs[0][0], v0MC.getMotherTrackId()};
    return motherVec;
}
