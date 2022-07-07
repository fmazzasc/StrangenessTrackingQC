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
#include "ReconstructionDataFormats/TrackTPCITS.h"
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
#include "cascadeStudyUtils.h"

#endif

using MCTrack = o2::MCTrack;
using CompClusterExt = o2::itsmft::CompClusterExt;
using ITSCluster = o2::BaseCluster<float>;


void checkCascadeDaughters()
{
    std::string source = "ITS-TPC";
    std::string path = "/data/fmazzasc/its_data/sim/xi/";
    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    std::vector<TString> kine_files;

    TH1D *h_decLength_bach = new TH1D("h_dec_rad_bach", "h_dec_rad_bach", 100, 0, 100);
    TH1D *h_decLength_v0_dau = new TH1D("h_dec_rad_v0_dat", "h_dec_rad_v0_dau", 100, 0, 100);
    TH1D *h_decLength_v0_dau_reco = new TH1D("h_dec_rad_v0_rec", "h_dec_rad_v0_rec", 100, 0, 100);

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
        // auto treeStrangeTracks = (TTree *)fStrangeTracks->Get("o2sim");
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

        // ITS tracks
        std::vector<o2::dataformats::TrackTPCITS> *ITSTPCtracks = nullptr;

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

        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeITSTPC->SetBranchAddress("TPCITS", &ITSTPCtracks);
        treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
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
                auto &track = MCtracks->at(mcI);
                if (std::abs(track.GetPdgCode()) == motherPDG)
                {
                    h_decLength_bach->Fill(calcDecLength(MCtracks, track, bachPDG));
                    h_decLength_v0_dau->Fill(calcDecLengthV0(MCtracks, track, firstV0dauPDG));
                }
            }
        }

        for (int frame = 0; frame < treeITSTPC->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeSecondaries->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeITSclus->GetEvent(frame) || !treeTPCTRD->GetEvent(frame))
                continue;

            const size_t sizeT = size_t(ITSTPCtracks->size());
            std::vector<std::vector<int>> refs;
            refs.resize(sizeT);

            // look for V0 from Cascade

            std::vector<std::array<int, 2>> v0FromCascade;

            for (unsigned int iTrack{0}; iTrack < map[source.data()]->size(); iTrack++)
            {
                auto &lab1 = map[source.data()]->at(iTrack);

                auto cascRef1 = matchCascDauToMC(mcTracksMatrix, lab1, firstV0dauPDG);
                if (cascRef1[0] == -1 || cascRef1[1] == -1)
                    continue;

                for (unsigned int jTrack{0}; jTrack < map[source.data()]->size(); jTrack++)
                {
                    auto &lab2 = map[source.data()]->at(jTrack);
                    auto cascRef2 = matchCascDauToMC(mcTracksMatrix, lab2, secondV0dauPDG);
                    if (cascRef2[0] != cascRef1[0] || cascRef2[1] != cascRef1[1])
                        continue;
                    v0FromCascade.push_back(cascRef2);
                }
            }
            LOG(info) << "Size of v0FromCascade: " << v0FromCascade.size();
            LOG(info) << "Looping on bachelors..";
            int counter = 0;
            for (unsigned int iTrack{0}; iTrack < map[source.data()]->size(); iTrack++)
            {
                auto &labBach = map[source.data()]->at(iTrack);
                for (auto &ref : v0FromCascade)
                {

                    auto refBach = matchBachToMC(mcTracksMatrix, labBach);
                    if (refBach[0] != ref[0] || refBach[1] != ref[1])
                        continue;

                    auto mcPart = mcTracksMatrix[ref[0]][ref[1]];

                    counter++;
                    h_decLength_v0_dau_reco->Fill(calcDecLengthV0(&mcTracksMatrix[ref[0]], mcPart, firstV0dauPDG));
                }
            }
            LOG(info) << "# cascades found: " << counter;
        }

        auto outFile = TFile("check_daughters_cascade.root", "RECREATE");
        h_decLength_bach->Write();
        h_decLength_v0_dau->Write();
        h_decLength_v0_dau_reco->Write();
        outFile.Close();
    }
}