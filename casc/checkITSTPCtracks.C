#if !defined(CLING) || defined(ROOTCLING)
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



void checkITSTPCtracks()
{
    TH1D *hXiCounter = new TH1D("hXiCounter", "hXiCounter", 1, 0, 1);

    TH2D *hPrFromV0Clusters = new TH2D("proton_clus", "; ; Layer", 3, -0.5, 2.5, 7, -0.5, 6.5);
    TH2D *hPiFromV0Clusters = new TH2D("pi_clus", "; ; ; Layer", 3, -0.5, 2.5, 7, -0.5, 6.5);
    TH2D *hBachClusters = new TH2D("bach_clus", "; ; ; Layer", 3, -0.5, 2.5, 7, -0.5, 6.5);

    TH1D *hPrVsRadius = new TH1D("proton_rad", "; Radius; Counts", 100, 0, 40);
    TH1D *hPrVsRadiusAB = new TH1D("proton_rad_ab", "; Radius; Counts", 100, 0, 40);
    TH1D *hPrVsRadiusFake = new TH1D("proton_rad_fake", "; Radius; Counts", 100, 0, 40);
    TH1D *hPrVsRadiusABFake = new TH1D("proton_rad_ab_fake", "; Radius; Counts", 100, 0, 40);

    TH1D *hPiVsRadius = new TH1D("pion_rad", "; Radius; Counts", 100, 0, 40);
    TH1D *hPiVsRadiusAB = new TH1D("pion_rad_ab", "; Radius; Counts", 100, 0, 40);
    TH1D *hPiVsRadiusFake = new TH1D("pion_rad_fake", "; Radius; Counts", 100, 0, 40);
    TH1D *hPiVsRadiusABFake = new TH1D("pi_rad_ab_fake", "; Radius; Counts", 100, 0, 40);

    TH1D *hBachVsRadius = new TH1D("bach_rad", "; Radius; Counts", 100, 0, 40);
    TH1D *hBachVsRadiusAB = new TH1D("bach_rad_ab", "; Radius; Counts", 100, 0, 40);
    TH1D *hBachVsRadiusFake = new TH1D("bach_rad_fake", "; Radius; Counts", 100, 0, 40);
    TH1D *hBachVsRadiusABFake = new TH1D("bach_rad_ab_fake", "; Radius; Counts", 100, 0, 40);

    std::map<int, TH2D *> mapClus = {{kPr, hPrFromV0Clusters}, {kPi, hPiFromV0Clusters}, {kBach, hBachClusters}};
    std::map<int, TH1D *> mapRad = {{kPr, hPrVsRadius}, {kPi, hPiVsRadius}, {kBach, hBachVsRadius}};
    std::map<int, TH1D *> mapRadAB = {{kPr, hPrVsRadiusAB}, {kPi, hPiVsRadiusAB}, {kBach, hBachVsRadiusAB}};
    std::map<int, TH1D *> mapRadFake = {{kPr, hPrVsRadiusFake}, {kPi, hPiVsRadiusFake}, {kBach, hBachVsRadiusFake}};
    std::map<int, TH1D *> mapRadABFake = {{kPr, hPrVsRadiusABFake}, {kPi, hPiVsRadiusABFake}, {kBach, hBachVsRadiusABFake}};

    const char *labHits[3] = {"Free Cluster", "Full track", "Wrongly tracked cluster"};
    for (auto iLab{0}; iLab < 3; ++iLab)
    {
        hPrFromV0Clusters->GetXaxis()->SetBinLabel(iLab + 1, labHits[iLab]);
        hPiFromV0Clusters->GetXaxis()->SetBinLabel(iLab + 1, labHits[iLab]);
        hBachClusters->GetXaxis()->SetBinLabel(iLab + 1, labHits[iLab]);
    }

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

    for (unsigned int i = 0; i < dirs.size(); i++)
    {
        auto &dir = dirs[i];
        auto &kine_file = kine_files[i];
        LOG(info) << "Processing " << dir;
        LOG(info) << "kine file " << kine_file;
        // Files
        auto fMCTracks = TFile::Open((TString(dir + "/") + kine_file));
        // auto fStrangeTracks = TFile::Open((dir + "/o2_strange_tracks.root").data());
        auto fSecondaries = TFile::Open((dir + "/o2_secondary_vertex.root").data());
        auto fITSTPC = TFile::Open((dir + "/o2match_itstpc.root").data());
        auto fITS = TFile::Open((dir + "/o2trac_its.root").data());
        auto fClusITS = TFile::Open((dir + "/o2clus_its.root").data());
        auto fTPC = TFile::Open((dir + "/tpctracks.root").data());

        // Geometry
        o2::base::GeometryManager::loadGeometry(dir + "/o2sim_geometry.root");

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");

        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSclus = (TTree *)fClusITS->Get("o2sim");
        auto treeTPC = (TTree *)fTPC->Get("tpcrec");

        // MC Tracks
        std::vector<o2::MCTrack> *MCtracks = nullptr;
        std::vector<o2::itsmft::Hit> *ITSHits = nullptr;

        // ITS tracks
        std::vector<o2::its::TrackITS> *ITStracks = nullptr;
        std::vector<o2::dataformats::TrackTPCITS> *ITSTPCtracks = nullptr;

        // Labels
        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labABvec = nullptr;

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
        treeITSTPC->SetBranchAddress("MatchABMCTruth", &labABvec);

        treeITS->SetBranchAddress("ITSTrack", &ITStracks);
        treeITS->SetBranchAddress("ITSTrackClusIdx", &ITSTrackClusIdx);
        treeITS->SetBranchAddress("ITSTrackClusIdx", &ITSTrackClusIdx);
        treeITSclus->SetBranchAddress("ITSClusterComp", &ITSclus);
        treeITSclus->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);

        // load geometry
        auto gman = o2::its::GeometryTGeo::Instance();
        gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));

        // fill MC matrix
        std::vector<std::vector<o2::MCTrack>> mcTracksMatrix;
        std::vector<std::vector<int>> tpcTracksMatrix; // -1 for no TPC track , 1 for TPC track

        auto nev = treeMCTracks->GetEntriesFast();

        mcTracksMatrix.resize(nev);
        tpcTracksMatrix.resize(nev);

        for (int n = 0; n < nev; n++)
        { // loop over MC events
            treeMCTracks->GetEvent(n);
            mcTracksMatrix[n].resize(MCtracks->size());
            tpcTracksMatrix[n].resize(MCtracks->size());

            for (unsigned int mcI{0}; mcI < MCtracks->size(); ++mcI)
            {
                mcTracksMatrix[n][mcI] = MCtracks->at(mcI);
                tpcTracksMatrix[n][mcI] = -1;

                auto &track = MCtracks->at(mcI);
                if (std::abs(track.GetPdgCode()) == motherPDG)
                {
                    hXiCounter->Fill(0.5);
                }
            }
        }

        for (int frame = 0; frame < treeITSTPC->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSclus->GetEvent(frame))
                continue;

            // revert cluster/track index
            std::vector<int> clusTrackIdx;
            clusTrackIdx.resize(ITSclus->size());

            // dummy fill w/ -1
            std::fill(clusTrackIdx.begin(), clusTrackIdx.end(), -1);

            for (unsigned int iTrack{0}; iTrack < ITStracks->size(); iTrack++)
            {
                auto &ITStrack = ITStracks->at(iTrack);
                auto firstClus = ITStrack.getFirstClusterEntry();
                auto ncl = ITStrack.getNumberOfClusters();
                for (int icl = 0; icl < ncl; icl++)
                {
                    clusTrackIdx[ITSTrackClusIdx->at(firstClus + icl)] = iTrack;
                }
            }

            // looping on clusters
            for (unsigned int iClus{0}; iClus < ITSclus->size(); iClus++)
            {
                auto &labCls = (clusLabArr->getLabels(iClus))[0];
                if (labCls.isFake())
                    continue;

                auto &clus = ITSclus->at(iClus);
                auto trackIdx = clusTrackIdx[iClus];
                auto layer = gman->getLayer(clus.getSensorID());

                auto clRefPr = matchCascDauToMC(mcTracksMatrix, labCls, firstV0dauPDG);
                auto clRefPi = matchCascDauToMC(mcTracksMatrix, labCls, secondV0dauPDG);
                auto clRefBach = matchBachToMC(mcTracksMatrix, labCls);

                std::array<int, 2> clRef;
                auto pdgRef = checkCascRef(clRefPr, clRefPi, clRefBach, clRef);
                if (pdgRef == -1)
                    continue;

                auto mcTrackXi = mcTracksMatrix[clRef[0]][clRef[1]];
                auto decL = pdgRef == kBach ? calcDecLength(&mcTracksMatrix[clRef[0]], mcTrackXi, bachPDG) : calcDecLengthV0(&mcTracksMatrix[clRef[0]], mcTrackXi, firstV0dauPDG);
                auto &hToFill = mapClus[pdgRef];

                // looking for its ITS track
                if (trackIdx != -1)
                {
                    auto &ITSlab = labITSvec->at(trackIdx);
                    ITSlab.isFake() ? hToFill->Fill(kFake, layer) : hToFill->Fill(kTracked, layer);
                }
                else
                {
                    hToFill->Fill(kFree, layer);
                }
            }

            // Study afterburner
            auto abSize = labABvec->size();
            for (unsigned int iTrack{0}; iTrack < labITSTPCvec->size(); iTrack++)
            {
                bool isAb = iTrack > labITSTPCvec->size() - abSize;
                auto &labTot = labITSTPCvec->at(iTrack);
                auto cascRef1 = matchCascDauToMC(mcTracksMatrix, labTot, firstV0dauPDG);
                auto cascRef2 = matchCascDauToMC(mcTracksMatrix, labTot, secondV0dauPDG);
                auto cascRefBach = matchBachToMC(mcTracksMatrix, labTot);

                std::array<int, 2> tpcRef;
                int pdgRef = checkCascRef(cascRef1, cascRef2, cascRefBach, tpcRef);
                if (pdgRef == -1)
                    continue;
            
                auto mcTrackXi = mcTracksMatrix[tpcRef[0]][tpcRef[1]];
                auto decL = pdgRef == kBach ? calcDecLength(&mcTracksMatrix[tpcRef[0]], mcTrackXi, bachPDG) : calcDecLengthV0(&mcTracksMatrix[tpcRef[0]], mcTrackXi, firstV0dauPDG);            

                if(labTot.isFake()){
                    isAb ? mapRadABFake[pdgRef]->Fill(decL) : mapRadFake[pdgRef]->Fill(decL);
                }
                else{
                    isAb ? mapRadAB[pdgRef]->Fill(decL) : mapRad[pdgRef]->Fill(decL);
                }

            }
        }
    }

    auto outFile = TFile("check_its_tpc_cascade.root", "recreate");
    outFile.cd();

    for(auto &pdg : {kPr, kPi, kBach})
    {
        mapClus[pdg]->Write();
        mapRad[pdg]->Write();
        mapRadAB[pdg]->Write();
        mapRadFake[pdg]->Write();
        mapRadABFake[pdg]->Write();
    }

    outFile.Close();
}

