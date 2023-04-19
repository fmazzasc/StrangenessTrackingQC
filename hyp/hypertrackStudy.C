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
#include "ReconstructionDataFormats/StrangeTrack.h"

#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using CompClusterExt = o2::itsmft::CompClusterExt;
using ITSCluster = o2::BaseCluster<float>;
using Vec3 = ROOT::Math::SVector<double, 3>;
using StrangeTrack = o2::dataformats::StrangeTrack;

const int motherPDG = 1010010030;
const int firstDaughterPDG = 1000020030;
const int secondDaughterPDG = 211;

double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);
double calcV0alpha(const V0 &v0);
double recomputeV0Pt(const V0 &v0);
double calcMass(const V0 &v0);
double calcMass(std::vector<o2::track::TrackParCovF> tracks);
double calcCtau(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);

std::vector<std::array<int, 2>> matchV0stoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec);
std::array<int, 2> matchV0DautoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, o2::dataformats::V0::GIndex dauID);
std::array<int, 2> matchITStracktoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel);
std::array<int, 2> matchCompLabelToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel);
std::vector<ITSCluster> getTrackClusters(const o2::its::TrackITS &ITStrack, const std::vector<ITSCluster> &ITSClustersArray, std::vector<int> *ITSTrackClusIdx);

void hypertrackStudy()
{
    // Output Histograms
    TH1D *hChi2Sgn = new TH1D("Chi2 Signal", "; #chi^{2}; Counts", 102, -2, 100);
    TH1D *hChi2Bkg = new TH1D("Chi2 Fake assoc", "; #chi^{2}; Counts", 102, -2, 100);

    TH1D *hSigBkg = new TH1D("Hypertracker eff", "; ; Efficiency", 2, 0, 2);
    TH2D *hMChisto = new TH2D("histo mc", ";#it{p}_{T} (GeV/#it{c}); Radius^2 (cm) ; Counts", 20, 1, 10, 30, 1, 900);
    TH2D *hV0histo = new TH2D("histo V0", ";#it{p}_{T} (GeV/#it{c}); Radius^2 (cm) ; Counts", 20, 1, 10, 30, 1, 900);
    TH2D *hHyperhisto = new TH2D("histo hyperV0", ";#it{p}_{T} (GeV/#it{c}); Radius^2 (cm) ; Counts", 20, 1, 10, 30, 1, 900);
    TH1D *hResV0histo = new TH1D("pT resolution before hypertracking", ";(#it{p}_{T}^{gen} - #it{p}_{T}^{rec})/#it{p}_{T}^{gen}; Counts", 20, -0.2, 0.2);
    TH1D *hResHyperhisto = new TH1D("pT resolution after hypertracking", ";(#it{p}_{T}^{gen} - #it{p}_{T}^{rec})/#it{p}_{T}^{gen}; Counts", 20, -0.2, 0.2);
    TH1D *hResV0histoR2 = new TH1D("R2 resolution before hypertracking", ";(R2^{gen} - R2^{rec})/R2^{gen}; Counts", 200, -0.2, 0.2);
    TH1D *hResHyperhistoR2 = new TH1D("R2 resolution after hypertracking", ";(R2^{gen} - R2^{rec})/R2^{gen}; Counts", 200, -0.2, 0.2);
    TH1D *hV0Counter = new TH1D("V0 counter", ";V0 counter; Counts", 1, 0.5, 1.5);

    TH1D *hRecHypCounter = new TH1D("Rec V0 hyp counter", "; ; Counts", 1, 0.5, 1.5);
    TH1D *hHypertrackerStats = new TH1D("hypertracker_stats", "; ; Counts", 3, 0.5, 3.5);
    TH1D *hHyperCounter = new TH1D("Hypertrack counter", ";Hypertrack counter; Counts", 1, 0.5, 1.5);
    TH1D *hFakeAssocCounter = new TH1D("Fake assoc counter", ";Fake assoc counter; Counts", 1, 0.5, 1.5);

    TH1D *hV0InvMass = new TH1D("v0_inv_mass", "; M ; Counts", 120, 2.96, 3.04);
    TH1D *hStrTrackInvMass = new TH1D("strtrack_inv_mass", "; M ; Counts", 120, 2.96, 3.04);

    TH1D *hGenHypRadius = new TH1D("gen_hyp_r", "; Radius (cm); Counts", 40, 0., 40.);
    TH1D *hGenHypCt = new TH1D("gen_hyp_ct", "; #it{c}t_{gen} (cm); Counts", 200, 0., 40.);

    TH1D *hGenHypMom = new TH1D("gen_hyp_pt", "; #it{p}_{T}^{gen}; Counts", 40, 1, 10);

    TH1D *hRecHypRadius = new TH1D("rec_hyp_r", "; Radius (cm); Counts", 40, 0., 40.);
    TH1D *hRecHypMom = new TH1D("rec_hyp_pt", "; #it{p}_{T}^{gen} (GeV/#it{c}); Counts", 40, 1, 10);

    TH1D *hRecHypRadiusTrackab = new TH1D("rec_hyp_r_trackable", "; Radius (cm); Counts", 40, 2., 40.);
    TH1D *hRecHypRadiusFakes = new TH1D("rec_hyp_r_fakes", "; Radius (cm); Counts", 40, 2., 40.);

    TH1D *hRecHypRadiusTracked = new TH1D("rec_hyp_r_tracked", "; Radius (cm); Counts", 40, 2., 40.);

    TH1D *hV0InvMassForRes = new TH1D("v0_inv_mass_res_study", "; M (GeV/c^{2}) ; Counts", 120, 2.96, 3.04);
    TH1D *hTrackedInvMassForRes = new TH1D("strtrack_inv_mass_res_study", "; M (GeV/c^{2}) ; Counts", 120, 2.96, 3.04);

    std::string path = "/data/fmazzasc/its_data/sim/hyp_2_body/";
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
        auto fHyperTracks = TFile::Open((dir + "/o2_strange_tracks.root").data());
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
        auto treeStrangeTracks = (TTree *)fHyperTracks->Get("o2sim");
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
        std::vector<float> *hyperChi2vec = nullptr;
        std::vector<o2::strangeness_tracking::ClusAttachments> *nAttachments = nullptr;

        // Secondary Vertices
        std::vector<V0> *v0vec = nullptr;
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

        treeSecondaries->SetBranchAddress("V0s", &v0vec);
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
                    auto &mcTrack = mcTracksMatrix[n][mcI];
                    hMChisto->Fill(mcTrack.GetPt(), calcDecLength(MCtracks, mcTrack, firstDaughterPDG));
                    hGenHypMom->Fill(mcTrack.GetPt());
                    hGenHypRadius->Fill(calcDecLength(MCtracks, mcTrack, firstDaughterPDG));
                    hGenHypCt->Fill(calcCtau(MCtracks, mcTrack, firstDaughterPDG));
                }
            }
        }

        // Starting matching  loop
        for (int frame = 0; frame < treeITS->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeSecondaries->GetEvent(frame) ||
                !treeStrangeTracks->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeITSclus->GetEvent(frame) || !treeTPCTRD->GetEvent(frame))
                continue;

            std::vector<std::array<int, 2>> V0sMCref = matchV0stoMC(mcTracksMatrix, map, v0vec);

            std::vector<bool> matchedHypertrackVec;
            matchedHypertrackVec.resize(strangeTrackVec->size());
            for (auto &&isMatched : matchedHypertrackVec)
                isMatched = false;

            for (unsigned int iV0vec = 0; iV0vec < v0vec->size(); iV0vec++)
            {
                hV0Counter->Fill(1);

                auto &v0MCref = V0sMCref[iV0vec];
                auto &v0 = (*v0vec)[iV0vec];

                hV0InvMass->Fill(calcMass(v0));

                if (v0MCref[0] == -1 || v0MCref[1] == -1)
                    continue;
                hRecHypCounter->Fill(1);

                auto &mcTrack = mcTracksMatrix[v0MCref[0]][v0MCref[1]];
                hV0histo->Fill(recomputeV0Pt(v0), v0.calcR2());
                hRecHypMom->Fill(recomputeV0Pt(v0));
                hRecHypRadius->Fill(sqrt(v0.calcR2()));

                // Matching ITS tracks to MC tracks and V0
                std::array<int, 2> ITSref = {-1, 1};
                o2::its::TrackITS ITStrack;
                int trackIdx{-1};
                std::array<std::array<int, 2>, 7> clsRef;

                int iTrack = -1;
                bool isMatched = false;
                o2::MCCompLabel labITSmatched;

                for (unsigned int iITStrack = 0; iITStrack < ITStracks->size(); iITStrack++)
                {
                    auto &labITS = (*labITSvec)[iITStrack];
                    auto &trackIdx = (*ITSTrackClusIdx)[iITStrack];

                    ITSref = matchITStracktoMC(mcTracksMatrix, labITS);
                    ITStrack = (*ITStracks)[iITStrack];

                    if (ITSref[0] == V0sMCref[iV0vec][0] && ITSref[1] == V0sMCref[iV0vec][1])
                    {

                        auto firstClus = ITStrack.getFirstClusterEntry();
                        auto ncl = ITStrack.getNumberOfClusters();
                        for (int icl = 0; icl < ncl; icl++)
                        {
                            auto &labCls = (clusLabArr->getLabels(ITSTrackClusIdx->at(firstClus + icl)))[0];
                            auto &clus = (*ITSclus)[(*ITSTrackClusIdx)[firstClus + icl]];
                            auto layer = gman->getLayer(clus.getSensorID());
                            clsRef[layer] = matchCompLabelToMC(mcTracksMatrix, labCls);
                        }
                        hHypertrackerStats->Fill(1);
                        isMatched = true;
                        iTrack = iITStrack;
                        labITSmatched = labITS;
                        break;
                    }
                }

                if (!isMatched)
                    continue;

                hRecHypRadiusTrackab->Fill(sqrt(v0.calcR2()));
                if (labITSmatched.isFake())
                    hRecHypRadiusFakes->Fill(sqrt(v0.calcR2()));

                LOG(info) << V0sMCref[iV0vec][0] << "  " << V0sMCref[iV0vec][1];

                // Matching hypertracks to MC tracks, V0s and ITS tracks
                bool isHypertracked = false;
                for (unsigned int iHyperVec = 0; iHyperVec < strangeTrackVec->size(); iHyperVec++)
                {
                    auto &strangeTrack = strangeTrackVec->at(iHyperVec);
                    if (strangeTrack.isCascade)
                        continue;
                    auto &hyperV0 = v0vec->at(strangeDecRefVec->at(iHyperVec));
                    auto &matchChi2 = strangeTrack.mMatchChi2;
                    auto &hyperITSref = strangeITSrefVec->at(iHyperVec);

                    if (hyperV0.getProngID(0) == v0.getProngID(0) && hyperV0.getProngID(1) == v0.getProngID(1) && hyperITSref == iTrack)
                    {

                        LOG(info) << "++++++++++++++++++++++++";
                        LOG(info) << V0sMCref[iV0vec][0] << "  " << V0sMCref[iV0vec][1] << ", iHypervec: " << iHyperVec << ", size" << strangeTrackVec->size();
                        LOG(info) << "Hypertrack found!: ITS track ref: " << hyperITSref;
                        isHypertracked = true;
                        matchedHypertrackVec[iHyperVec] = true;
                        hRecHypRadiusTracked->Fill(sqrt(v0.calcR2()));
                        hHyperCounter->Fill(1);
                        hHypertrackerStats->Fill(2);
                        hHyperhisto->Fill(strangeTrack.mMother.getPt(), hyperV0.calcR2());
                        hResHyperhisto->Fill((strangeTrack.mMother.getPt() - mcTrack.GetPt()) / mcTrack.GetPt());
                        hResV0histo->Fill((recomputeV0Pt(v0) - mcTrack.GetPt()) / mcTrack.GetPt());
                        hResV0histoR2->Fill((sqrt(v0.calcR2()) - calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG)) / calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG));
                        auto hyperR2 = strangeTrack.decayVtx[0] * strangeTrack.decayVtx[0] + strangeTrack.decayVtx[1] * strangeTrack.decayVtx[1];

                        hResHyperhistoR2->Fill((sqrt(hyperR2) - calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG)) / calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG));
                        hChi2Sgn->Fill(matchChi2);
                        break;
                    }
                }

                if (!isHypertracked)
                {
                    LOG(info) << "------------------";
                    LOG(info) << "No hypertrack found, but ITS track found";
                    LOG(info) << "processing frame " << frame << ", dir: " << dir;
                    LOG(info) << "V0 pos: " << v0.getProngID(0) << " V0 neg: " << v0.getProngID(1) << ", ITS track ref: " << iTrack;
                    LOG(info) << "V0 Eta: " << v0.getEta() << " V0 phi" << v0.getPhi() << " ITS eta: " << ITStrack.getEta() << " ITS phi: " << ITStrack.getPhi();

                    LOG(info) << "V0 pos ref: " << v0.getProngID(0) << ", neg ref: " << v0.getProngID(1);
                    LOG(info) << "ITS ref: " << iTrack;

                    LOG(info) << "Number of hits: " << ITStrack.getNClusters();
                    LOG(info) << "ITS Track ref: " << ITSref[0] << " " << ITSref[1] << " , PDG: " << mcTracksMatrix[ITSref[0]][ITSref[1]].GetPdgCode();
                    LOG(info) << "+++++++";
                    hHypertrackerStats->Fill(3);

                    for (unsigned int i{0}; i < 7; i++)
                    {
                        if (ITStrack.hasHitOnLayer(i))
                        {
                            LOG(info) << "ITS track has hit on layer " << i << ", is fake: " << ITStrack.isFakeOnLayer(i);
                            if (clsRef[i][0] != -1 && clsRef[i][1] != -1)
                                LOG(info) << "Cluster ref: " << clsRef[i][0] << " " << clsRef[i][1] << " , PDG: " << mcTracksMatrix[clsRef[i][0]][clsRef[i][1]].GetPdgCode();
                            if (clsRef[i][0] != ITSref[0])
                            {
                                LOG(info) << "EvID mismatch: " << clsRef[i][0] << " " << ITSref[0];
                                continue;
                            }
                            if (clsRef[i][1] != ITSref[1])
                            {

                                auto motherID = mcTracksMatrix[clsRef[i][0]][clsRef[i][1]].getMotherTrackId();
                                if (motherID != -1)
                                    LOG(info) << "Mother cluster ref: " << clsRef[i][0] << " " << motherID << " , PDG: " << mcTracksMatrix[clsRef[i][0]][motherID].GetPdgCode();
                            }
                        }
                    }
                }
            }

            // LOG(info) << " Start studying fake associations";
            for (unsigned int iHyperVec = 0; iHyperVec < strangeTrackVec->size(); iHyperVec++)
            {

                auto &strangeTrack = strangeTrackVec->at(iHyperVec);
                if (strangeTrack.isCascade)
                {
                    continue;
                }
                auto &hyperChi2 = strangeTrack.mMatchChi2;
                auto &clusAttachments = nAttachments->at(iHyperVec);
                auto &ITStrack = ITStracks->at(strangeITSrefVec->at(iHyperVec));
                auto &ITStrackLab = labITSvec->at(strangeITSrefVec->at(iHyperVec));

                auto clusAttArr = clusAttachments.arr;
                bool isMatched = matchedHypertrackVec[iHyperVec];
                auto &hyperV0 = v0vec->at(strangeDecRefVec->at(iHyperVec));

                std::vector<o2::track::TrackParCovF> dauTracks = {strangeTrack.mDaughterFirst, strangeTrack.mDaughterSecond};

                auto v0PosRef = matchV0DautoMC(mcTracksMatrix, map, hyperV0.getProngID(0));
                auto v0NegRef = matchV0DautoMC(mcTracksMatrix, map, hyperV0.getProngID(1));
                auto ITStrackRef = matchCompLabelToMC(mcTracksMatrix, ITStrackLab);

                if (isMatched)
                {
                    hTrackedInvMassForRes->Fill(calcMass(dauTracks));
                    hStrTrackInvMass->Fill(calcMass(dauTracks));
                    hV0InvMassForRes->Fill(calcMass(hyperV0));
                    continue;
                }

                // LOG(info) << "**************";
                // LOG(info) << "ITS track position: " << strangeITSrefVec->at(iHyperVec);
                // LOG(info) << "V0 pos: " << hyperV0.getProngID(0) << " V0 neg: " << hyperV0.getProngID(1) << " V0pt: " << hyperV0.getPt() << " ITSpt: " << ITStrack.getPt() << ", IS matched:" << isMatched << ", iHyperVec: " << iHyperVec << ", size" << strangeTrackVec->size();
                // LOG(info) << "V0 pos lab: " << v0PosRef[0] << " " << v0PosRef[1] << ", V0 neg lab: " << v0NegRef[0] << " " << v0NegRef[1] << ", ITS ref: " << ITStrackRef[0] << " " << ITStrackRef[1];

                // LOG(info) << "V0 pos lab: " << v0PosRef[0] << " " << v0PosRef[1] << ", V0 neg lab: " << v0NegRef[0] << " " << v0NegRef[1] << ", ITS ref: " << ITStrackRef[0] << " " << ITStrackRef[1];
                if (v0PosRef[0] >= 0 && v0PosRef[1] >= 0 && v0NegRef[0] >= 0 && v0NegRef[1] >= 0)
                {
                    // LOG(info) << "V0 pos mother: " << mcTracksMatrix[v0PosRef[0]][v0PosRef[1]].GetPdgCode();
                    // LOG(info) << "V0 neg mother: " << mcTracksMatrix[v0NegRef[0]][v0NegRef[1]].GetPdgCode();
                }
                else
                {
                    // LOG(info) << "Mother not found!";
                }

                hChi2Bkg->Fill(hyperChi2);
                hFakeAssocCounter->Fill(1);
                auto firstClus = ITStrack.getFirstClusterEntry();
                auto ncl = ITStrack.getNumberOfClusters();
                for (int icl = 0; icl < ncl; icl++)
                {
                    auto &labCls = (clusLabArr->getLabels(ITSTrackClusIdx->at(firstClus + icl)))[0];
                    auto &clus = (*ITSclus)[(*ITSTrackClusIdx)[firstClus + icl]];
                    auto layer = gman->getLayer(clus.getSensorID());
                    std::array<int, 2> clsRef = matchCompLabelToMC(mcTracksMatrix, labCls);
                    // if (clsRef[0] > -1 && clsRef[1] > -1)
                    //     LOG(info) << "Layer: " << layer << "PDG: " << mcTracksMatrix[clsRef[0]][clsRef[1]].GetPdgCode() << ", Attached to: " << clusAttArr[layer];
                    // else
                    //     LOG(info) << "Layer: " << layer << ", No valid cluster ref";
                }
            }
        }
    }

    auto outFile = TFile("hypertrack_study_2.root", "recreate");
    hChi2Sgn->Write();
    hChi2Bkg->Write();

    hV0histo->Write();
    hResV0histo->Write();
    hResV0histoR2->Write();
    hHyperhisto->Write();
    hResHyperhisto->Write();
    hResHyperhistoR2->Write();
    hMChisto->Write();

    hV0Counter->Write();
    hHypertrackerStats->Write();
    hHyperCounter->Write();
    hFakeAssocCounter->Write();
    hRecHypCounter->Write();

    hV0InvMass->Write();
    hStrTrackInvMass->Write();

    auto cv = TCanvas("inv_mass_hyp", "", 1000, 1000);
    hV0InvMass->GetXaxis()->SetTitle("M(GeV/#it{c}^{2})");
    hV0InvMass->GetYaxis()->SetTitle("Normalised Counts");

    hV0InvMass->SetLineWidth(2);

    hStrTrackInvMass->SetLineWidth(2);
    hStrTrackInvMass->SetLineColor(kRed);
    hStrTrackInvMass->DrawNormalized();
    hV0InvMass->DrawNormalized("same");
    hStrTrackInvMass->SetLineColor(kRed);
    auto leg = new TLegend(0.5, 0.5, 0.8, 0.8);
    leg->AddEntry(hV0InvMass, "Before tracking");
    leg->AddEntry(hStrTrackInvMass, "After tracking");
    leg->Draw();
    cv.Write();

    auto cv2 = TCanvas("inv_mass_res_study", "", 1000, 1000);
    hV0InvMassForRes->GetXaxis()->SetTitle("M(GeV/#it{c}^{2})");
    hV0InvMassForRes->GetYaxis()->SetTitle("Normalised Counts");
    hV0InvMassForRes->SetLineWidth(2);
    hTrackedInvMassForRes->SetLineWidth(2);

    hTrackedInvMassForRes->DrawNormalized();
    hV0InvMassForRes->DrawNormalized("same");
    hTrackedInvMassForRes->SetLineColor(kRed);
    auto leg2 = new TLegend(0.5, 0.5, 0.8, 0.8);
    leg2->AddEntry(hV0InvMassForRes, "Before tracking");
    leg2->AddEntry(hTrackedInvMassForRes, "After tracking");
    leg2->Draw();
    cv2.Write();

    hRecHypRadius->Sumw2();
    hRecHypRadius->Divide(hGenHypRadius);
    auto cv3 = TCanvas("reco_efficiency_r", "", 1000, 1000);
    hRecHypRadius->GetXaxis()->SetTitle("Radius (cm)");
    hRecHypRadius->GetYaxis()->SetTitle("Efficiency");
    hRecHypRadius->SetLineColor(kRed);
    hRecHypRadius->Draw();
    cv3.Write();

    hRecHypMom->Write();
    hGenHypMom->Write();

    hRecHypMom->Sumw2();
    hRecHypMom->Divide(hGenHypMom);
    auto cv4 = TCanvas("reco_efficiency_pt", "", 1000, 1000);
    hRecHypMom->GetYaxis()->SetTitle("Efficiency");
    hRecHypMom->SetLineColor(kRed);
    hRecHypMom->Draw("pe");
    cv4.Write();

    auto cv0 = TCanvas("trackable_tracked_r", "", 1000, 1000);
    hRecHypRadiusTracked->GetXaxis()->SetTitle("Radius (cm)");
    hRecHypRadiusTracked->GetYaxis()->SetTitle("Efficiency_{trackable}");
    hRecHypRadiusTracked->SetLineColor(kRed);
    hRecHypRadiusTracked->SetLineWidth(2.);
    hRecHypRadiusTrackab->Draw();
    hRecHypRadiusTracked->Draw("same");
    cv0.Write();

    hRecHypRadiusTracked->Divide(hRecHypRadiusTrackab);
    auto cv1 = TCanvas("trackable_efficiency_r", "", 1000, 1000);
    hRecHypRadiusTracked->GetXaxis()->SetTitle("Radius (cm)");
    hRecHypRadiusTracked->GetYaxis()->SetTitle("Efficiency_{trackable}");
    hRecHypRadiusTracked->SetLineColor(kRed);
    hRecHypRadiusTracked->SetLineWidth(2.);
    hRecHypRadiusTracked->Draw();
    cv1.Write();

    hRecHypRadiusFakes->Write();
    hRecHypRadiusFakes->Sumw2();
    hRecHypRadiusFakes->Divide(hRecHypRadiusTrackab);
    auto cv5 = TCanvas("fakes_efficiency_r", "", 1000, 1000);
    hRecHypRadiusFakes->GetXaxis()->SetTitle("Radius (cm)");
    hRecHypRadiusFakes->GetYaxis()->SetTitle("#fakes / #trackables");
    hRecHypRadiusFakes->SetLineColor(kRed);
    hRecHypRadiusFakes->SetLineWidth(2.);
    hRecHypRadiusFakes->Draw("pe");
    cv5.Write();

    hGenHypCt->Write();
    hRecHypRadiusTrackab->Write();

    outFile.Close();
}
std::vector<std::array<int, 2>> matchV0stoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec)
{
    std::vector<std::array<int, 2>> outArray;
    outArray.resize(v0vec->size());
    int count_V0 = 0;
    for (unsigned int iV0vec = 0; iV0vec < v0vec->size(); iV0vec++)
    {
        std::vector<int> motherIDvec;
        std::vector<int> daughterIDvec;
        std::vector<int> evIDvec;

        outArray[iV0vec] = {-1, -1};
        auto &v0 = (*v0vec)[iV0vec];

        for (unsigned int iV0 = 0; iV0 < 2; iV0++)
        {
            if (map[v0.getProngID(iV0).getSourceName()])
            {
                auto labTrackType = map[v0.getProngID(iV0).getSourceName()];
                auto lab = labTrackType->at(v0.getProngID(iV0).getIndex());

                int trackID, evID, srcID;
                bool fake;
                lab.get(trackID, evID, srcID, fake);
                if (lab.isValid())
                {
                    auto motherID = mcTracksMatrix[evID][trackID].getMotherTrackId();
                    motherIDvec.push_back(motherID);
                    daughterIDvec.push_back(trackID);
                    evIDvec.push_back(evID);
                }
            }
        }

        if (motherIDvec.size() < 2)
            continue;
        if (motherIDvec[0] != motherIDvec[1] || evIDvec[0] != evIDvec[1])
            continue;

        if (motherIDvec[0] < 0 || motherIDvec[0] > 10000)
            continue;

        int pdg0 = mcTracksMatrix[evIDvec[0]][daughterIDvec[0]].GetPdgCode();
        int pdg1 = mcTracksMatrix[evIDvec[0]][daughterIDvec[1]].GetPdgCode();

        if (!(std::abs(pdg0) == firstDaughterPDG && std::abs(pdg1) == secondDaughterPDG) && !(std::abs(pdg0) == secondDaughterPDG && std::abs(pdg1) == firstDaughterPDG))
            continue;

        // std::cout << "Mother PDG: " << mcTracksMatrix[evIDvec[0]][motherIDvec[0]].GetPt() << std::endl;
        outArray[iV0vec] = {evIDvec[0], motherIDvec[0]};
        count_V0++;
    }
    std::cout << "Number of V0s: " << count_V0 << std::endl;
    return outArray;
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

std::array<int, 2> matchV0DautoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, o2::dataformats::V0::GIndex dauID)
{
    std::array<int, 2> outArray{-1, -1};

    if (map[dauID.getSourceName()])
    {
        auto labTrackType = map[dauID.getSourceName()];
        auto lab = labTrackType->at(dauID.getIndex());

        int trackID, evID, srcID;
        bool fake;
        lab.get(trackID, evID, srcID, fake);
        if (lab.isValid())
        {
            auto motherID = mcTracksMatrix[evID][trackID].getMotherTrackId();
            outArray = {evID, motherID};
        }
    }

    return outArray;
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
            return sqrt(decLength);
        }
    }
    return -1;
}

double calcCtau(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
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
                                 (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) +
                             (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ()) *
                                 (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ());
            return sqrt(decLength) * 2.99131 / motherTrack.GetP();
        }
    }
    return -1;
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

double recomputeV0Pt(const V0 &v0)
{
    double alpha = calcV0alpha(v0);
    std::array<float, 3> fPmom, fNmom = {0, 0, 0};
    v0.getProng(0).getPxPyPzGlo(fPmom);
    v0.getProng(1).getPxPyPzGlo(fNmom);

    if (alpha > 0)
        return sqrt((2 * fPmom[0] + fNmom[0]) * (2 * fPmom[0] + fNmom[0]) + (2 * fPmom[1] + fNmom[1]) * (2 * fPmom[1] + fNmom[1]));

    return sqrt((fPmom[0] + 2 * fNmom[0]) * (fPmom[0] + 2 * fNmom[0]) + (fPmom[1] + 2 * fNmom[1]) * (fPmom[1] + 2 * fNmom[1]));
}

double calcMass(const V0 &v0)
{
    std::vector<o2::dataformats::V0::Track> dauTracks = {v0.getProng(0), v0.getProng(1)};
    std::vector<int> dauCharges = {2, 1};
    std::vector<float> dauMasses = {2.80839160743, 0.13957};

    if (calcV0alpha(v0) < 0)
        std::swap(dauCharges[0], dauCharges[1]);
    if (calcV0alpha(v0) < 0)
        std::swap(dauMasses[0], dauMasses[1]);

    TLorentzVector moth, prong;
    std::array<float, 3> p;
    for (int i = 0; i < 2; i++)
    {
        auto &track = dauTracks[i];
        auto &mass = dauMasses[i];
        track.getPxPyPzGlo(p);
        int charge = dauCharges[i];
        prong.SetVectM({charge * p[0], charge * p[1], charge * p[2]}, mass);
        moth += prong;
    }
    return moth.M();
}

double calcMass(std::vector<o2::track::TrackParCovF> tracks)
{
    std::vector<float> dauMasses = {2.80839160743, 0.13957};
    TLorentzVector moth, prong;
    std::array<float, 3> p;
    for (unsigned int i = 0; i < tracks.size(); i++)
    {
        auto &track = tracks[i];
        auto &mass = dauMasses[i];
        track.getPxPyPzGlo(p);
        prong.SetVectM({p[0], p[1], p[2]}, mass);
        moth += prong;
    }
    return moth.M();
}