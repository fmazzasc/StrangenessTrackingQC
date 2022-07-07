#if !defined(CLING) || defined(ROOTCLING)
#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "DataFormatsITS/TrackITS.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "TSystemDirectory.h"
#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using Vec3 = ROOT::Math::SVector<double, 3>;

const int motherPDG = 1010010030;
const int firstDaughterPDG = 1000020030;
const int secondDaughterPDG = 211;

// const int motherPDG = 3122;
// const int firstDaughterPDG = 2212;
// const int secondDaughterPDG = 211;

// const int motherPDG = 310;
// const int firstDaughterPDG = 211;
// const int secondDaughterPDG = -211;

o2::its::TrackITS *getITSTrack(int motherEvID, int motherTrackID, TTree *ITStree, std::vector<o2::MCCompLabel> *ITSlabel, std::vector<o2::its::TrackITS> *ITStrack);
void doMatching(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, TTree *treeDetectors, std::vector<o2::MCCompLabel> *labDetectors, TH1D *histo);
double calcMass(const V0 &v0, double dauMass[2], int dauCharges[2]);
double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);
double calcRadius(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);
double calcV0alpha(const V0 &v0);
bool checkV0Decay(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int firstDauPDG, int secDauPDG);
double calcMass(const V0 &v0, double dauMass[2], int dauCharges[2]);

void v0Study()
{
    // std::string path = "/data/fmazzasc/its_data/sim/lambdas/";
    std::string path = "/data/fmazzasc/its_data/sim/hyp_ab_new/";

    // std::string path = "alice/run_sim/";

    double bins[2] = {2.96, 3.04};
    double motherMass = 2.99131;
    double dauMass[2] = {2.80839160743, 0.13957};
    int dauCharges[2] = {2, 1};

    if (std::abs(motherPDG) == 3122)
    {
        motherMass = 1.115683;
        bins[0] = 1.0;
        bins[1] = 1.2;
        dauMass[0] = 0.938272;
        dauMass[1] = 0.13957;
        dauCharges[0] = 1;
        dauCharges[1] = 1;
    }
    if (motherPDG == 310)
    {
        motherMass = 0.493677;
        bins[0] = 0.4;
        bins[1] = 0.6;
        dauMass[0] = 0.13957;
        dauMass[1] = 0.13957;
        dauCharges[0] = 1;
        dauCharges[1] = 1;
    }

    int injectedParticles = 0;

    std::vector<TH1D *> hists(5);
    hists[0] = new TH1D("recoPDGits", "Reconstructed ITS PDG;;Efficiency", 3, 0, 3);
    hists[1] = new TH1D("recoPDGtpc", "Reconstructed TPC PDG;;Efficiency", 3, 0, 3);
    hists[2] = new TH1D("recoPDGitsTPC", "Reconstructed ITS-TPC PDG;;Efficiency", 3, 0, 3);
    hists[3] = new TH1D("recoPDGtpcTOF", "Reconstructed TPC-TOF PDG;;Efficiency", 3, 0, 3);
    hists[4] = new TH1D("recoPDGitsTPCTOF", "Reconstructed ITS-TPC-TOF PDG;;Efficiency", 3, 0, 3);

    TH1D *histInvMass = new TH1D("V0 invariant mass", "; V0 Mass (GeV/c^{2}); Counts", 30, bins[0], bins[1]);

    TH1D *histGenRadius = new TH1D("Gen Radius", "; Gen Radius (cm); Counts", 300, 0, 90);
    TH1D *histRecRadius = new TH1D("Rec Radius", "; Rec Radius (cm); Counts", 300, 0, 90);

    TH1D *histGenDecLength = new TH1D("Gen Dec Length", "; Gen Dec Length (cm); Counts", 300, 0, 90);
    TH1D *histGenLifetime = new TH1D("Gen ct", "; Gen ct (cm); Counts", 300, 0, 90);
    TH1D *histRecDecLength = new TH1D("Rec Dec Length", "; Gen Dec Length (cm); Counts", 300, 0, 90);

    TH1D *histGeneratedV0s = new TH1D("# of generated V0s", ";; Counts", 1, 0, 1);

    TH2D *histV0radiusRes = new TH2D("V0 radius resolution", "; Gen Radius (cm); Gen - Rec / Gen; Counts", 400, 0, 90, 20, -1, 1);

    TH1D *histITSHits = new TH1D("V0 candidate ITS hits", " Number of ITS hits; #Hits; Counts", 7, 0.5, 7.5);
    TH1D *histITScounter = new TH1D("V0 candidate ITS hits and tracks counter", ";; Counts/(# of V0s)", 3, 0, 3);

    // Looping over all the directories inside the path
    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    std::vector<TString> kine_files;
    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (file.substr(0, 2) == "tf")
        {
            dirs.push_back(file);
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

        // Files
        auto fMCTracks = TFile::Open((TString(path + dir + "/") + kine_file));
        auto fSecondaries = TFile::Open((path + dir + "/o2_secondary_vertex.root").data());
        auto fITSTPC = TFile::Open((path + dir + "/o2match_itstpc.root").data());
        auto fTPCTOF = TFile::Open((path + dir + "/o2match_tof_tpc.root").data());
        auto fITSTPCTOF = TFile::Open((path + dir + "/o2match_tof_itstpc.root").data());

        auto fITS = TFile::Open((path + dir + "/o2trac_its.root").data());
        auto fTPC = TFile::Open((path + dir + "/tpctracks.root").data());

        auto fITSclus = TFile::Open((path + dir + "/o2clus_its.root").data());

        if(!fMCTracks || !fSecondaries || !fITSTPC || !fTPCTOF || !fITSTPCTOF || !fITS || !fTPC || !fITSclus)
        {
            LOG(error) << "Could not open one of the files";
            continue;
        }

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeSecondaries = (TTree *)fSecondaries->Get("o2sim");
        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeTPC = (TTree *)fTPC->Get("tpcrec");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
        auto treeITSTPCTOF = (TTree *)fITSTPCTOF->Get("matchTOF");
        auto treeTPCTOF = (TTree *)fTPCTOF->Get("matchTOF");

        // Tracks
        std::vector<o2::MCTrack> *MCtracks = nullptr;
        std::vector<V0> *v0vec = nullptr;
        std::vector<o2::its::TrackITS> *ITStracks = nullptr;
        std::vector<o2::itsmft::ROFRecord> *rofArr = nullptr;

        // Labels
        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTOFvec = nullptr;

        treeSecondaries->SetBranchAddress("V0s", &v0vec);
        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);

        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
        treeTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTOFvec);
        treeITSTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTOFvec);
        treeITS->SetBranchAddress("ITSTrack", &ITStracks);
        treeITS->SetBranchAddress("ITSTracksROF", &rofArr);

        std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}, {"TPC", labTPCvec}, {"ITS-TPC", labITSTPCvec}, {"TPC-TOF", labTPCTOFvec}, {"ITS-TPC-TOF", labITSTPCTOFvec}};

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
                if (std::abs(MCtracks->at(mcI).GetPdgCode()) == motherPDG && MCtracks->at(mcI).isPrimary())
                {
                    auto &mcTrack = mcTracksMatrix[n][mcI];
                    if(!checkV0Decay(MCtracks, mcTrack, firstDaughterPDG, secondDaughterPDG))
                    {
                        continue;
                    }
                    injectedParticles++;
                    histGeneratedV0s->Fill(0.5);
                    double L = calcDecLength(MCtracks, mcTrack, firstDaughterPDG);
                    histGenDecLength->Fill(L);
                    histGenLifetime->Fill(L * motherMass / mcTrack.GetP());
                    histGenRadius->Fill(calcRadius(MCtracks, mcTrack, firstDaughterPDG));
                }
            }
        }

        doMatching(mcTracksMatrix, treeITS, labITSvec, hists[0]);
        doMatching(mcTracksMatrix, treeTPC, labTPCvec, hists[1]);
        doMatching(mcTracksMatrix, treeITSTPC, labITSTPCvec, hists[2]);
        doMatching(mcTracksMatrix, treeTPCTOF, labTPCTOFvec, hists[3]);
        doMatching(mcTracksMatrix, treeITSTPCTOF, labITSTPCTOFvec, hists[4]);

        treeSecondaries->GetEntry();
        treeITS->GetEntry();
        treeTPC->GetEntry();
        treeITSTPC->GetEntry();
        treeTPCTOF->GetEntry();
        treeITSTPCTOF->GetEntry();

        for (auto &v0 : *v0vec)
        {
            std::vector<int> motherIDvec;
            std::vector<int> daughterIDvec;
            std::vector<int> evIDvec;

            for (int iV0 = 0; iV0 < 2; iV0++)
            {

                if (map[v0.getProngID(iV0).getSourceName()])
                {

                    auto labTrackType = map[v0.getProngID(iV0).getSourceName()];
                    auto lab = labTrackType->at(v0.getProngID(iV0).getIndex());

                    int trackID, evID, srcID;
                    bool fake;
                    lab.get(trackID, evID, srcID, fake);
                    if (!lab.isNoise() && lab.isValid() && lab.isCorrect())
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

            auto motherTrack = mcTracksMatrix[evIDvec[0]][motherIDvec[0]];

            if (!motherTrack.isPrimary() || std::abs(motherTrack.GetPdgCode()) != motherPDG)
                continue;

            auto genRad = calcRadius(&mcTracksMatrix[evIDvec[0]], motherTrack, firstDaughterPDG);

            if (calcV0alpha(v0) < 0)
            {
                std::swap(dauMass[0], dauMass[1]);
                std::swap(dauCharges[0], dauCharges[1]);
            }

            histInvMass->Fill(calcMass(v0, dauMass, dauCharges));
            histRecDecLength->Fill(TMath::Sqrt(v0.calcR2() + v0.getZ() * v0.getZ()));

            auto recRad = TMath::Sqrt(v0.calcR2());
            histRecRadius->Fill(recRad);
            histV0radiusRes->Fill(genRad, (recRad - genRad) / genRad);

            counter++;
            std::cout << "---------------------------------" << std::endl;
            std::cout << "Counter: " << counter << std::endl;
            std::cout << evIDvec[0] << ", " << motherIDvec[0] << ", " << motherIDvec[1] << std::endl;
            std::cout << "Common mother found, PDG: " << mcTracksMatrix[evIDvec[0]][motherIDvec[0]].GetPdgCode() << std::endl;
            std::cout << "Daughter 0, PDG: " << pdg0 << ", Pt: " << mcTracksMatrix[evIDvec[0]][daughterIDvec[0]].GetPt() << std::endl;
            std::cout << "Daughter 0, Rec Pt: " << v0.getProng(0).getPt() << ", Track type: " << v0.getProngID(0).getSourceName() << std::endl;
            std::cout << "Daughter 1, PDG: " << pdg1 << ", Pt: " << mcTracksMatrix[evIDvec[0]][daughterIDvec[1]].GetPt() << std::endl;
            std::cout << "Daughter 1, Rec Pt: " << v0.getProng(1).getPt() << ", Track type: " << v0.getProngID(1).getSourceName() << std::endl;

            if (motherTrack.leftTrace(0))
            {
                histITScounter->Fill(0);
                std::cout << "ITS sees mother hits! " << std::endl;
                o2::its::TrackITS *motherITStrack = getITSTrack(evIDvec[0], motherIDvec[0], treeITS, labITSvec, ITStracks);

                if (motherITStrack != nullptr)
                {
                    histITSHits->Fill(motherITStrack->getNClusters());
                    motherITStrack->getNFakeClusters() > 0 ? histITScounter->Fill(2) : histITScounter->Fill(1);
                }
            }
        }
    }
    auto outFile = TFile("v0study.root", "recreate");

    const char *labHits[3] = {"Has at least 1 ITS hit", "Has ITS track w/o fake", "Has ITS track w/ fake"};
    for (auto iLab{0}; iLab < 3; ++iLab)
    {
        histITScounter->SetBinContent(iLab + 1, histITScounter->GetBinContent(iLab + 1) / counter);
        histITScounter->GetXaxis()->SetBinLabel(iLab + 1, labHits[iLab]);
    }

    const char *labels[3] = {std::to_string(motherPDG).data(), std::to_string(firstDaughterPDG).data(), std::to_string(secondDaughterPDG).data()};
    for (auto iH{0}; iH < 5; ++iH)
    {
        for (auto iLab{0}; iLab < 3; ++iLab)
        {
            hists[iH]->GetXaxis()->SetBinLabel(iLab + 1, labels[iLab]);
        }
        for (int iBin = 1; iBin < hists[iH]->GetNbinsX() + 1; iBin++)
        {
            hists[iH]->SetBinContent(iBin, hists[iH]->GetBinContent(iBin) / injectedParticles);
        }
        hists[iH]->Write();
    }

    TH1D *histoEffvsRadius = (TH1D *)histRecRadius->Clone("histoEffvsRadius");
    histoEffvsRadius->Divide(histGenRadius);
    histoEffvsRadius->GetYaxis()->SetTitle("Efficiency");
    histoEffvsRadius->GetXaxis()->SetTitle("V0 Radius (cm)");
    histoEffvsRadius->Write();

    histInvMass->Write();
    histGenLifetime->Write();
    histRecRadius->Write();
    histRecDecLength->Write();
    histGenDecLength->Write();
    histGenRadius->Write();
    histV0radiusRes->Write();
    histITSHits->Write();
    histITScounter->Write();
    histGeneratedV0s->Write();
    outFile.Close();
}

o2::its::TrackITS *getITSTrack(int motherEvID, int motherTrackID, TTree *ITStree, std::vector<o2::MCCompLabel> *ITSlabel, std::vector<o2::its::TrackITS> *ITStrack)
{
    o2::its::TrackITS *motherTrack{nullptr};

    for (int frame = 0; frame < ITStree->GetEntriesFast(); frame++)
    {
        if (!ITStree->GetEvent(frame) || !ITStree->GetEvent(frame))
            continue;
        if (!ITStree->GetEvent(frame))
        {
            continue;
        }
        for (unsigned int iTrack{0}; iTrack < ITSlabel->size(); ++iTrack)
        {
            auto lab = ITSlabel->at(iTrack);
            int trackID, evID, srcID;
            bool fake;
            lab.get(trackID, evID, srcID, fake);
            if (!lab.isNoise() && lab.isValid())
            {
                if (evID == motherEvID and trackID == motherTrackID)
                {
                    std::cout << "Matching indexes: " << evID << "   " << motherTrackID << std::endl;
                    motherTrack = &ITStrack->at(iTrack);
                    std::cout << "ITS sees mother track! " << std::endl;
                    return motherTrack;
                }
            }
        }
    }
    return motherTrack;
};

void doMatching(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, TTree *treeDetectors,
                std::vector<o2::MCCompLabel> *labDetectors,
                TH1D *histo)
{
    for (int frame = 0; frame < treeDetectors->GetEntriesFast(); frame++)
    {
        if (!treeDetectors->GetEvent(frame))
        {
            continue;
        }
        for (unsigned int iTrack{0}; iTrack < labDetectors->size(); ++iTrack)
        {

            auto lab = labDetectors->at(iTrack);
            int trackID, evID, srcID;
            bool fake;
            lab.get(trackID, evID, srcID, fake);
            if (!lab.isNoise() && lab.isValid() && lab.isCorrect() && !fake)
            {
                auto motherID = mcTracksMatrix[evID][trackID].getMotherTrackId();
                auto dauPDG = std::abs(mcTracksMatrix[evID][trackID].GetPdgCode());
                auto momPDG = std::abs(mcTracksMatrix[evID][motherID].GetPdgCode());
                if (dauPDG == motherPDG || momPDG == motherPDG)
                {
                    int ent = -1;
                    switch (dauPDG)
                    {
                    case secondDaughterPDG:
                        ent = 2;
                        break;
                    case firstDaughterPDG:
                        ent = 1;
                        break;
                    case motherPDG:
                        ent = 0;
                        break;
                    }
                    histo->Fill(ent);
                }
            }
        }
    }
};

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

double calcRadius(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
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

bool checkV0Decay(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int firstDauPDG, int secDauPDG)
{
    bool dau1 = false;
    bool dau2 = false;
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    if(idStart == -1 || idStop == -1)
    {
        return false;
    }
    for (auto iD{idStart}; iD <= idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        // LOG(info) << "Dau PDG: " << dauTrack.GetPdgCode();
        if (std::abs(dauTrack.GetPdgCode()) == firstDauPDG)
        {
            dau1 = true;
        }
        if (std::abs(dauTrack.GetPdgCode()) == secDauPDG)
        {
            dau2 = true;
        }
    }
    return dau1 && dau2;
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
                                 (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) +
                             (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ()) *
                                 (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ());
            return sqrt(decLength);
        }
    }
    return -1;
}

double calcMass(const V0 &v0, double dauMass[2], int dauCharges[2])
{
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
