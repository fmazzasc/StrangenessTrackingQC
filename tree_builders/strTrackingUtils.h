

#if !defined(CLING) || defined(ROOTCLING)
#include "ITSMFTSimulation/Hit.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/IOUtils.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"

#include "CommonDataFormat/RangeReference.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TSystemDirectory.h"
#include "TTree.h"
#include <TLorentzVector.h>
#include <gsl/gsl>

#endif

using MCTrack = o2::MCTrack;
using CompClusterExt = o2::itsmft::CompClusterExt;
using ITSCluster = o2::BaseCluster<float>;

class strTrackingUtils
{
public:
    ///  Default constructor
    strTrackingUtils(int motherPdg, int firstDauPdg, int secondDauPdg) : mMotherPDG(motherPdg), mFirstDauPDG(firstDauPdg), mSecondDauPDG(secondDauPdg)
    {
        setMotherMass();
        mV0PDG = mFirstDauPDG;    // first daughter is always the V0 (casc only)
        mBachPDG = mSecondDauPDG; // second daughter is always the bachelor (casc only)
    }

    ///  Destructor
    ~strTrackingUtils() = default;

    int mMotherPDG = -1;
    int mFirstDauPDG = -1;
    int mSecondDauPDG = -1;
    int mV0PDG = -1;
    int mBachPDG = -1;

    float mMotherMass = -1;

    void setMotherMass()
    {
        mMotherMass = -1;
        auto part = TDatabasePDG::Instance()->GetParticle(mMotherPDG);
        if (part)
            mMotherMass = part->Mass();

        if (std::abs(mMotherPDG) == 1010010030) // hypertriton
            mMotherMass = 2.99131;
    }

    std::array<int, 2> matchCascDauToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix,
                                        o2::MCCompLabel &lab, int dauPDG)
    {
        std::array<int, 2> outArray{-1, -1};

        int trackID, evID, srcID;
        bool fake;
        lab.get(trackID, evID, srcID, fake);
        if (lab.isValid())
        {
            auto motherID = mcTracksMatrix[evID][trackID].getMotherTrackId();
            if (motherID >= 0 &&
                std::abs(mcTracksMatrix[evID][trackID].GetPdgCode()) == dauPDG)
            {
                auto v0MomPDG = mcTracksMatrix[evID][motherID].GetPdgCode();
                auto v0MomID = mcTracksMatrix[evID][motherID].getMotherTrackId();
                if (std::abs(v0MomPDG) == mV0PDG && v0MomID >= 0)
                {
                    auto cascPDG = mcTracksMatrix[evID][v0MomID].GetPdgCode();
                    if (std::abs(cascPDG) == mMotherPDG)
                    {
                        outArray[0] = evID;
                        outArray[1] = v0MomID;
                    }
                }
            }
        }

        return outArray;
    }

    std::array<int, 2> matchBachToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel &lab)
    {
        std::array<int, 2> outArray{-1, -1};

        int trackID, evID, srcID;
        bool fake;
        lab.get(trackID, evID, srcID, fake);
        if (lab.isValid())
        {
            auto motherID = mcTracksMatrix[evID][trackID].getMotherTrackId();
            if (motherID >= 0)
            {
                auto cascPDG = mcTracksMatrix[evID][motherID].GetPdgCode();
                if (std::abs(cascPDG) == mMotherPDG)
                {
                    outArray[0] = evID;
                    outArray[1] = motherID;
                }
            }
        }

        return outArray;
    }

    double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack)
    {
        auto idStart = motherTrack.getFirstDaughterTrackId();
        auto idStop = motherTrack.getLastDaughterTrackId();

        if (idStart == -1 || idStop == -1)
            return -1;
        for (auto iD{idStart}; iD <= idStop; ++iD)
        {
            auto dauTrack = MCTracks->at(iD);
            if (std::abs(dauTrack.GetPdgCode()) == mFirstDauPDG)
            {
                auto decLength = (dauTrack.GetStartVertexCoordinatesX() -
                                  motherTrack.GetStartVertexCoordinatesX()) *
                                     (dauTrack.GetStartVertexCoordinatesX() -
                                      motherTrack.GetStartVertexCoordinatesX()) +
                                 (dauTrack.GetStartVertexCoordinatesY() -
                                  motherTrack.GetStartVertexCoordinatesY()) *
                                     (dauTrack.GetStartVertexCoordinatesY() -
                                      motherTrack.GetStartVertexCoordinatesY());
                return sqrt(decLength);
            }
        }
        return -1;
    }

    double calcR2(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack)
    {
        auto idStart = motherTrack.getFirstDaughterTrackId();
        auto idStop = motherTrack.getLastDaughterTrackId();

        if (idStart == -1 || idStop == -1)
            return -1;
        for (auto iD{idStart}; iD <= idStop; ++iD)
        {
            auto dauTrack = MCTracks->at(iD);
            if (std::abs(dauTrack.GetPdgCode()) == mFirstDauPDG)
            {
                auto decLength = (dauTrack.GetStartVertexCoordinatesX() -
                                  motherTrack.GetStartVertexCoordinatesX()) *
                                     (dauTrack.GetStartVertexCoordinatesX() -
                                      motherTrack.GetStartVertexCoordinatesX()) +
                                 (dauTrack.GetStartVertexCoordinatesY() -
                                  motherTrack.GetStartVertexCoordinatesY()) *
                                     (dauTrack.GetStartVertexCoordinatesY() -
                                      motherTrack.GetStartVertexCoordinatesY());
                return decLength;
            }
        }
        return -1;
    }

    double calcLifetime(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack)
    {
        auto idStart = motherTrack.getFirstDaughterTrackId();
        auto idStop = motherTrack.getLastDaughterTrackId();

        if (idStart == -1 || idStop == -1)
            return -1;
        for (auto iD{idStart}; iD <= idStop; ++iD)
        {
            auto dauTrack = MCTracks->at(iD);
            if (std::abs(dauTrack.GetPdgCode()) == mFirstDauPDG)
            {
                auto decLength = (dauTrack.GetStartVertexCoordinatesX() -
                                  motherTrack.GetStartVertexCoordinatesX()) *
                                     (dauTrack.GetStartVertexCoordinatesX() -
                                      motherTrack.GetStartVertexCoordinatesX()) +
                                 (dauTrack.GetStartVertexCoordinatesY() -
                                  motherTrack.GetStartVertexCoordinatesY()) *
                                     (dauTrack.GetStartVertexCoordinatesY() -
                                      motherTrack.GetStartVertexCoordinatesY()) +
                                 (dauTrack.GetStartVertexCoordinatesZ() -
                                  motherTrack.GetStartVertexCoordinatesZ()) *
                                     (dauTrack.GetStartVertexCoordinatesZ() -
                                      motherTrack.GetStartVertexCoordinatesZ());
                return sqrt(decLength) * mMotherMass / motherTrack.GetP();
            }
        }
        return -1;
    }

    double calcDecLengthV0Casc(std::vector<MCTrack> *MCTracks,
                               const MCTrack &motherTrack, int dauV0PDG)
    {
        auto idStart = motherTrack.getFirstDaughterTrackId();
        auto idStop = motherTrack.getLastDaughterTrackId();

        if (idStart == -1 || idStop == -1)
            return -1;
        for (auto iD{idStart}; iD <= idStop; iD++)
        {
            auto &v0Track = MCTracks->at(iD);

            auto jdStart = v0Track.getFirstDaughterTrackId();
            auto jdStop = v0Track.getLastDaughterTrackId();
            if (std::abs(v0Track.GetPdgCode()) == mV0PDG)
            {
                if (jdStart == -1 || jdStop == -1)
                    return -1;

                for (auto jD{jdStart}; jD <= jdStop; jD++)
                {
                    auto dauTrack = MCTracks->at(jD);
                    if (std::abs(dauTrack.GetPdgCode()) == dauV0PDG)
                    {
                        auto decLength = (dauTrack.GetStartVertexCoordinatesX() -
                                          motherTrack.GetStartVertexCoordinatesX()) *
                                             (dauTrack.GetStartVertexCoordinatesX() -
                                              motherTrack.GetStartVertexCoordinatesX()) +
                                         (dauTrack.GetStartVertexCoordinatesY() -
                                          motherTrack.GetStartVertexCoordinatesY()) *
                                             (dauTrack.GetStartVertexCoordinatesY() -
                                              motherTrack.GetStartVertexCoordinatesY());
                        return sqrt(decLength);
                    }
                }
            }
        }
        return -1;
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

    std::array<int, 2> matchITStracktoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel)

    {
        std::array<int, 2> outArray = {-1, -1};
        int trackID, evID, srcID;
        bool fake;
        ITSlabel.get(trackID, evID, srcID, fake);
        // LOG(info) << "pdg: " << mcTracksMatrix[evID][trackID].GetPdgCode();
        if (ITSlabel.isValid() &&
            std::abs(mcTracksMatrix[evID][trackID].GetPdgCode()) == mMotherPDG)
        {
            outArray = {evID, trackID};
        }

        return outArray;
    }

    std::vector<ITSCluster>
    getTrackClusters(const o2::its::TrackITS &ITStrack,
                     const std::vector<ITSCluster> &ITSClustersArray,
                     std::vector<int> *ITSTrackClusIdx)
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

    // double calcMass(std::vector<o2::track::TrackParCovF> tracks)
    // {
    //     TLorentzVector moth, prong;
    //     std::array<float, 3> p;
    //     for (unsigned int i = 0; i < tracks.size(); i++)
    //     {
    //         auto &track = tracks[i];
    //         auto mass = tracks.size() == 2 ? kFirstDauMasses[i] : kDauMasses[i];
    //         track.getPxPyPzGlo(p);
    //         prong.SetVectM({p[0], p[1], p[2]}, mass);
    //         moth += prong;
    //     }
    //     return moth.M();
    // }
};