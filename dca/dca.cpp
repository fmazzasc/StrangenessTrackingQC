// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.



#include "TFile.h"
#include "TNtuple.h"
#include "TH3F.h"
#include "SimulationDataFormat/MCTrack.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "MathUtils/Cartesian.h"
#include "GPUCommonArray.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsITS/TrackITS.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "Steer/MCKinematicsReader.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "SimulationDataFormat/DigitizationContext.h"
#include "StrangenessTracking/StrangenessTracker.h"
#include <vector>
#include <fstream>

using namespace o2;
using namespace o2::framework;

using GIndex = o2::dataformats::VtxTrackIndex;


template <typename T>
void BinLogAxisY(T h) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis* axis = h->GetYaxis();
  int bins = axis->GetNbins();
  double from = axis->GetXmin();
  double to = axis->GetXmax();
  double *newBins = new double[bins + 1];

  newBins[0] = from;
  double factor = std::pow(to / from, 1. / bins);

  for (int i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i - 1];
  }
  axis->Set(bins, newBins);
  delete[] newBins;
}

template <typename TrackT>
std::vector<TrackT> *fetchTracks(const char *filename, const char *treename, const char *branchname)
{
  TFile file(filename, "OPEN");
  auto tree = (TTree *)file.Get(treename);
  auto br = tree->GetBranch(branchname);
  std::vector<TrackT> *tracks = nullptr;
  br->SetAddress(&tracks);
  br->GetEntry(0);
  file.Close();
  return tracks;
}

// add vertices/collisions
void dca()
{
  int tf =2;

  // open the file for vertices
  TFile f(Form("/Users/andrea/alice/data/xi_new/tf%d/o2_primary_vertex.root", tf), "OPEN");

  auto t = (TTree *)f.Get("o2sim");

  const auto grp = o2::parameters::GRPObject::loadFrom(Form("/Users/andrea/alice/data/xi_new/tf%d/o2sim_grp.root",tf)); 
  o2::base::GeometryManager::loadGeometry(Form("/Users/andrea/alice/data/xi_new/tf%d/o2sim_geometry-aligned.root",tf));

  //fetch the tracks (these names are not following any convention!!)
  auto itstracks = fetchTracks<o2::its::TrackITS>(Form("/Users/andrea/alice/data/xi_new/tf%d/o2trac_its.root",tf), "o2sim", "ITSTrack");
  auto strangetracks = fetchTracks<o2::dataformats::StrangeTrack>(Form("/Users/andrea/alice/data/xi_new/tf%d/o2_strange_tracks.root",tf), "o2sim", "StrangeTracks");

  if (t)
  {
    auto br = t->GetBranch("PrimaryVertex");
    std::vector<o2::dataformats::PrimaryVertex> *vertices = nullptr;
    br->SetAddress(&vertices);
    br->GetEntry(0);

    // this referes to actual tracks
    auto indexbr = t->GetBranch("PVTrackIndices");
    std::vector<GIndex> *vertexTrackIDs = nullptr;
    indexbr->SetAddress(&vertexTrackIDs);
    indexbr->GetEntry(0);

    // this makes the connection of vertex to track indices
    auto v2totrackrefbr = t->GetBranch("PV2TrackRefs");
    std::vector<o2::dataformats::VtxTrackRef> *v2trackref = nullptr;
    v2totrackrefbr->SetAddress(&v2trackref);
    v2totrackrefbr->GetEntry(0);


    if (vertices && vertexTrackIDs)
    {

      o2::base::Propagator::initFieldFromGRP(grp);
      auto propagator = o2::base::Propagator::Instance();
      TFile outfile(Form("outputTF%d.root",tf), "RECREATE");
      TH3 *hImpactParameterXY, *hImpactParameterZ;
     
      hImpactParameterXY = new TH3F("hImpactParameterXY", ";Vertex contributors;#it{p}_{T} (GeV/#it{c});Impact parameter XY (cm)", 100, 0.5f, 100.5f, 200, 0.05f, 10.f, 1000, -1.f, 1.f);
      BinLogAxisY(hImpactParameterXY);
      hImpactParameterZ= new TH3F("hImpactParameterZ", ";Vertex contributors;#it{p}_{T} (GeV/#it{c});Impact parameter Z (cm)", 100, 0.5f, 100.5f, 200, 0.05f, 10.f, 1000, -1.f, 1.f);
      BinLogAxisY(hImpactParameterZ);

      std::vector<int> ITSidx;
      for (auto& sTrack : *strangetracks) {
        ITSidx.push_back(sTrack.mITSRef);
      }

      int countVertices = 0;

      std::ofstream prova;
      prova.open("prova3.txt");


      // for (unsigned int index = 0; index < 30; ++index)  ////loop on primary vertices 
      for (unsigned int index = 0; index < vertices->size(); ++index)  ////loop on primary vertices 
      {
        prova<<"\nVertex "<<countVertices+1<<"/"<< vertices->size()<<"\n"<<std::endl;
        std::cout<<"\nVertex "<<countVertices+1<<"/"<< vertices->size()<<"\n"<<std::endl;
        // get the track for each vertex and fill the tracks table
        // now go over tracks via the indices
        auto &v = vertices->at(index);
        prova<<"    v  "<<v<<std::endl;
        auto &trackref = v2trackref->at(index);
        prova<<"    trackref ---  "<<trackref<<std::endl;
        int start = trackref.getFirstEntry();
        prova<<"    start  "<<start<<std::endl;
        int ntracks = trackref.getEntries();    ///////number of thracks for the selected vertex
        prova<<"11111111111111111  ntracks "<<ntracks<<std::endl;

        for (int ti = 0; ti < ntracks; ++ti) //// loop on tracks for each vertex
        {
          auto trackindex = (*vertexTrackIDs)[start + ti];
          prova<<"Vertex "<<countVertices+1<<"/"<< vertices->size()<<"    Track "<<ti+1<<"/"<<ntracks<<"    Trackindex"<<trackindex<<"    "<<trackindex.getIndex()<<std::endl;


          // now we need to fetch the actual track and fill the table
          const auto source = trackindex.getSource();
          

          o2::track::TrackParCov *track = nullptr;
          int nCls = 0;
          if (source == o2::dataformats::VtxTrackIndex::Source::ITS) ////taking the its tracks so we already exclude the other tracks (to be controlled if it is needed)
          {

            
            unsigned int intIndex = trackindex.getIndex();
            auto ref = std::find(ITSidx.begin(), ITSidx.end(), intIndex); //// find the iterator index in the strange vector when the strange track index is equal to the its track index
            if (ref == ITSidx.end()) continue;  ////if the track is not found continue
            long stIndex=ref - ITSidx.begin();  ///// the position of the index has to be passed to strangetracks, not the index
            fmt::print("Track {}, ITS index {}, ST index {}/{}\n", ti, intIndex, stIndex,strangetracks->size());
            // auto& sTrack = (*strangetracks)[stIndex].mMother;
            auto& sTrack = (*strangetracks)[stIndex].mMother;
            track = &sTrack;  
          }
          else
          {
            continue;
          }
          gpu::gpustd::array<float, 2> dca{-999.f, -999.f};
          if (!propagator->propagateToDCA(v.getXYZ(), *track, -5.f, 2.f, o2::base::PropagatorF::MatCorrType::USEMatCorrLUT, &dca)) {
            continue;
          }
          hImpactParameterXY->Fill(v.getNContributors(), track->getPt(), dca[0]);
          hImpactParameterZ->Fill(v.getNContributors(), track->getPt(), dca[1]);
        }
        countVertices++;

      }
      
      hImpactParameterXY->Write();
      hImpactParameterZ->Write();
      prova<<"File Written"<<std::endl;
      f.Close();
    }
  }
}
