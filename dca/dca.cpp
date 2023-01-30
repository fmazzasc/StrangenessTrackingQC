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
#include "TH2F.h"
#include "TH1F.h"

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
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "Steer/MCKinematicsReader.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "SimulationDataFormat/DigitizationContext.h"
#include "StrangenessTracking/StrangenessTracker.h"
#include <vector>
#include <fstream>
#include "TSystemDirectory.h"
#include "TMath.h"

using namespace o2;
using namespace o2::framework;

using GIndex = o2::dataformats::VtxTrackIndex;

template <typename T>
void BinLogAxisY(T h)
{
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = h->GetYaxis();
  int bins = axis->GetNbins();
  double from = axis->GetXmin();
  double to = axis->GetXmax();
  double *newBins = new double[bins + 1];

  newBins[0] = from;
  double factor = std::pow(to / from, 1. / bins);

  for (int i = 1; i <= bins; i++)
  {
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
void dca(std::string path="/data/fmazzasc/its_data/sim/omega/", std::string outfile_str = "dca_omega.root")
{
  TSystemDirectory dir("MyDir", path.data());
  auto files = dir.GetListOfFiles();
  std::vector<std::string> dirs;
  std::vector<TString> kine_files;

  for (auto fileObj : *files)
  {
    std::string file = ((TSystemFile *)fileObj)->GetName();
    if (file.substr(0, 2) == "tf")
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

  TH2 *h2ImpactParameterXY, *h2ImpactParameterZ;
  h2ImpactParameterXY = new TH2F("h2ImpactParameterXY", ";#it{p}_{T} (GeV/#it{c});Impact parameter XY (cm)", 200, 0.05f, 10.f, 100, -.02f, .02f);
  h2ImpactParameterZ = new TH2F("h2ImpactParameterZ", ";#it{p}_{T} (GeV/#it{c});Impact parameter Z (cm)", 200, 0.05f, 10.f, 100, -.02f, .02f);

  TH1 *h1ImpactParameterXY, *h1ImpactParameterZ;
  h1ImpactParameterXY = new TH1F("h1ImpactParameterXY", ";Impact parameter XY (cm)", 50, -.02f, .02f);
  h1ImpactParameterZ = new TH1F("h1ImpactParameterZ", ";Impact parameter Z (cm)", 50, -.02f, .02f);

  int counter = 0;
  for (unsigned int i = 0; i < dirs.size(); i++)
  {
    LOG(info) << "Processing " << dirs[i];

    // open the file for vertices
    TFile f((dirs[i] + "/" + "o2_primary_vertex.root").data(), "OPEN");
    auto t = (TTree *)f.Get("o2sim");

    // fetch the tracks (these names are not following any convention!!)
    auto itstracks = fetchTracks<o2::its::TrackITS>((dirs[i] + "/" + "o2trac_its.root").data(), "o2sim", "ITSTrack");
    auto strangetracks = fetchTracks<o2::strangeness_tracking::StrangeTrack>((dirs[i] + "/" + "o2_strange_tracks.root").data(), "o2sim", "StrangeTracks");
    auto cascadetracks = fetchTracks<o2::dataformats::Cascade>((dirs[i] + "/" + "o2_secondary_vertex.root").data(), "o2sim", "Cascades");

    const auto grp = o2::parameters::GRPObject::loadFrom(dirs[i] + "/" + "o2sim_grp.root");
    o2::base::GeometryManager::loadGeometry(dirs[i] + "/" + "o2sim_geometry-aligned.root");
    o2::base::Propagator::initFieldFromGRP(grp);
    auto propagator = o2::base::Propagator::Instance();

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

        // BinLogAxisY(hImpactParameterZ);

        std::vector<int> ITSidx;
        for (auto &sTrack : *strangetracks)
        {
          ITSidx.push_back(sTrack.mITSRef);
        }

        int countVertices = 0;

        // for (unsigned int index = 0; index < 30; ++index)  ////loop on primary vertices
        for (unsigned int index = 0; index < vertices->size(); ++index) ////loop on primary vertices
        {

          // get the track for each vertex and fill the tracks table
          // now go over tracks via the indices
          auto &v = vertices->at(index);

          auto &trackref = v2trackref->at(index);

          int start = trackref.getFirstEntry();

          int ntracks = trackref.getEntries(); ///////number of thracks for the selected vertex

          for (int ti = 0; ti < ntracks; ++ti) //// loop on tracks for each vertex
          {
            auto trackindex = (*vertexTrackIDs)[start + ti];

            // now we need to fetch the actual track and fill the table
            const auto source = trackindex.getSource();

            o2::track::TrackParCov *track = nullptr;
            int nCls = 0;
            if (source == o2::dataformats::VtxTrackIndex::Source::ITS) ////taking the its tracks so we already exclude the other tracks (to be controlled if it is needed)
            {

              unsigned int intIndex = trackindex.getIndex();
              auto ref = std::find(ITSidx.begin(), ITSidx.end(), intIndex); //// find the iterator index in the strange vector when the strange track index is equal to the its track index
              if (ref == ITSidx.end())
                continue;                          ////if the track is not found continue
              long stIndex = ref - ITSidx.begin(); ///// the position of the index has to be passed to strangetracks, not the index

              // auto& sTrack = (*strangetracks)[stIndex].mMother;
              if ((*strangetracks)[stIndex].mPartType != 1)
                continue; ////if the track is not a cascade track continue
              auto &sTrack = (*strangetracks)[stIndex].mMother;
              auto &itsTrack = (*itstracks)[intIndex];
              track = &itsTrack;
              auto &casc = (*cascadetracks)[(*strangetracks)[stIndex].mDecayRef];
              LOG(info) << "\n--------------------------------------------";
              LOG(info) << "Vertex " << countVertices + 1 << "/" << vertices->size();
              LOG(info) << "Casc Vertex " << int(casc.getVertexID());
              fmt::print("Track {}, ITS index {}, ST index {}/{}\n", ti, intIndex, stIndex, strangetracks->size());

              if (casc.getVertexID() != countVertices)
                continue;
            }
            else
            {
              continue;
            }
            gpu::gpustd::array<float, 2> dca{-999.f, -999.f};
            if (!propagator->propagateToDCA(v.getXYZ(), *track, propagator->getNominalBz(), 2.f, o2::base::PropagatorF::MatCorrType::USEMatCorrLUT, &dca))
            {
              continue;
            }
            h2ImpactParameterXY->Fill(track->getPt(), dca[0]);
            h2ImpactParameterZ->Fill(track->getPt(), dca[1]);
            h1ImpactParameterXY->Fill(dca[0]);
            h1ImpactParameterZ->Fill(dca[1]);
          }
          countVertices++;
        }
      }
    }
  }
  TFile outfile(outfile_str.data(), "RECREATE");
  h2ImpactParameterXY->Write();
  h2ImpactParameterZ->Write();
  h1ImpactParameterXY->Write();
  h1ImpactParameterZ->Write();
  outfile.Close();
}
