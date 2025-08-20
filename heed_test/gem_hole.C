#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>

#include <cstdlib>
#include <iostream>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewField.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);

  // Load the field map.
    ComponentComsol fm;
    fm.Initialise("mesh.mphtxt", "dielectrics.dat", "field.txt", "mm");
    fm.EnableMirrorPeriodicityX();
    fm.EnableMirrorPeriodicityY();
    fm.PrintRange();

  // Enable Saving
    TFile* file = new TFile("../analysis/data/hole.root", "RECREATE");
    TTree* electronTree = new TTree("electrons", "GEM shower electrons");
    TTree* eventTree = new TTree("events", "GEM shower events");

    double gain, exf, eyf, ezf, etf, ee, exs, eys, ezs, ets;

    electronTree->Branch("xf", &exf, "xf/D");
    electronTree->Branch("yf", &eyf, "yf/D");
    electronTree->Branch("zf", &ezf, "zf/D");
    electronTree->Branch("tf", &etf, "tf/D");

    electronTree->Branch("xs", &exs, "xs/D");
    electronTree->Branch("ys", &eys, "ys/D");
    electronTree->Branch("zs", &ezs, "zs/D");
    electronTree->Branch("ts", &ets, "ts/D");

    electronTree->Branch("energy", &ee, "energy/D");
    eventTree->Branch("gain", &gain, "gain/D");

  // Dimensions of the GEM [cm]
    constexpr double pitch = 0.014;

    ViewField fieldView(&fm);
    // constexpr bool plotField = true;
    // if (plotField) {
    //   // Set the normal vector of the viewing plane (xz plane).
    //   fieldView.SetPlane(-1, 0, 0, 0, 0, 0);
    //   // Set the plot limits in the current viewing plane.
    //   fieldView.SetArea(-0.5 * pitch, -0.04, 0.5 * pitch, 0.015);
    //   fieldView.SetVoltageRange(-400., 400.);
    //   TCanvas* cf = new TCanvas("cf", "", 600, 600);
    //   cf->SetLeftMargin(0.16);
    //   fieldView.SetCanvas(cf);
    //   fieldView.PlotContour();
    // }

  // Setup the gas.
    MediumMagboltz gas("ar", 70., "co2", 30.);
    gas.SetTemperature(293.15);
    gas.SetPressure(760.);
    gas.SetMaxElectronEnergy(1e4);
    gas.Initialise(true);

  // Set the Penning transfer efficiency.
    constexpr double rPenning = 0.57;
    constexpr double lambdaPenning = 0.;
    gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");

  // Associate the gas with the corresponding field map material.
    fm.SetGas(&gas);
    fm.PrintMaterials();

  // Create the sensor.
    Sensor sensor(&fm);
    sensor.SetArea(-5 * pitch, -5 * pitch, -0.2, 5 * pitch, 5 * pitch, 0.05);



  AvalancheMicroscopic aval(&sensor);

  AvalancheMC drift(&sensor);
  drift.SetDistanceSteps(2.e-4);

  ViewDrift driftView;
  constexpr bool plotDrift = true;
  if (plotDrift) {
    aval.EnablePlotting(&driftView);
    drift.EnablePlotting(&driftView);
  }

  constexpr unsigned int nEvents = 2000;
  for (unsigned int i = 0; i < nEvents; ++i) {
    std::cout << i << "/" << nEvents << "\n";
    // Randomize the initial position.
    const double x0 = (RndmUniform()-0.5)*0.005;
    const double y0 = (RndmUniform()-0.5)*0.005; //* pitch + RndmUniform() * pitch;
    const double z0 = 0.0045;
    const double t0 = 0.;
    const double e0 = 0.1;
    aval.AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    int ne = 0, ni = 0;

    gain = 0;
    aval.GetAvalancheSize(ne, ni);
    for (const auto& electron : aval.GetElectrons()) {
      // const auto& p0 = electron.path.front();
      // drift.DriftIon(p0.x, p0.y, p0.z, p0.t);

      const auto& final = electron.path.back();
      const auto& start = electron.path.front();
      if (final.z < -0.01){
        exf = final.x; eyf = final.y; ezf = final.z; etf = final.t; ee = final.energy;
        exs = start.x; eys = start.y; ezs = start.z; ets = start.t;
        electronTree->Fill();
        gain++;
      }
      
    }
    
    if(gain>0) {eventTree->Fill();}
  }
  file->Write();
  file->Close();

  if (false) {
    TCanvas* cd = new TCanvas();
    constexpr bool plotMesh = true;
    if (plotMesh) {
      ViewFEMesh* meshView = new ViewFEMesh(&fm);
      std::cout << "asdsadadasda\n";
      meshView->SetArea(-2 * pitch, -2 * pitch, -0.04, 2 * pitch, 2 * pitch,
                        0.04);
      meshView->SetCanvas(cd);
      // x-z projection.
      meshView->SetPlane(0, 0, -1, 0, 0, 0);
      meshView->SetFillMesh(true);
      // Set the color of the kapton and the metal.
      meshView->SetColor(1, kYellow + 3);
      meshView->SetColor(2, kGray);
      // meshView->EnableAxes();
      meshView->SetViewDrift(&driftView);
      meshView->Plot(true);
    } else {
      driftView.SetPlane(0, 0, -1, 0, 0, 0);
      driftView.SetArea(-0.014, -0.04, 2 * pitch, 0.04);
      driftView.SetCanvas(cd);
      constexpr bool twod = false;
      driftView.Plot(twod);
    }
  }
  // app.Run();
}
