#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>

#include <cstdlib>
#include <iostream>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ComponentUser.hh"

#include <TFile.h>
#include <TTree.h>

using namespace Garfield;

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);

  TH1F* hIonTime = new TH1F("hIonTime", "Ion Arrival Times;Time [ns];Counts", 100, 500, 600);
  
  TFile* file = new TFile("3730V.root", "RECREATE");
  TTree* T = new TTree("T", "Tree");
  double ionTime = 0;
  T->Branch("ionTime", &ionTime, "ionTime/D");


  TH1::SetDefaultSumw2();
  TH1::AddDirectory(kFALSE);
  hIonTime->SetCanExtend(TH1::kAllAxes);

  // Setup the gas.
  MediumMagboltz gas("ar", 70., "co2", 30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.Initialise(true);
  // Set the Penning transfer efficiency.
  constexpr double rPenning = 0.51;
  constexpr double lambdaPenning = 0.;
  gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");

  // Load the field map.
  // ComponentAnsys123 fm;
  // fm.Initialise("ELIST.lis", "NLIST.lis", "MPLIST.lis", "PRNSOL.lis", "mm");
  // fm.EnableMirrorPeriodicityX();
  // fm.EnableMirrorPeriodicityY();
  // fm.PrintRange();

  ComponentUser fm;
  fm.SetArea(-10, -10, -10, 10, 10 , 10 );
  fm.SetMedium(&gas);
  auto efield = [](const double x, const double y, const double z,
    double& ex, double& ey, double& ez) {
    ey = ex = 0.;
    ez = 3730.;
  };
  fm.SetElectricField(efield);

  // Associate the gas with the corresponding field map material.
  // fm.SetGas(&gas);
  // fm.PrintMaterials();
  // fm.Check();

  // Dimensions of the GEM [cm]
  constexpr double pitch = 0.014;


  // Create the sensor.
  Sensor sensor(&fm);
  sensor.SetArea(-5, -5, 0, 5, 5, 1.1);

  AvalancheMicroscopic aval(&sensor);

  AvalancheMC drift(&sensor);
  drift.SetDistanceSteps(1.e-5);

  ViewDrift driftView;
  constexpr bool plotDrift = true;

  // Count the total number of ions produced the back-flowing ions.
  unsigned int nTotal = 0;
  unsigned int nBF = 0;
  constexpr unsigned int nEvents = 500;
  for (unsigned int i = 0; i < nEvents; ++i) {
    std::cout << i << "/" << nEvents << "\n";
    // Randomize the initial position.
    const double x0 = 0.;
    const double y0 = 0.;
    const double z0 = 1.;
    const double t0 = 0.;
    const double e0 = 0.;
    aval.DriftElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    for (const auto& electron : aval.GetElectrons()) {
      const auto& p0 = electron.path[0];

      hIonTime->Fill(electron.path.back().t);
      ionTime = electron.path.back().t;
      T->Fill();
      std::cout << "Electron status: " << electron.status << "\n";
      std::cout << "Drifted ion z from: " << electron.path.front().z << " to "
                << electron.path.back().z << "\n";
      std::cout << "Drifted ion x from: " << electron.path.front().x << " to "
                << electron.path.back().x << "\n";
      std::cout << "Drifted ion y from: " << electron.path.front().y << " to "
                << electron.path.back().y << "\n";
      std::cout << "Time of  " 
                << electron.path.back().t << "\n";
    }
  }
  TCanvas* c1 = new TCanvas("c1", "Ion Arrival Times", 800, 600);
  hIonTime->Draw();
  c1->Update();
  file->Write();
  file->Close();
  app.Run();
}
