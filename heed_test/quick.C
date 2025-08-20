#include <iostream>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <random>

#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/Random.hh"
#include "Garfield/DriftLineRKF.hh"


using namespace Garfield;


// Setup Global Variables
static TFile* file = new TFile("../analysis/data/angle_res/0-40.root", "RECREATE");
static TTree* clusterTree = new TTree("clusters", "Heed Ionization Clusters");
static TTree* electronTree = new TTree("electrons", "GF++ Avalanche Electrons");
static TTree* eventTree = new TTree("events", "Heed Events");
static TTree* secondaryElectronTree = new TTree("sec_electrons", "Secondary ELectrons");
static ComponentComsol fm;
static double pitch = 0.014; // Pitch of the GEM readout in cm

static std::vector<double> esx, esy, esz, est;
std::vector<double> efx, efy, efz, eft, avg_delta_energy;
std::vector<double> ecz, ece, ecne;
static std::vector<int> nelectrons, sec_gain;
static int eventID, status, n_clusters;
static double place_holder, Angle;
static double true_x, true_y; // Set to the position in the middle of the drift chamber
static double start_x, start_y; // Set to position when first entering drift chamber
static std::vector<std::vector<TrackHeed::Electron>> events;
static std::vector<TVector3>* paths = new std::vector<TVector3>();

static bool randomize = true;

static const double xy_spread = 0.0118;
static const double t_spread = 2.2;
static const double t_mean = 28.5;
static const double gain_mean = 14;
static const double gain_t = 14;


static const double sense_l = 1000;

static const double layer1_z = 0; 
static const double layer2_z = -0.2;
static const double layer3_z = -0.4;
static const double readout_z = -0.6;
static const double drift_electrode_z = 0.3;

static const std::vector<double> layer_gains = {0.0001, 0.012 ,sqrt(2)*0.012, sqrt(3)*0.012};

static int num_electrons;
static const std::string gasfile = "ar80co220.gas";

int samplePolya(double theta, double meanGain) {
    double shape = 1.0 + theta;
    double scale = meanGain / shape;

    static std::random_device rd;
    static std::mt19937 gen(rd());

    std::gamma_distribution<> gammaDist(shape, scale);
    double sample = gammaDist(gen);

    return static_cast<int>(std::round(sample));  // or use std::floor(sample)

}

int det_layer(TrackHeed::Electron electron){

  int layer;
  if (electron.z > layer1_z) {
    layer = 3; // electron between 1st gem and drift electrode
  } else if (electron.z > layer2_z) {
    layer = 2; // electron between 1st and 2nd gem
  } else if (electron.z > layer3_z) {
    layer = 1; // electron between 2nd and 3rd gem
  } else {
    layer = 0; // electrons between 3rd gem and readout
  }
  return layer;
}

int main(int argc, char* argv[]) {
  // TApplication app("app", &argc, argv);

  // Setup Output
    electronTree->Branch("xf", &efx);
    electronTree->Branch("yf", &efy);
    electronTree->Branch("zf", &efz);
    electronTree->Branch("tf", &eft);
    electronTree->Branch("path", &paths);

    clusterTree->Branch("n_clusters", &n_clusters, "n_clusters/I");
    clusterTree->Branch("ecz", &ecz);
    clusterTree->Branch("ece", &ece);
    clusterTree->Branch("ecne", &ecne);

    secondaryElectronTree->Branch("x", &esx);
    secondaryElectronTree->Branch("y", &esy);
    secondaryElectronTree->Branch("z", &esz);
    secondaryElectronTree->Branch("t", &est);
    secondaryElectronTree->Branch("gain", &sec_gain);

    eventTree->Branch("theta", &Angle, "theta/D");
    eventTree->Branch("true_x", &true_x, "true_x/D");
    eventTree->Branch("true_y", &true_y, "true_y/D");

    eventTree->Branch("start_x", &start_x, "start_x/D");
    eventTree->Branch("start_y", &start_y, "start_y/D");

  // Setup field
  
    // Load the field map.
    fm.Initialise("mesh3.mphtxt", "dielectrics.dat", "field3.txt", "mm");
    fm.EnableMirrorPeriodicityX();
    fm.EnableMirrorPeriodicityY();

    fm.PrintRange();
    const double pitch = 0.014;
    const double d = pitch;

    ViewField fieldView(&fm);
    constexpr bool plotField = true;

  // Setup sensor
  
    Sensor sensor(&fm);
    sensor.SetArea(-sense_l*pitch, -sense_l*pitch, -0.6, sense_l*pitch, sense_l*pitch, 0.4);
    sensor.SetTimeWindow(-0.01, 90, 1000);
  
  
  // Setup Gas
    MediumMagboltz gas("ar", 80., "co2", 20.);
    bool calc_gas = true;
    // Set temperature [K] and pressure [Torr].
    gas.SetTemperature(293.15);
    gas.SetPressure(760.);
    gas.SetMaxElectronEnergy(1e5);
    gas.LoadGasFile(gasfile);
    
    // if(calc_gas){
    //   // gas.SetFieldGrid(-5.e3, 5.e3, 20, true);
    //   // gas.GenerateGasTable(7);
    //   // gas.WriteGasFile("gasgas.gas");
    //   gas.Initialise();
    // }else{
    //   gas.LoadGasFile("gasgas.gas");

    // }
    gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
    
    fm.SetGas(&gas);
    // gas.WriteGasFile("gas_file.gas");



  // Setup Tracking and Energy
    ViewDrift driftView;
    TrackHeed track(&sensor);
    track.EnablePlotting(&driftView);

    // Set the particle type and momentum [eV/c].
    track.SetParticle("electron");
    track.SetMomentum(10.e6);
    track.EnableDeltaElectronTransport();
    track.EnableElectricField();
    // track.EnableCoulombScattering();
  

  // Setup Avalanche
    
    AvalancheMC drift(&sensor);
    // Set the step size [cm].
    drift.SetDistanceSteps(1.e-4);
    

  // Track and Avalanche for n Events

    const int nEvents = 4000;
    double theta = 53;
    double phi = 0;

    for(int i=0; i<nEvents; i++){
    
      double x0 = 0.0, y0 = 0.0, z0 = 0.299, t0 = 0.;
      if(true){
        x0 = 1*RndmUniform() - 0.5;
        y0 = 1*RndmUniform() - 0.5;
        theta = 40*RndmUniform();
      }
      double dx = sin(theta * M_PI / 180.0)*cos(phi * M_PI / 180.0), 
             dy = sin(theta * M_PI / 180.0)*sin(phi * M_PI / 180.0), 
             dz = -cos(theta * M_PI / 180.0);
      
      std::cout << "Event " << i << " start position: (" << x0 << ", " << y0 << ", " << z0 << ")" << std::endl;

      track.NewTrack(x0, y0, z0, t0, dx, dy, dz);
      
      
    
      bool bad = true;
      for(auto& c : track.GetClusters()){
        if (c.z < layer3_z){bad=false;}
      }

      if(track.GetClusters().size() < 1){
        std::cout << "Detected Bad Track, No Ionizations. Redoing Event" << std::endl;
        i--;
        continue;
      }
      if(bad){
        std::cout << "Detected Bad Track, Track Disrupted. Redoing Event" << std::endl;
        i--;
        continue;
      }

      eventID = i;
      std::vector<TrackHeed::Electron> electrons;

      ece.clear();
      ecz.clear();
      ecne.clear();

      n_clusters = track.GetClusters().size();
      
      for (const auto& cluster : track.GetClusters()) {
        double tot_delta_energy = 0;
        
        int ne = 0;
        for (const auto& electron : cluster.electrons) {
          // Get the coordinates of the electron.
          tot_delta_energy += electron.e;
          ne++;
          electrons.push_back(electron);
        }

        nelectrons.push_back(ne);
        if (ne > 0) { avg_delta_energy.push_back(tot_delta_energy / ne);}

        ecne.push_back(ne);
        ecz.push_back(cluster.z);
        ece.push_back(cluster.energy);
        
      }

      clusterTree->Fill();
      
      // bool bad = true;
      // for(const auto& e : electrons){
      //   if(e.z < layer3_z){
      //     bad = false;
      //   }
      // }
      // if(bad){
      //   std::cout << "Detected Bad Track, Redoing Event" << std::endl;
      //   i--;
      //   continue;
      // }


      std::cout<<"Event " << i << " Finished " << "\n";
      Angle = theta;

      start_x = x0;
      start_y = y0;

      true_x = x0 + 0.15* tan(theta * M_PI / 180.0) * cos(phi * M_PI / 180.0);
      true_y = y0 + 0.15* tan(theta * M_PI / 180.0) * sin(phi * M_PI / 180.0);

      eventTree->Fill();
      std::cout<<"pushback"<<std::endl;
      events.push_back(electrons);

    }
    std::cout<< "HEED Simulation Finished" << std::endl;

    
    int id = 0;
    for(auto& event : events){
      std::cout<<"Sampling Event " << id << std::endl;
      efx.clear(); efy.clear(); esx.clear();esz.clear();est.clear();esy.clear();
      sec_gain.clear();
      for(auto& e : event){
        esx.push_back(e.x);
        esy.push_back(e.y);
        esz.push_back(e.z);
        est.push_back(e.t);

        int layer = det_layer(e);
        double g = 1;
        for(int i=layer; i>0; i--){
          g *= samplePolya(gain_t, gain_mean)+1;
        }
        sec_gain.push_back(g);
        // std::mt19937 generator( (int) RndmUniform()*100);
        // std::normal_distribution<double> distx(e.x, layer_gains[layer]);
        // std::normal_distribution<double> disty(e.y, layer_gains[layer]);
        // // std::cout<<e.x << e.y<<std::endl;

        // for(int i=0; i<g; i++){
        //   // sefx.push_back(e.x);
        //   // sefy.push_back(disty(generator));

        //   efx.push_back(e.x);
        //   // efy.push_back(-0.1);
        // }
        // efx = std::move(sefx);
        // efy = std::move(sefy);
        
      }
      electronTree->Fill();
      secondaryElectronTree->Fill();
      
      id++;
    }

    std::cout<<"Random Sampling Finished, Writing to ROOT" << std::endl;
    file->Write();
    file->Close();
    std::cout<<"Written to Disk\n";



  return 0;
}












