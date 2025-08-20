void xy_dist() {
  // Load the ROOT file
  TFile *file = TFile::Open("./data/hole_characterization.root");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open file!" << std::endl;
    return;
  }

  // Get the TTree from the file
  TTree *tree = (TTree*)file->Get("x"); // Change "mytree" to your actual TTree name
  if (!tree) {
    std::cerr << "TTree not found!" << std::endl;
    return;
  }

  // Variable to hold branch data
  double value;
  tree->SetBranchAddress("mybranch", &value); // Change "mybranch" to your actual branch name

  // Create a histogram
  int nbins = 100;
  double xmin = 0, xmax = 100; // Adjust these to your expected data range
  TH1D *hist = new TH1D("hist", "Histogram with Gaussian Fit;X;Counts", nbins, xmin, xmax);
  hist->Sumw2(); // Enable storage of sum of squares for error bars

  // Fill the histogram
  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    hist->Fill(value);
  }

  // Create canvas to draw
  TCanvas *canvas = new TCanvas("canvas", "Histogram with Fit", 800, 600);
  hist->Draw("E"); // "E" option draws error bars

  // Fit with Gaussian
  TF1 *gaussFit = new TF1("gaussFit", "gaus", xmin, xmax);
  hist->Fit(gaussFit, "R"); // "R" restricts fit to function range

  // Draw the fit
  gaussFit->SetLineColor(kRed);
  gaussFit->Draw("SAME");

  // Optional: Print fit parameters
  std::cout << "Mean: " << gaussFit->GetParameter(1)
            << ", Sigma: " << gaussFit->GetParameter(2) << std::endl;

  // Keep the canvas open
  canvas->Update();
}