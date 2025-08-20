void plot_paths() {
  // Open the ROOT file
  TFile* file = TFile::Open("./data/other.root");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open file!" << std::endl;
    return;
  }

  // Get the TTree
  TTree* tree = (TTree*)file->Get("electrons");
  if (!tree) {
    std::cerr << "TTree 'paths' not found!" << std::endl;
    return;
  }

  // Set up branch
  std::vector<TVector3>* path = nullptr;
  tree->SetBranchAddress("path", &path);

  // Create canvas
  TCanvas* canvas = new TCanvas("c", "Electron Paths in XZ Plane", 800, 600);
  canvas->DrawFrame(-0.1, -0.6, 0.2, 0.4, "Electron Paths;X (cm);Z (cm)");

  // Loop over all entries and draw XZ paths
  int nEntries = tree->GetEntries();
  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (!path || path->empty()) continue;

    TPolyLine* line = new TPolyLine(path->size());
    for (size_t j = 0; j < path->size(); ++j) {
      line->SetPoint(j, (*path)[j].X(), (*path)[j].Z());  // XZ only
    }
    line->SetLineColor(kBlue + i % 5);  // color variation
    line->Draw("same");
  }

  canvas->Update();
}
