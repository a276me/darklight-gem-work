void test() {
  // Create a canvas
  TCanvas* c1 = new TCanvas("c1", "Example Histogram", 800, 600);

  // Create a histogram: 100 bins from 0 to 10
  TH1F* h = new TH1F("h", "Random Gaussian;X;Entries", 100, 0, 10);

  // Fill the histogram with 10,000 Gaussian random numbers (mean=5, sigma=1)
  for (int i = 0; i < 10000; ++i) {
    h->Fill(gRandom->Gaus(5, 1));
  }

  // Draw the histogram
  h->Draw();

  // Save the canvas as a PNG
  c1->SaveAs("histogram.png");
}