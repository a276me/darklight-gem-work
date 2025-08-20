import ROOT

# Open the ROOT file
file = ROOT.TFile.Open("./data/hole_2mm.root")
if not file or file.IsZombie():
    print("Cannot open file!")
    exit()

# Get the TTree
tree = file.Get("electrons")  # Change to your actual TTree name
if not tree:
    print("TTree not found!")
    exit()

# Create histogram
nbins = 50
xmin = -0.06
xmax = 0.06
hist = ROOT.TH1D("hist", "X Distribution;X (cm);Counts", nbins, xmin, xmax)
hist.Sumw2()  # Enable proper error bar handling


# Create a variable to hold the branch value
value = ROOT.std.vector('double')()  # works if branch is a vector
# OR, if scalar branch:
# value = array('d', [0.0])
# tree.SetBranchAddress("mybranch", value)

# Use Draw method for convenience
tree.Draw("yf >> hist", "", "goff")  # "goff" = no graphics

# Fit with Gaussian
gauss = ROOT.TF1("gauss", "gaus", xmin, xmax)
fit_result = hist.Fit(gauss, "RS")  # "R" = restrict to range

# Extract fit stats
chi2 = fit_result.Chi2()
ndf = fit_result.Ndf()
reduced_chi2 = chi2 / ndf if ndf > 0 else float('inf')
# print(f"Chi2 / NDF = {chi2:.2f} / {ndf} = {chi2 / ndf:.2f}")


# Residuals histogram
resid = ROOT.TH1D("resid", "Residuals;X;Data - Fit", nbins, xmin, xmax)
for i in range(1, nbins + 1):
    x = hist.GetBinCenter(i)
    y_data = hist.GetBinContent(i)
    y_fit = gauss.Eval(x)
    resid.SetBinContent(i, y_data - y_fit)
    resid.SetBinError(i, hist.GetBinError(i))

canvas = ROOT.TCanvas("canvas", "X Distribution", 800, 800)
top_pad = ROOT.TPad("top", "Top Pad", 0, 0.4, 1, 1.0)      # Main histogram taller
bottom_pad = ROOT.TPad("bottom", "Bottom Pad", 0, 0.05, 1, 0.4)  # Residual shorter

top_pad.Draw()
bottom_pad.Draw()

# Upper plot
top_pad.cd()
hist.Draw("E")
gauss.SetLineColor(ROOT.kRed)
gauss.Draw("SAME")
ROOT.gPad.SetBottomMargin(0.01)
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.DrawLatex(0.15, 0.85, f"#chi^{{2}}/ndf = {reduced_chi2:.2f}")

# Lower plot
bottom_pad.cd()
resid.Draw("E")
ROOT.gPad.SetTopMargin(0.01)

canvas.Update()


ROOT.gApplication.Run()
