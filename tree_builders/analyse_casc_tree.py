import uproot
import pandas as pd
import ROOT
from utils import *
# ROOT.gROOT.LoadMacro('RooCustomPdfs/RooDSCBShape.cxx++')
# from ROOT import RooDSCBShape

kBlueC = ROOT.TColor.GetColor('#1f78b4')
kOrangeC  = ROOT.TColor.GetColor("#ff7f00")


############################################
mass_casc = ROOT.TH1F("mass_casc", "mass_casc", 100, 1.26, 1.44)
mass_casc_bkg = ROOT.TH1F("mass_casc_bkg", "mass_casc_bkg", 100, 1.26, 1.44)

mass_tracked = ROOT.TH1F("mass_tracked", "mass_tracked", 100, 1.26, 1.44)
##only for true hypertritons
dca_z_casc = ROOT.TH1F("dca_z_casc", "dca_z_casc", 100, -.5, .5)
dca_z_tracked = ROOT.TH1F("dca_z_tracked", "dca_z_tracked", 100, -.05, .05)
dca_xy_casc = ROOT.TH1F("dca_xy_casc", "dca_xy_casc", 100, -.5, .5)
dca_xy_tracked = ROOT.TH1F("dca_xy_tracked", "dca_xy_tracked", 100, -.05, .05)

mom_gen = ROOT.TH1F("mom_gen", "mom_gen", 40, 0, 10)
mom_casc = ROOT.TH1F("mom_casc", "mom_casc", 40, 0, 10)
mom_trackable = ROOT.TH1F("mom_trackable", "mom_trackable", 40, 0, 10)
mom_tracked = ROOT.TH1F("mom_tracked", "mom_tracked", 40, 0, 10)
dec_rad_casc = ROOT.TH1F("dec_rad_casc", "dec_rad_casc", 40, 0, 40)
dec_rad_trackable = ROOT.TH1F("dec_rad_trackable", "dec_rad_trackable", 40, 0, 40)
dec_rad_tracked = ROOT.TH1F("dec_rad_tracked", "dec_rad_tracked", 40, 0, 40)

mom_reso = ROOT.TH2F("mom_reso", "mom_reso", 40, 0, 10, 40, -0.1, 0.1)

############################################


df_casc = uproot.open("TrackedCascTree_xi_test.root")["CascTree"].arrays(library="pd")
print(df_casc)

df_casc_true = df_casc.query("gPt > 0") ## remove background
df_casc_bkg = df_casc.query("gPt < 0") ## get background

fill_th1_hist(df_casc_true, "gPt", mom_gen)
df_casc_true.query("isReconstructed == 1", inplace=True) ## select reconstructed cascades
df_casc_true.eval("decR = sqrt(gR2)", inplace=True)
df_casc_true.eval("reso = (gPt - recoPt)/gPt", inplace=True)


df_trackable = df_casc_true.query("hasITSTrack == 1")
df_tracked = df_casc_true.query("isStTracked == 1 and isFakeMatching == 0")


fill_th1_hist(df_casc_true, "recoMass", mass_casc)
fill_th1_hist(df_casc_bkg, "recoMass", mass_casc_bkg)


fill_th1_hist(df_tracked, "recoMass", mass_tracked)


fill_th1_hist(df_casc_true, "recoDCAZ", dca_z_casc)
fill_th1_hist(df_casc_true, "recoDCAXY", dca_xy_casc)
fill_th2_hist(df_casc_true, "gPt", "reso", mom_reso)
fill_th1_hist(df_tracked, "trackedDCAZ", dca_z_tracked)
fill_th1_hist(df_tracked, "trackedDCAXY", dca_xy_tracked)

fill_th1_hist(df_casc_true, "gPt", mom_casc)
fill_th1_hist(df_trackable, "gPt", mom_trackable)
fill_th1_hist(df_tracked, "gPt", mom_tracked)

fill_th1_hist(df_casc_true, "decR", dec_rad_casc)
fill_th1_hist(df_trackable, "decR", dec_rad_trackable)
fill_th1_hist(df_tracked, "decR", dec_rad_tracked)

########## fitting invariant masses


### sum w2 for all the histograms
mass_casc.Sumw2()
mass_casc_bkg.Sumw2()
mass_tracked.Sumw2()
dca_z_casc.Sumw2()
dca_z_tracked.Sumw2()
dca_xy_casc.Sumw2()
dca_xy_tracked.Sumw2()
mom_casc.Sumw2()
mom_trackable.Sumw2()
mom_tracked.Sumw2()
dec_rad_casc.Sumw2()
dec_rad_trackable.Sumw2()
dec_rad_tracked.Sumw2()



mass_tracked.SetStats(0)
mass_tracked.SetStats(0)
mass_tracked.GetYaxis().SetTitle("Counts")
mass_tracked.GetXaxis().SetTitle("#it{M}(^{3}He + #pi^{-} and c.c.)   (GeV/#it{c}^{2})")
mass_casc.GetYaxis().SetTitle("Counts")
mass_casc.GetXaxis().SetTitle("#it{M}(^{3}He + #pi^{-} and c.c.)   (GeV/#it{c}^{2})")

# mass = ROOT.RooRealVar('m', '#it{M}(^{3}He + #pi^{-} and c.c.)', 1.26, 1.44, 'GeV/c^{2}')
# mu = ROOT.RooRealVar('mu', 'hypernucl mass', 1.26, 1.44, 'GeV/c^{2}')
# sigma = ROOT.RooRealVar('sigma', 'hypernucl width', 0.001, 0.004, 'GeV/c^{2}')
# a1 = ROOT.RooRealVar('a1', 'a1', 0, 3.)
# a2 = ROOT.RooRealVar('a2', 'a2', 0.3, 0.8)

# c0 = ROOT.RooRealVar('c0', 'constant c0', -100., 100)
# c1 = ROOT.RooRealVar('c1', 'constant c1', -100., 100)
# c2 = ROOT.RooRealVar('c2', 'constant c2', -100., 100)

# n1 = ROOT.RooRealVar('n1', 'n1', 1, 10.)
# n2 = ROOT.RooRealVar('n2', 'n2', 1, 10.)
# signal = ROOT.RooDSCBShape('cb', 'cb', mass, mu, sigma, a1, n1, a2, n2)

# background1 = ROOT.RooChebychev('bkg', 'pol3 bkg', mass, ROOT.RooArgList(c0,c1, c2))
# background2 = ROOT.RooChebychev('bkg', 'pol1 bkg', mass, ROOT.RooArgList(c0))

# n = ROOT.RooRealVar('n', 'n const', 0.01, 1)
# # define the fit funciton and perform the actual fit
# fit_function1 = ROOT.RooAddPdf('total_pdf1', 'signal + background 1', ROOT.RooArgList(signal, background1), ROOT.RooArgList(n))
# fit_function2 = ROOT.RooAddPdf('total_pdf2', 'signal + background 2', ROOT.RooArgList(signal, background2), ROOT.RooArgList(n))


outfile = ROOT.TFile("histos.root", "recreate")

# _, _, frame1 = fit_and_plot(mass_casc, mass, fit_function1, signal, background1, sigma, mu, n)
# _, _, frame2 = fit_and_plot(mass_tracked, mass, fit_function2, signal, background2, sigma, mu, n)


## rescale background and signal
# currently, one signal event per bkg is generated. to be rescaled: 0.8*1e-7




cv1 = ROOT.TCanvas()
# frame1.Draw()
leg1 = ROOT.TLegend(0.52,0.52,0.93,0.76)
leg1.AddEntry("data","{}_{#Lambda}^{3}H + {}_{#bar{#Lambda}}^{3}#bar{H}", "PE")
leg1.AddEntry("fit_func","Signal + Background", "L")
leg1.AddEntry("bkg","Background", "L")
leg1.SetBorderSize(0)
leg1.SetFillStyle(0)
leg1.SetTextFont(42)
leg1.Draw()
cv1.Write()

cv2 = ROOT.TCanvas()
# frame2.Draw()
cv2.Write()

mass_casc_bkg.Write()
mass_tracked.Write()
mass_casc.Write()
dca_z_tracked.Write()
dca_z_casc.Write()
dca_xy_tracked.Write()
dca_xy_casc.Write()

pinfo = ROOT.TPaveText(0.7,0.7,0.9,0.9,"NDC")
# pinfo.AddText("This thesis, MC simulation")
pinfo.AddText("pp, #sqrt{#it{s}_{NN}} = 13.6 TeV")
pinfo.SetBorderSize(0)
pinfo.SetFillStyle(0)
pinfo.SetTextAlign(11)
pinfo.SetTextFont(42)



## efficiency plots
eff_rec_gen = mom_casc.Clone("eff_rec_gen")
eff_rec_gen.Divide(mom_gen)
eff_rec_gen.SetStats(0)
eff_rec_gen.GetYaxis().SetTitle("Reconstruction efficiency")
eff_rec_gen.GetXaxis().SetTitle("#it{p}_{T}   (GeV/#it{c})")
eff_rec_gen.SetTitle("")
eff_rec_gen.Write()






eff_algo = mom_tracked.Clone("eff_algo_pt")
eff_algo.Divide(mom_trackable)
eff_algo.SetStats(0)
eff_algo.GetYaxis().SetTitle("Algorithm efficiency")
eff_algo.GetXaxis().SetTitle("#it{p}_{T}   (GeV/#it{c})")
eff_algo.SetTitle("")

cv3 = ROOT.TCanvas()
eff_algo.Draw("pe")
pinfo.Draw()
cv3.Write()

eff_algo_rad = dec_rad_tracked.Clone("eff_algo_rad")
eff_algo_rad.Divide(dec_rad_trackable)
eff_algo_rad.SetStats(0)
eff_algo_rad.GetYaxis().SetTitle("Algorithm efficiency")
eff_algo_rad.GetXaxis().SetTitle("Decay radius (cm)")
eff_algo_rad.SetTitle("")
eff_algo_rad.Write()

cv4 = ROOT.TCanvas()
eff_algo_rad.Draw("pe")
pinfo.Draw()
cv4.Write()

eff_found_rad = dec_rad_tracked.Clone("eff_found_rad")
eff_found_rad.Divide(dec_rad_casc)
eff_found_rad.SetStats(0)
eff_found_rad.GetYaxis().SetTitle("Tracked / Found")
eff_found_rad.GetXaxis().SetTitle("Decay radius (cm)")
eff_found_rad.SetTitle("")
eff_found_rad.Write()

cv4 = ROOT.TCanvas()
eff_found_rad.Draw("pe")
pinfo.Draw()
cv4.Write()




eff_found = mom_tracked.Clone("eff_found_pt")
eff_found.Divide(mom_casc)
eff_found.SetStats(0)
eff_found.GetYaxis().SetTitle("Tracked / Found")
eff_found.GetXaxis().SetTitle("#it{p}_{T}   (GeV/#it{c})")
eff_found.SetTitle("")

cv5 = ROOT.TCanvas()
eff_found.Draw()
pinfo.Draw()
cv5.Write()

mom_reso.Write()


print(df_casc[['gPt', 'recoPt']])