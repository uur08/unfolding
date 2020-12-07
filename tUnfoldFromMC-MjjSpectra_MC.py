from ROOT import TFile, TF1, TUnfold, TUnfoldDensity
from ROOT import gDirectory, TH1F, TH1D, TH2D, TGraph, TSpline3, TH2F
from ROOT import gROOT, TCanvas, gStyle, TPad, TLine, TLegend
from ROOT import TMath, TText, TLatex, gApplication, gPad
from array import *
import math
import sys
import ROOT
import subprocess

gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("5.1f")

label = sys.argv[1]
print label
if label == "MC":
   MCfile = TFile("output-MC-1.root")
elif label == "MCherwig":
    MCfile = TFile("TTbar-Events-Herwig.root")
elif label == "MCfxfx":
    MCfile = TFile("TTbar-Events-CUEP8M2T4-FxFx.root")
elif label == "MCcuetm2":
    MCfile = TFile("TTbar-Events-CUETP8M2-UpJER-CentralBTagging.root")
else:
    print "Allowed arguments are MCcuetm1, MCherwig, MCfxfx, MCcuetm2"
    quit()

#dataFile = TFile("cleanfakes_massSpectra2016.root")
global bineta
bineta = ["0.0", "0.5", "1.0", "1.5", "2.0", "2.5"]

for i in range(0, 5):

    MCfile.cd("Standard/Eta_" + bineta[i] + "-" + bineta[i + 1] + "/mc")
    dLevel = gDirectory.Get("hdjmass")
    RM = gDirectory.Get("matrix_gen_reco")
    Gen = gDirectory.Get("hdjmass_gen")
    Miss =  gDirectory.Get("miss")
    Fake =  gDirectory.Get("fake")
     
    #dLevel = dataFile.Get('massSpectrum_'+bineta[i]+'-'+bineta[i+1])  # data reco level
    #Fake.Scale(Fake_norm.Integral()/Fake.Integral())
    
    print "Unfolding dijet mass spectrum for rapidity bin: " + bineta[i] + "-" + bineta[
        i + 1
    ]

    genbinning = [
	[160,  200,  249,  306,  372,  449,  539, 641,  756,  887, 1029, 1187, 1361, 1556, 1769, 2008, 2273, 2572, 2915, 3306, 3754, 4244, 4805, 5374, 6094, 6908, 7861, 8929, 10050],
	[160,  200,  249,  306,  372,  449,  539, 641,  756,  887, 1029, 1187, 1361, 1556, 1769, 2008, 2273, 2572, 2915, 3306, 3754, 4244, 4805, 5374, 6094, 6908, 7861, 8929, 10050],
	[160,  200,  249,  306,  372,  449,  539, 641,  756,  887, 1029, 1187, 1361, 1556, 1769, 2008, 2273, 2572, 2915, 3306, 3754, 4244, 4805, 5374, 6094, 6908, 7861, 8929, 10050],
	[160,  200,  249,  306,  372,  449,  539, 641,  756,  887, 1029, 1187, 1361, 1556, 1769, 2008, 2273, 2572, 2915, 3306, 3754, 4244, 4805, 5374, 6094, 6908, 7861, 8929, 10050],
	[160,  200,  249,  306,  372,  449,  539, 641,  756,  887, 1029, 1187, 1361, 1556, 1769, 2008, 2273, 2572, 2915, 3306, 3754, 4244, 4805, 5374, 6094, 6908, 7861, 8929, 10050]
      ]

    nbinsGen = len(genbinning[i]) - 1
    print "gen bin: ", nbinsGen
    correlationMatrix = TH2D(
        "correlationMatrix" + str(i) + "bin", "correlationMatrix" +
        str(i) + "bin", nbinsGen, array("d", genbinning[i]),
        nbinsGen,
        array("d", genbinning[i]),
    )

    recobinning = [
	[160,  200,  249,  306,  372,  449,  539, 641,  756,  887, 1029, 1187, 1361, 1556, 1769, 2008, 2273, 2572, 2915, 3306, 3754, 4244, 4805, 5374, 6094, 6908, 7861, 8929, 10050],
	[160,  200,  249,  306,  372,  449,  539, 641,  756,  887, 1029, 1187, 1361, 1556, 1769, 2008, 2273, 2572, 2915, 3306, 3754, 4244, 4805, 5374, 6094, 6908, 7861, 8929, 10050],
	[160,  200,  249,  306,  372,  449,  539, 641,  756,  887, 1029, 1187, 1361, 1556, 1769, 2008, 2273, 2572, 2915, 3306, 3754, 4244, 4805, 5374, 6094, 6908, 7861, 8929, 10050],
	[160,  200,  249,  306,  372,  449,  539, 641,  756,  887, 1029, 1187, 1361, 1556, 1769, 2008, 2273, 2572, 2915, 3306, 3754, 4244, 4805, 5374, 6094, 6908, 7861, 8929, 10050],
	[160,  200,  249,  306,  372,  449,  539, 641,  756,  887, 1029, 1187, 1361, 1556, 1769, 2008, 2273, 2572, 2915, 3306, 3754, 4244, 4805, 5374, 6094, 6908, 7861, 8929, 10050]
      ]

    nbinsRec = len(recobinning[i]) - 1
    print "reco bin: ", nbinsRec
    dCovMatrix = TH2D(
        "dCovMatrix" + str(i) + "bin",
        "dCovMatrix" + str(i) + "bin",
        nbinsRec,
        array("d", recobinning[i]),
        nbinsRec,
        array("d", recobinning[i]),
    )
    
    dCovMatrix.Sumw2()
    ### Ufl and Ofl
    RMx = RM.ProjectionX("RMx", 0, -1);
    if RMx.GetBinContent(0) != 0 or RMx.GetBinContent(nbinsGen+1):
        print "rm has underflow",RMx.GetBinContent(0),"rm has overflow",RMx.GetBinContent(nbinsGen+1)
    
    ### Correcting reco dist... (no need, since TUnfold provides SubstractBackground)
    #fake_fraction = dLevel.Clone()
    #fake_fraction.Divide(Fake, dLevel, 1, 1, "b")
    #dLevel.Multiply(1-fake_fraction)       
    
    ### Determining miss fraction...will be used after unfolding
    miss_frac = Miss.Clone()
    Gen.Add(Miss, -1);
    miss_frac.Divide(Miss,Gen, 1, 1, "b");
    
    # SET ERRORS OF DETECTOR LEVEL SPECTRA WITH COVARIANCE MATRIX OF THE MEASUREMENT
    for t in range(1, nbinsRec + 1):

        if dLevel.GetBinError(t) == 0:
            print "bin number ", t, " bin value ", dLevel.GetBinLowEdge(
                t
            ), " bin error ", dLevel.GetBinError(
                t
            ), " bin content ", dLevel.GetBinContent(
                t
            )
            continue
        dCovMatrix.SetBinContent(t, t, math.pow(dLevel.GetBinError(t),2))
        #print  t, ".Diagonal value is ",dCovMatrix.GetBinContent(t,t)," corresponding data input bin value is ", dLevel.GetBinLowEdge(t), " content of data " ,dLevel.GetBinContent(t)

    unfold = TUnfoldDensity(
        RM, TUnfold.kHistMapOutputHoriz, TUnfold.kRegModeNone, TUnfold.kEConstraintNone
    )
    # SET INPUT HERE
    unfold.SetInput(dLevel, 0, 0, 0)
    unfold.SubtractBackground(Fake,"bgr",1.0, 0.05);
    noReg = True

    if noReg:
        unfold.DoUnfold(0.00)

        print "( " + str(unfold.GetChi2A()) + "+" + str(
            unfold.GetChi2L()
        ) + ") / " + str(unfold.GetNdf())
        folderName = "NoReg_" + str(unfold.GetTau()) + "_y" + label
        subprocess.call(["mkdir", folderName])
        if i == "0.0":
            subprocess.call(["rm", folderName + "/unfoldedMjjSpectra0.0.root"])
    
    unfoldedHisty1 = unfold.GetOutput("HistoOutput" + str(i) + "bin")
    unfoldedErrory1 = unfold.GetEmatrixInput("unfolding stat error matrix")
    
    for ii in range(1, nbinsGen + 1):
        ei = math.sqrt(unfoldedErrory1.GetBinContent(ii, ii))
        unfoldedHisty1.SetBinError(ii, math.sqrt(
            unfoldedErrory1.GetBinContent(ii, ii)))

        if ei == 0:
#            print "x error"
#            print "bin number ", ii, " error ", ei
            continue
        for j in range(1, nbinsGen + 1):
            ej = math.sqrt(unfoldedErrory1.GetBinContent(j, j))
            if ej == 0:
#                print "y error"
#                print "bin number ", j, " error ", ej
                continue
            correlationMatrix.SetBinContent(
                ii, j, unfoldedErrory1.GetBinContent(ii, j) / ei / ej
            )
    for k in range(1, unfoldedHisty1.GetNbinsX()): 
        content = unfoldedHisty1.GetBinContent(k);
        if content < 0 : 
            error = unfoldedHisty1.GetBinError(k);
            print k,".bin has negatif content "  
        factor = 1
        factor += miss_frac.GetBinContent(k)
        content *= factor
        unfoldedHisty1.SetBinContent(k, content);

    dOutCovMatrix = dCovMatrix.Clone()
    # Covariance after unfolding
    for t in range(1, nbinsRec + 1):

        if unfoldedHisty1.GetBinError(t) == 0:
            print "bin number ", t, " bin value ", unfoldedHisty1.GetBinLowEdge(
                t
            ), " bin error ", unfoldedHisty1.GetBinError(
                t
            ), " bin content ", unfoldedHisty1.GetBinContent(
                t
            )
            continue
        dOutCovMatrix.SetBinContent(t, t, math.pow(unfoldedHisty1.GetBinError(t),2))
        #print  t, ".Diagonal value is ",dCovMatrix.GetBinContent(t,t)," corresponding data input bin value is ", dLevel.GetBinLowEdge(t), " content of data " ,dLevel.GetBinContent(t)


    ### Covariance matrices of uncorrelated systematics. Sources: RM,background,inefficiency 
    RM_Uncorr = unfold.GetEmatrixSysUncorr("dOutCovMatrix","RM_unc_cov")
    Bkg_Uncorr = unfold.GetEmatrixSysBackgroundUncorr("bgr","Bkg_unc_cov")
    Miss_Uncorr = Bkg_Uncorr.Clone()
    Cov_Uncorr_sum = Bkg_Uncorr.Clone()
    for k in range(1,Miss.GetNbinsX()):
        miss_E = Miss.GetBinError(k)
        if miss_E == 0:
            continue
        Miss_Uncorr.SetBinContent(k,k,math.pow(miss_E,2))
    Cov_Uncorr_sum.Add(RM_Uncorr)
    Cov_Uncorr_sum.Add(Bkg_Uncorr)
    Cov_Uncorr_sum.Add(Miss_Uncorr)
    
    ### Taking unc vectors. Sources: RM,background,inefficiency
        
    ### Background
    Bkg_shift = unfold.GetDeltaSysBackgroundScale("bgr","bgr_up","bgr_up")
    #Bkg_shift_d = unfold.GetDeltaSysBackgroundScale("bgr","bgr_down","bgr_down")
    #Bkg_shift_d.Scale(-1)

    Bkg_up = Fake.Clone()
    Bkg_down = Fake.Clone()

    for k in range(1,Fake.GetNbinsX()):
        bkg_nominal = Fake.GetBinContent(k)
        bkg_shift = Bkg_shift.GetBinContent(k)

        Bkg_up.SetBinContent(k, bkg_nominal + bkg_shift)
        Bkg_down.SetBinContent(k, bkg_nominal - bkg_shift)
    
    ### Miss
    Miss_shift = ROOT.TH1D(miss_frac)
    Miss_shift.Multiply(unfoldedHisty1)
    Miss_shift.Scale(0.05)
    #Miss_shift_d = ROOT.TH1D(Miss_shift_u)
    #Miss_shift_d.Scale(-1)

    Miss_up = Miss.Clone()
    Miss_down = Miss.Clone()

    for k in range(1,Miss.GetNbinsX()):
        miss_nominal = Miss.GetBinContent(k)
        miss_shift = Miss_shift.GetBinContent(k) 

        Miss_up.SetBinContent(k, miss_nominal + miss_shift)
        Miss_down.SetBinContent(k, miss_nominal - miss_shift)

    lx = unfold.GetLxMinusBias("Lx" + str(i) + "bin")
    c2 = TCanvas()
    lx.Draw()
    c2.Print(folderName + "/Lx" + str(i) + "Bin.pdf")

    c3 = TCanvas()
    correlationMatrix.GetXaxis().SetMoreLogLabels()
    correlationMatrix.GetYaxis().SetMoreLogLabels()
    correlationMatrix.GetYaxis().SetNoExponent()
    correlationMatrix.GetXaxis().SetNoExponent()

    correlationMatrix.SetTitle("#tau : " + str(unfold.GetTau()))
    correlationMatrix.Draw("colz text45")
    c3.SetLogy()
    c3.SetLogx()

    c3.Print(folderName + "/" + str(i) + "binsCorr" +
             str(unfold.GetTau()) + ".pdf")

    f = TFile(
        folderName + "/unfoldedSpectra-Mjj-Tau_" +
        str(unfold.GetTau()) + ".root",
        "update",
    )

    unfoldedHisty1.Write()
    #dLevel.Write("hdjmass_bin" + str(i))
    #dLevelMC.Write("hdjmass_gen_bin" + str(i))
    Fake.Write("fake" + str(i))
    Miss.Write("miss" + str(i))
    Miss_shift.Write("miss_shift" + str(i))
    #Fake_data.Write("P8_fake" + str(i))
    #Fake_MC.Write("HW_fake" + str(i))
    #RM.Write("normalizedRM" + str(i))
    #RM_MC.Write("mcRM" + str(i))
    dCovMatrix.Write("Input_cov" + str(i))
    #correlationMatrix.Write()
    #RM_Uncorr.Write("RM_cov" + str(i))
    #Bkg_unc.Write("Bkg_cov" + str(i))
    Bkg_up.Write("BKG_up"+str(i))
    Bkg_down.Write("BKG_down"+str(i))
    Miss_up.Write("Miss_up"+str(i))
    Miss_down.Write("Miss_down"+str(i))
    #MC_up.Write("Total_shift_up"+str(i))
    #MC_down.Write("Total_shift_down"+str(i))
    lx.Write()
    f.Close()
    print "Finished with rapidty bin " + bineta[i] + "-" + bineta[i + 1]
    print("\n")
