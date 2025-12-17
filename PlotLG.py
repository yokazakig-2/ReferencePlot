from ROOT import TH2F, TCanvas, TFile
import numpy as np

hlists=[]
inf=TFile("output.root","READ")
denh=inf.Get("xy25m_13mj1w0")

hfinal=TH2F("hE","",100,-50,50,100,-50,50)

for ii in range(0,10):
    ph=inf.Get("xy25m_13mj1w0_E0"+str(ii))
    ph.Divide(ph,denh,1,1)
    hlists.append(ph)

for xbin in range(1,101):
    for ybin in range(1,101):
        maxr = 0.
        maxii = 0
        if denh.GetBinContent(xbin,ybin) == 0: continue
        for rqE in range(1,10):
            if maxr < hlists[rqE].GetBinContent(xbin,ybin):
                maxr = hlists[rqE].GetBinContent(xbin,ybin)
                maxii = rqE
        hfinal.SetBinContent(xbin,ybin,float(maxii)*0.12+2.34)
        if ybin == 75 and maxii > 6: print(str(maxii))

contourlevel=np.linspace(2.28,3.48,11)
print(str(contourlevel))

c=TCanvas()
c.SetRightMargin(0.15)
hfinal.SetStats(0)
hfinal.SetMinimum(2.28)
hfinal.SetMaximum(3.48)
hfinal.GetXaxis().SetTitle("x [mm]")
hfinal.GetYaxis().SetTitle("y [mm]")
hfinal.GetZaxis().SetTitle("[GeV]")
hfinal.SetContour(10,contourlevel)
hfinal.Draw("COLZ")
c.SaveAs("Edis.pdf")
