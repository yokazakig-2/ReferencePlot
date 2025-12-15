from ROOT import TGraph, TCanvas, kRed, kBlue, TLegend, TText
import ctypes

#Input CSV file
DataDir="Data/"
InputFileName="ARTBL_2025_June_Beam_Rate_AR5p0GeV.csv"
inf=open(DataDir+InputFileName,"r")
lines=inf.readlines()

OUTDir="PDF/"
OUTT0Name=InputFileName
OUTT0Name=OUTT0Name.replace(".csv","_T0rate.pdf")
OUTT04Name=InputFileName
OUTT04Name=OUTT04Name.replace(".csv","_T04rate.pdf")

gro=TGraph(100)
gro.SetName("Open")
gro.SetTitle("T0 rate [Hz]")
grc=TGraph(100)
grc.SetName("Close")
gro4=TGraph(100)
gro4.SetName("Open")
gro4.SetTitle("T0&T4 rate [Hz]")
grc4=TGraph(100)
grc4.SetName("Close")

np=0

for aline in lines:
	if "#" in aline: continue
	gro.SetPoint(np,float(aline.split(",")[2]), float(aline.split(",")[3]))
	grc.SetPoint(np,float(aline.split(",")[2]), float(aline.split(",")[5]))
	gro4.SetPoint(np,float(aline.split(",")[2]), float(aline.split(",")[4]))
	grc4.SetPoint(np,float(aline.split(",")[2]), float(aline.split(",")[6]))
	np+=1

gro.Set(np)
grc.Set(np)
gro4.Set(np)
grc4.Set(np)

c=TCanvas()
c.SetLeftMargin(0.15)
c.SetRightMargin(0.05)
c.SetGridx()
c.SetGridy()
gro.SetMarkerStyle(8)
gro.SetMarkerColor(kBlue)
gro.SetLineColor(kBlue)
grc.SetMarkerStyle(8)
grc.SetMarkerColor(kRed)
grc.SetLineColor(kRed)
gro4.SetMarkerStyle(8)
gro4.SetMarkerColor(kBlue)
gro4.SetLineColor(kBlue)
grc4.SetMarkerStyle(8)
grc4.SetMarkerColor(kRed)
grc4.SetLineColor(kRed)

gro.GetXaxis().SetTitle("[GeV]")
gro.GetYaxis().SetTitle("T0 rate [Hz]")
gro.Draw("APL")
grc.Draw("PL sames")
tl = TLegend(0.75,0.7,0.95,0.9)
tl.SetFillColor(0)
tl.AddEntry(gro,"Gap open","pl")
tl.AddEntry(grc,"Gap close","pl")
tl.Draw()
labels=[]
for i in range(gro.GetN()):
	xp = ctypes.c_double()
	yp = ctypes.c_double()
	gro.GetPoint(i, xp, yp)
	xval = xp.value
	if xval < 2: xval = xval-0.2
	if xval > 2: xval = xval+0.2
	yval = yp.value+900
	label = TText(xval,yval,f"{yp.value:.1f}")
	label.SetTextSize(0.03)
	label.SetTextAlign(21)
	label.SetTextColor(kBlue)
	labels.append(label)
	print(str(i)+" "+str(xp.value)+" "+str(yp.value))
	grc.GetPoint(i, xp, yp)
	xval = xp.value
	if xval < 2: xval = xval+0.2
	if xval > 2: xval = xval-0.2
	yval = yp.value
	if xval == 2. : yval = yval+1000.
	elif xval >= 1. and xval <= 5.5: yval = yval-1000
	elif xval > 5.5:
		yval = yval - 500.
		xval = xval + 0.5
	else :
		xval = xval + 0.2
		yval = yval + 300.
	label2 = TText(xval,yval,f"{yp.value:.1f}")
	label2.SetTextSize(0.03)
	label2.SetTextAlign(21)
	label2.SetTextColor(kRed)
	labels.append(label2)

print(str(labels))
for ilabel in labels:
	ilabel.Draw()

c.Update()
c.SaveAs(OUTDir+OUTT0Name)
gro4.GetXaxis().SetTitle("[GeV]")
gro4.GetYaxis().SetTitle("T0&T4 rate [Hz]")
gro4.Draw("APL")
grc4.Draw("PL sames")
tl.Draw()
labels=[]
for i in range(gro.GetN()):
	xp = ctypes.c_double()
	yp = ctypes.c_double()
	gro4.GetPoint(i, xp, yp)
	xval = xp.value
	if xval < 2: xval = xval-0.2
	if xval > 2: xval = xval+0.2
	yval = yp.value+300
	if xval < 0. :
		xval = xval - 0.1
		yval = yval - 150.
	label = TText(xval,yval,f"{yp.value:.1f}")
	label.SetTextSize(0.03)
	label.SetTextAlign(21)
	label.SetTextColor(kBlue)
	labels.append(label)
	grc4.GetPoint(i, xp, yp)
	xval = xp.value
	if xval < 2.5: xval = xval+0.2
	if xval > 2.5: xval = xval-0.2
	yval = yp.value
	if xval == 2.5 : yval = yval+300.
	elif xval >= 1.: yval = yval-300
	elif xval > 0.45:
		xval = xval + 0.2
		yval = yval + 30.
	elif xval > 0.35:
		xval = xval+0.1
		yval = yval-50.
	else:
		xval = xval+0.1
		yval = yval-150.
	print(str(i)+" "+str(xp.value)+" "+str(yp.value)+" "+str(xval))
	label2 = TText(xval,yval,f"{yp.value:.1f}")
	label2.SetTextSize(0.03)
	label2.SetTextAlign(21)
	label2.SetTextColor(kRed)
	labels.append(label2)

print(str(labels))
for ilabel in labels:
	ilabel.Draw()

c.Update()
c.SaveAs(OUTDir+OUTT04Name)
