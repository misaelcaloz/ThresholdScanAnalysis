						#
						#
						#
		#################################################################################
		# Call script: python -i threshold_CCPDv2.py DATAFILENAME ROOTFILENAME.root		#
		#																				#
		# 	Script to plot the threshold distribution from								#
		#	threshold scan using python functions 										#
		#################################################################################
						#
						#
						#
						#
# --------packages 

import sys
sys.path.insert(0, "/Users/misaelcaloz/Documents/GitHub/FunctionsRepository")
#sys.path.insert(0, "/Users/misaelcaloz/Documents/CERN master thesis/python scripts/functions")
sys.path.append("/usr/local/lib/root")
#import ROOT
from ROOT import gROOT, TCanvas, TF1, TGraph, TLegend, TMath, TMultiGraph, TFile, TH2F, TH1F, TDirectory,TH2D, TH1D, TLatex, gPad, TH2I, TLine, TAxis, kTRUE, gStyle, TKey, THStack, kRed, kOrange, kBlue, kGreen
#from ROOT import *
from array import array
import numpy
import math
import pylab
from functions import *

# -------- Define sys.argv

FileDataThrScan1 = sys.argv[1]
FileDataThrScan2 = sys.argv[2]
rootFileName = sys.argv[3]



# -------- Create .root file

outFile=TFile(rootFileName,"RECREATE")

# -------- Get histogram 1

FileThrScan1 = ROOT.TFile(FileDataThrScan1,"read")

print FileThrScan1.ls()

CanvThrScan1 = TCanvas("CanvThrScan1","CanvThrScan1")

thresh2Dmap1 = TH2F("thresh2Dmap1", "Thresh2Dmap;#Row;#Column;Threshold [e]",24,0,24, 36,12,48)
thresh2Dmap1 = FileThrScan1.Get("thresh2D")
thresh2Dmap1.Draw("colz")
CanvThrScan1.Update()


# -------- Get histogram 2

FileThrScan2 = ROOT.TFile(FileDataThrScan2,"read")

print FileThrScan2.ls()

CanvThrScan2 = TCanvas("CanvThrScan2","CanvThrScan2")

thresh2Dmap2 = TH2F("thresh2Dmap2", "Thresh2Dmap;#Row;#Column;Threshold [e]",24,0,24, 36,12,48)
thresh2Dmap2 = FileThrScan2.Get("thresh2D")
thresh2Dmap2.Draw("colz")

CanvThrScan2.Update()


# ------ Create correlation

Correlation_plot = TH2F("Correlation_plot", "Correlation;Temperature = 0 degree celsius ;Temperature = 20 degree celsius;Threshold [e]",50,0,2000,50,0,2000)
Correlation_plot_2Dmap = TH2F("Correlation_plot_2Dmap", "Correlation 2D map;#Row;#Col;Difference [e]",24,0,24, 36,12,48)

for r in range(1,25):
	for c in range(1,37):
		row = r-1
		col = c+11
		BinHist1 = thresh2Dmap1.GetBinContent(r,c)
		BinHist2 = thresh2Dmap2.GetBinContent(r,c)
		
		Diff_BinHist = BinHist1-BinHist2
		
		Correlation_plot.Fill(BinHist1,BinHist2)
		Correlation_plot_2Dmap.Fill(row,col,Diff_BinHist)
		
		print BinHist1
		print BinHist2



# ------- Customize plots

# -- Correlation plot

#Diag_Line =  TF1("Diag_Line","x",0,1700)

x = [0.00000001,1750]

x_array=array("d",x)

gr = TGraph(2,x_array,x_array)

Correlation_plot.SetAxisRange(-1.,2000.,"X")
Correlation_plot.SetAxisRange(-1.,2000.,"Y")
Correlation_plot.SetTitleSize(0.025,"xyz")

Correlation_plot.SetTitleOffset(1.3,"z")
Correlation_plot.SetTitleOffset(1.4,"y")


Correlation_plot.SetLabelSize(0.02,"X")
Correlation_plot.SetLabelSize(0.02,"Y")
Correlation_plot.GetZaxis().SetLabelSize(0.025)


# -- 2D map

#Correlation_plot_2Dmap.SetAxisRange(5.,150.,"Z")

Correlation_plot_2Dmap.SetTitleSize(0.025,"xyz")
Correlation_plot_2Dmap.SetTitleOffset(1.7,"z")

Correlation_plot_2Dmap.GetXaxis().SetNdivisions(32)
Correlation_plot_2Dmap.GetYaxis().SetNdivisions(64)
Correlation_plot_2Dmap.SetLabelSize(0.02,"X")
Correlation_plot_2Dmap.SetLabelSize(0.02,"Y")
Correlation_plot_2Dmap.GetZaxis().SetLabelSize(0.025)

x_axis_Correlation_plot_2Dmap = Correlation_plot_2Dmap.GetXaxis()
x_axis_Correlation_plot_2Dmap.CenterLabels(kTRUE)
y_axis_Correlation_plot_2Dmap = Correlation_plot_2Dmap.GetYaxis()
y_axis_Correlation_plot_2Dmap.CenterLabels(kTRUE)

Correlation_plot_2Dmap.SetMaximum(1000)
Correlation_plot_2Dmap.SetMinimum(-1000)



# -------- Create canvas




CanvCorrelation = ROOT.TCanvas("Correlation","Correlation",1200,600)
ROOT.SetOwnership(CanvCorrelation,False) # TODO: DOESN'T WORK...	
CanvCorrelation.Divide(2,1)
CanvCorrelation.cd(1)	
gStyle.SetOptStat("men")


NRGBs=5
NCont =255
stops = [0.00,0.34,0.61,1.00,1.00]
red = [0.00,0.00,0.87,1.00,1.00]
green = [0.00,1.00,0.80,0.20,0.00]
blue = [1.00,1.00,0.20,0.00,0.00]

stopsArray = array('d',stops)
redArray = array('d',red)
greenArray = array('d',green)
blueArray = array('d',blue)

ROOT.TColor.CreateGradientColorTable(NRGBs,stopsArray,redArray,greenArray,blueArray,NCont)
gStyle.SetNumberContours(NCont)
gPad.SetLogz(1)
#gStyle.SetOptStat(0)
gPad.SetRightMargin(0.15)




Correlation_plot.Draw("colz")
gr.Draw("same")
gStyle.SetStatX(0.85)
gStyle.SetStatY(0.95)
gPad.Update()





CanvCorrelation.cd(2)
gStyle.SetOptStat("e")
#gPad.SetLogz(0)
Correlation_plot_2Dmap.Draw("colz")
gStyle.SetStatX(0.9)
gStyle.SetStatY(0.95)
gPad.SetRightMargin(0.15)

gPad.Update()

list_TLine3=[]
for i in range(0,25):
	line_vert =  TLine(i,12,i,48)
	line_vert.Draw()
	list_TLine3.append(line_vert)

for i in range(12,49):
	line_hor =  TLine(0,i,24,i)
	line_hor.Draw()
	list_TLine3.append(line_hor)
	
gPad.SetLogz(0)

gPad.Update()




#k = TKey()
#k = FileThrScan1.GetKey("thresh2D")
#histo = k.ReadObj().FindObject("Threshold 2D plot")
##
#c1 = TCanvas("c1_test","c1_test")
#c1.cd()
#histo.Draw()
#
#c2 = TCanvas("c1","c1")
#
#for col in range(142,166):
#	for raw in range(1,13):
#		noise=histo.GetBinContent(raw,col)
#		noiseocc_hist_091.Fill(noise)
#		if noise == 0:
#			print raw,col
#
#		noiseocc_map.Fill(raw,col,noise)
#
#c1.cd()
#noiseocc_hist_091.Draw()
#
#
#noiseocc_map.Draw("colz")


#
#
#
#
##-------- conversion factor injection -> electrons
#
#injToElectrons=1660./.39	# WARNING: different depending the version of the chip: v2: 1660/0.39 v4:1660/0.25
#
## -------- Fill dictionaries with function AnalyseThresholdScan() in functions.pyc
#
#DIC_AnalyseThresholdScan = AnalyseThresholdScan(FileDataThrScan)
#TDAC_value_dic = DIC_AnalyseThresholdScan[0]
#Threshold_value_dic = DIC_AnalyseThresholdScan[1]
#Sigma_value_dic = DIC_AnalyseThresholdScan[2]
#Chi2_value_dic = DIC_AnalyseThresholdScan[3]
#Scurve_plot_dic = DIC_AnalyseThresholdScan[4]
#Data_pointsX_dic = DIC_AnalyseThresholdScan[5]
#Data_pointsY_dic = DIC_AnalyseThresholdScan[6]
#
## -------- Set histograms and fitting function
#
#myfit = TF1("myfit", "[0]+0.5*TMath::Erf([1]*([2]+x))",-1.,3000.)
#TDAC2D = TH2F("TDAC2D", "TDAC 2D plot;#Row;#Column;Threshold [e]",24,0,24, 36,12,48)
#TDAC1D = TH1F("TDACDist", "TDAC distribution;Threshold [e];nb of pixels", 16,0, 16)
#thresh2D = TH2F("thresh2D", "Threshold 2D plot;#Row;#Column;Threshold [e]",24,0,24, 36,12,48)
#thresh1D = TH1F("threshDist", "Threshold distribution;Threshold [e];nb of pixels", 300,-.1, 4000)
#sigma2D = TH2F("sigma2D", "Sigma 2D plot;#Row;#Column;Sigma [e]",24,0,24, 36,12,48)
#sigma1D = TH1F("sigmaDist", "Sigma distribution;Sigma [e];nb of pixels", 300,-.1, 500)
#chi2_1D = TH1F("chi2Dist", "Chi2 distribution; ;nb of pixels", 300,0, 1)
#chi2_2D = TH2F("chi2_2D", "Chi2 2D plot;#Row;#Column;chi2",24,0,24, 36,12,48)
#correlation = TH2F("correlation", "Threshold vs TDAC;TDAC;Threshold [e]",16,0,16, 50,0,2000)
#eyeDiagram = TH2D("eyeDiagram","S-curves eye diagram;injection [e];Probe",69*2/3,0,0.69*injToElectrons*2/3,128,0,1.5) 
#
##'$\ detected / simulated$'",69*2/3,0,0.69*injToElectrons*2/3,128,0,1.5)
#
## ------- Create Scurves for each pixel
#
#outFile.cd("Graphs")
#
#CanvScurves = TCanvas("S-curves_all")
#CanvScurves.SetGrid()
#CanvScurves.SetFillColor(0)
#CanvScurves.cd()
#
#cnt = 0
#for r in range(24):
#	for c in range(12,48):
#		graph = Scurve_plot_dic["r"+str(r)+"_c"+str(c)+""]
#		Xaxis = graph.GetXaxis()
#		Xaxis.SetLimits(0,2000)
#		graph.GetHistogram().SetMaximum(2)          
#		graph.GetHistogram().SetMinimum(0)
#		outFile.cd("Graphs")
#		graph.Write()
#		
#		# select data
#		
#		if Threshold_value_dic["r"+str(r)+"_c"+str(c)+""] > 800:
#			outFile.cd("High thresh")
#			graph.Write()
#		if Sigma_value_dic["r"+str(r)+"_c"+str(c)+""] > 100:
#			outFile.cd("High sigma")
#			graph.Write()			
#		if Chi2_value_dic["r"+str(r)+"_c"+str(c)+""] > 0.01:
#			outFile.cd("High Chi2")
#			graph.Write()			
#		
#		
#		if cnt ==0:	
#			graph.Draw("AC")
#		else:
#			graph.Draw("C")		# method with Draw("same") also works
#		
#		cnt += 1
#		graphList.append(graph)
#gPad.Update()
#
## ------- Create histograms
#
#for r in range(24):
#	for c in range(12,48):
#		TDAC_value = TDAC_value_dic["r"+str(r)+"_c"+str(c)+""]
#		Threshold_value = Threshold_value_dic["r"+str(r)+"_c"+str(c)+""]
#		Sigma_value = Sigma_value_dic["r"+str(r)+"_c"+str(c)+""]
#		Chi2_value = Chi2_value_dic["r"+str(r)+"_c"+str(c)+""]
#		Data_pointsX_list = Data_pointsX_dic["r"+str(r)+"_c"+str(c)+""]
#		Data_pointsY_list = Data_pointsY_dic["r"+str(r)+"_c"+str(c)+""]
#		
#		# ---- Fill eyeDiagram 
#		
#		for i in range(len(Data_pointsX_list)):
#			Data_pointsX_value = Data_pointsX_list[i]
#			Data_pointsY_value = Data_pointsY_list[i]
#			
#			eyeDiagram.Fill(Data_pointsX_value,Data_pointsY_value)	
#		
#		# ----- Fill other histograms
#		TDAC1D.Fill(TDAC_value)
#		TDAC2D.Fill(r,c,TDAC_value)
#		thresh1D.Fill(Threshold_value)
#		thresh2D.Fill(r,c,Threshold_value)
#		sigma1D.Fill(Sigma_value)
#		sigma2D.Fill(r,c,Sigma_value)
#		chi2_1D.Fill(Chi2_value)
#		chi2_2D.Fill(r,c,Chi2_value)
#		correlation.Fill(TDAC_value,Threshold_value)
#
#outFile.cd()
#
## ------- Customize plots
#
## - eyeDiagram
#
#eyeDiagram.SetAxisRange(0.,2000.,"X")
#eyeDiagram.SetAxisRange(0.,1.4,"Y")
#
## - thresh 2D
#
#thresh2D.SetTitleSize(0.025,"xyz")
#thresh2D.SetTitleOffset(1.3,"z")
#
#thresh2D.SetAxisRange(-1,1500.,"Z")
#thresh2D.GetXaxis().SetNdivisions(32)
#thresh2D.GetYaxis().SetNdivisions(64)
#thresh2D.SetLabelSize(0.02,"X")
#thresh2D.SetLabelSize(0.02,"Y")
#thresh2D.GetZaxis().SetLabelSize(0.025)
#
#x_axis_thresh2D = thresh2D.GetXaxis()
#x_axis_thresh2D.CenterLabels(kTRUE)
#y_axis_thresh2D = thresh2D.GetYaxis()
#y_axis_thresh2D.CenterLabels(kTRUE)
#
## - thresh 1D
#
#thresh1D.SetAxisRange(0.,2000.,"X")
#thresh1D.SetFillColor(38)
#thresh1D.SetFillColorAlpha(38, 0.70)
#
## - Correlation
#
#correlation.GetXaxis().SetNdivisions(16)
#x_axis_correlation = correlation.GetXaxis()
#x_axis_correlation.CenterLabels(kTRUE)
#
## - TDAC 2D
#
#TDAC2D.GetZaxis().SetRangeUser(-1, 16)
#
## - sigma 1D
#
#sigma1D.Rebin()
#
## ------- Create canvas
#
#Canv4Plots = ROOT.TCanvas("Thresh_2D","Thresh_2D",800,600)
#ROOT.SetOwnership(Canv4Plots,False) # TODO: DOESN'T WORK...
#
#Canv4Plots.Divide(2,2)
#Canv4Plots.cd(1)
#gPad.SetGrid()
#gStyle.SetOptStat("rmen")
#thresh1D.Draw()
#gPad.Update()
#
#
#Canv4Plots.cd(2)
#gStyle.SetOptStat("e")
#gStyle.SetStatX(0.9)
#gStyle.SetStatY(0.95)
#
#thresh2D.GetZaxis().SetTitleOffset(1.3)
#gPad.SetRightMargin(0.15)
#
#thresh2D.Draw("colz")
#TDAC2D.Draw("sametext")
#TDAC2D.SetMarkerSize(0.9)
#list_TLine=[]
#for i in range(0,25):
#	line_vert =  TLine(i,12,i,48)
#	line_vert.Draw()
#	list_TLine.append(line_vert)
#
#for i in range(12,49):
#	line_hor =  TLine(0,i,24,i)
#	line_hor.Draw()
#	list_TLine.append(line_hor)
#	
#gPad.Update()
#
#Canv4Plots.cd(3)
#gPad.SetGrid()
#gStyle.SetOptStat("men")
#gStyle.SetStatX(0.9)
#gStyle.SetStatY(0.9)
#correlation.Draw("colztext")
#
#gPad.Update()
#
#Canv4Plots.cd(4)
#gPad.SetGrid()
#gStyle.SetOptStat("n")
#eyeDiagram.Draw("colz")
#
#gPad.Update()
#
## ------- Create threshold 2D map 
#
#
#Canv2DThresh = ROOT.TCanvas("2D_Thresh","Thresh_2D",800,600)
#ROOT.SetOwnership(Canv4Plots,False) # TODO: DOESN'T WORK...
#
#Canv2DThresh.cd()
#gStyle.SetOptStat("e")
#gStyle.SetStatX(0.9)
#gStyle.SetStatY(0.95)
#
#thresh2D.GetZaxis().SetTitleOffset(1.3)
#gPad.SetRightMargin(0.15)
#
#thresh2D.Draw("colz")
#TDAC2D.Draw("sametext")
#TDAC2D.SetMarkerSize(0.9)
#list_TLine=[]
#for i in range(0,25):
#	line_vert =  TLine(i,12,i,48)
#	line_vert.Draw()
#	list_TLine.append(line_vert)
#
#for i in range(12,49):
#	line_hor =  TLine(0,i,24,i)
#	line_hor.Draw()
#	list_TLine.append(line_hor)
#	
#gPad.Update()
#
#Canv2DThresh.Close()
#
## ------- Write plots 
#
#outFile.cd()
#
#thresh1D.Write()
#thresh2D.Write()
#TDAC1D.Write()
#TDAC2D.Write()
#correlation.Write()
#sigma1D.Write()
#sigma2D.Write()
#chi2_1D.Write()
#chi2_2D.Write()
#eyeDiagram.Write()
#CanvScurves.Write()
#Canv4Plots.Write()
#Canv2DThresh

#CanvScurves.Close()

#raw_input("Press Enter to continue...")


