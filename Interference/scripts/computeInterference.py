import ROOT
import pandas as pd
import pickle
import math
import os, sys
import json
from optparse import OptionParser
import glob
import re
from collections import OrderedDict as od

from tools.interferenceModel import *

print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG INTERFERENCE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
def leave():
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG INTERFERENCE (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
  exit()

def get_options():
  parser = OptionParser()
  parser.add_option("--xvar", dest='xvar', default='CMS_hgg_mass', help="Observable to fit")
  parser.add_option("--inputWSDir", dest='inputWSDir', default='', help="Input flashgg WS directory")
  parser.add_option("--inputSigM", dest='inputSigM', default='', help="Input signal model")
  parser.add_option("--ext", dest='ext', default='', help="Extension")
  parser.add_option("--proc", dest='proc', default='', help="Signal process")
  parser.add_option("--cat", dest='cat', default='', help="RECO category")
  parser.add_option("--year", dest='year', default='2016', help="Year")
  parser.add_option('--width', dest='width', default='001', help="Signal width")
  parser.add_option('--massPoints', dest='massPoints', default='250,300,350,400,450,500', help="Mass points to fit")
  parser.add_option('--minMass', dest='minMass', default='200', help="Mass range lower boundary")
  parser.add_option('--maxMass', dest='maxMass', default='600', help="Mass range upper boundary")
  return parser.parse_args()
(opt,args) = get_options()

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

w_indicator = 'W' if 'p' in opt.width else 'kMpl'
lowW = '001' if opt.proc=='rsg' else '0p014'
lowW_str = w_indicator + lowW
nomW_str = w_indicator + opt.width
MHLow = opt.minMass
MHHigh = opt.maxMass
masses = opt.massPoints.split(",")
MHNominal = masses[len(masses)//2]

nominalWSFileName = glob.glob("%s/output*M%s_%s*%s.root"%(opt.inputWSDir,MHNominal,nomW_str,opt.proc))[0]
f0 = ROOT.TFile(nominalWSFileName,"read")
inputWS0 = f0.Get(inputWSName__)
xvar = inputWS0.var(opt.xvar)
xvar.setRange(int(MHLow), int(MHHigh))
xvarFit = xvar.Clone()

if 'p' in opt.width:
  width = opt.width.replace('p','.')
  width = f"({float(width)/100})"
else:
  width = "%s.%s"%(opt.width[0],opt.width[1:])

intfm = InterferenceModel(opt.proc,opt.cat,opt.ext,opt.year,sqrts__,xvar,opt.massPoints,width)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE: to output workspace
foutDir = "%s/outdir_%s/computeIntf/output"%(iwd__,opt.ext)
foutName = "%s/outdir_%s/computeIntf/output/CMS-HGG_intfm_%s_%s.root"%(iwd__,opt.ext,opt.cat,opt.year)
print("\n --> Saving output workspace to file: %s"%foutName)
if not os.path.isdir(foutDir): os.system("mkdir %s"%foutDir)
fout = ROOT.TFile(foutName,"RECREATE")
outWS = ROOT.RooWorkspace("%s"%(intfWSName__),"%s"%(intfWSName__))
intfm.save(outWS)
outWS.Write()
fout.Close()


def plotInterference(ifm,_range= 0.1,_binwidth=1.):

  mass = 500
  ifm.Vars['MH'].setVal(mass)

  sig_pdf = ifm.Pdfs["Msig"]
  sig_norm = ifm.Functions["Msig_norm"]

  bkg_pdf = ifm.Pdfs["Mbkg"]
  bkg_norm = ifm.Functions["Mbkg_norm"]

  sbi_pdf = ifm.Pdfs["SBI"]
  sbi_norm = ifm.Functions["SBI_norm"]

  graph_sig = ROOT.TGraph()
  graph_bkg = ROOT.TGraph()
  graph_sbi = ROOT.TGraph()

  ymin = 0
  ymax = 0
  for i, x in enumerate(range(400, 1201)):
      ifm.xvar.setVal(x)
      y_sig = sig_pdf.getVal(ROOT.RooArgSet(ifm.xvar)) * sig_norm.getVal()
      graph_sig.SetPoint(i, x, y_sig)

      y_bkg = bkg_pdf.getVal(ROOT.RooArgSet(ifm.xvar)) * bkg_norm.getVal()
      graph_bkg.SetPoint(i, x, y_bkg)

      y_sbi = sbi_pdf.getVal(ROOT.RooArgSet(ifm.xvar)) * sbi_norm.getVal()
      graph_sbi.SetPoint(i, x, y_sbi)
      ymax = np.max([ymax,y_sig,y_bkg,y_sbi])

  c = ROOT.TCanvas()
  graph_sig.SetLineColor(ROOT.kGreen)
  graph_sig.SetLineWidth(2)
  graph_bkg.SetLineColor(ROOT.kRed)
  graph_bkg.SetLineWidth(2)
  graph_sbi.SetLineColor(ROOT.kBlue)
  graph_sbi.SetLineWidth(2)
  graph_sbi.SetLineStyle(2)

  graph_sig.Draw("AL")
  graph_sig.SetMaximum(1.1 * ymax)   # 10% headroom
  graph_sig.SetMinimum(0.)
  graph_bkg.Draw("L SAME")
  graph_sbi.Draw("L SAME")

  leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)

  leg.AddEntry(graph_sbi, "S+B+I", "l")
  leg.AddEntry(graph_bkg, "Background", "l")
  leg.AddEntry(graph_sig, "Signal", "l")

  leg.Draw()

  graph_sig.GetXaxis().SetLimits(400, 1200)
  graph_sig.SetTitle('Interference process')
  c.Update()

  c.SaveAs("./interference_model.pdf")


plotInterference(intfm)