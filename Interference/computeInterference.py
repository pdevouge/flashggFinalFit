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

MH = ROOT.RooRealVar("MH","m_{H}", int(MHLow), int(MHHigh))
MH.setUnit("GeV")
MH.setConstant(True)

if 'p' in opt.width:
  width = opt.width.replace('p','.')
  width = f"({float(width)/100})"
else:
  width = "%s.%s"%(opt.width[0],opt.width[1:])

intfm = InterferenceModel(opt.proc,opt.cat,opt.ext,opt.year,sqrts__,xvar,MH,opt.massPoints,width)


def plotInterference(ifm,_range= 0.1,_binwidth=1.):

  mass = 750

  canv = ROOT.TCanvas()
  canv.SetLeftMargin(0.15)
  range_m, range_p = 600,900
  ifm.MH.setVal(int(mass))
  ifm.MH.setConstant(True)
  frame = ifm.xvar.frame()
  frame.GetYaxis().SetRangeUser(-2, 2)
  # ifm.Pdfs['interference'].plotOn(frame, ROOT.RooFit.Name("pdf"), LineColor=ROOT.kBlue, LineStyle=1, LineWidth=2)
  # ifm.Functions['I_re'].plotOn(frame, ROOT.RooFit.Name("Re"), LineColor=ROOT.kRed, LineStyle=1, LineWidth=2)
  ifm.Functions['I_im'].plotOn(frame, ROOT.RooFit.Name("Im"), LineColor=ROOT.kGreen, LineStyle=1, LineWidth=2)
  frame.Draw()

  canv.SaveAs("./interference_model.pdf")


plotInterference(intfm)