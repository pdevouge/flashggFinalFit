import ROOT
import numpy as np
import pandas as pd
import pickle
import math
import os, sys
import json
from optparse import OptionParser
import glob
import re
from collections import OrderedDict as od

from interferenceModel import *
from plottingTools import *

print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG INTERFERENCE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
def leave():
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG INTERFERENCE (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
  exit()

def get_options():
  parser = OptionParser()
  parser.add_option("--xvar", dest='xvar', default='CMS_hgg_mass', help="Observable to fit")
  parser.add_option("--inputWSDir", dest='inputWSDir', default='', help="Input flashgg WS directory")
  parser.add_option("--ext", dest='ext', default='', help="Extension")
  parser.add_option("--proc", dest='proc', default='', help="Signal process")
  parser.add_option("--cat", dest='cat', default='', help="RECO category")
  parser.add_option("--year", dest='year', default='2016', help="Year")
  parser.add_option('--width', dest='width', default='001', help="Signal width")
  parser.add_option('--massPoints', dest='massPoints', default='250,300,350,400,450,500', help="Mass points to fit")
  parser.add_option('--minMass', dest='minMass', default='200', help="Mass range lower boundary")
  parser.add_option('--maxMass', dest='maxMass', default='600', help="Mass range upper boundary")
  parser.add_option('--doPlots', dest='doPlots', default=False, action="store_true", help="Save an S / B / S+B+I overlay next to the output workspace")
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

print(" --> Building interference model for (proc,cat,year) = (%s,%s,%s), width %s, range [%s, %s]"
      %(opt.proc,opt.cat,opt.year,width,MHLow,MHHigh))

intfm = InterferenceModel(opt.proc,opt.cat,opt.ext,opt.year,sqrts__,xvar,MHLow,MHHigh,opt.massPoints,width)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE: to output workspace
foutDir = "%s/results/outdir_%s/computeIntf/output"%(iwd__,opt.ext)
foutName = "%s/CMS-HGG_intfm_%s_%s_%s_%s.root"%(foutDir,opt.ext,opt.proc,opt.year,opt.cat)
print("\n --> Saving output workspace to file: %s"%foutName)
if not os.path.isdir(foutDir): os.system("mkdir %s"%foutDir)
fout = ROOT.TFile(foutName,"RECREATE")
outWS = ROOT.RooWorkspace("%s"%(intfWSName__),"%s"%(intfWSName__))
intfm.save(outWS)
outWS.Write()
fout.Close()


if opt.doPlots:
  plotDir = "%s/results/outdir_%s/computeIntf/plots"%(iwd__,opt.ext)
  if not os.path.isdir(plotDir): os.system("mkdir %s"%plotDir)
  tag = "%s_%s_%s_%s"%(opt.ext,opt.proc,opt.year,opt.cat)
  lumi = float(lumiMap[opt.year])*1000
  plotInterferenceTemplates(intfm, float(MHNominal), lumi, plotDir, _extension=tag, _logy=True)
  plotInterferenceTemplates(intfm, float(MHNominal), lumi, plotDir, _extension=tag, _logy=False)
  plotInterferenceTerm(intfm, float(MHNominal), lumi, plotDir, _extension=tag)
  plotInterferencePhase(intfm, float(MHNominal), plotDir, _extension=tag)
