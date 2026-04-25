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

from commonTools import *
from commonObjects import *
from signalTools import *
from replacementMap import globalReplacementMap
from XSBRMap import *
from simultaneousFit import *
from interferenceModel import *
from finalModel import *
from plottingTools import *

print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG SIGNAL FITTER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
def leave():
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG SIGNAL FITTER (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
  exit()

def get_options():
  parser = OptionParser()
  parser.add_option("--xvar", dest='xvar', default='CMS_hgg_mass', help="Observable to fit")
  parser.add_option("--inputWSDir", dest='inputWSDir', default='', help="Input flashgg WS directory")
  parser.add_option("--ext", dest='ext', default='', help="Extension")
  parser.add_option("--proc", dest='proc', default='', help="Signal process")
  parser.add_option("--cat", dest='cat', default='', help="RECO category")
  parser.add_option("--year", dest='year', default='2016', help="Year")
  parser.add_option("--analysis", dest='analysis', default='STXS', help="Analysis handle: used to specify replacement map and XS*BR normalisations")
  parser.add_option('--width', dest='width', default='001', help="Signal width")
  parser.add_option('--massPoints', dest='massPoints', default='250,300,350,400,450,500', help="Mass points to fit")
  parser.add_option('--minMass', dest='minMass', default='200', help="Mass range lower boundary")
  parser.add_option('--maxMass', dest='maxMass', default='600', help="Mass range upper boundary")
  parser.add_option('--skipBeamspotReweigh', dest='skipBeamspotReweigh', default=True, action="store_true", help="Skip beamspot reweigh to match beamspot distribution in data")
  parser.add_option('--doPlots', dest='doPlots', default=False, action="store_true", help="Produce Signal Fitting plots")
  parser.add_option("--doVoigtian", dest='doVoigtian', default=False, action="store_true", help="Use Voigtians instead of Gaussians for signal models with Higgs width as parameter")
  parser.add_option('--skipResolutionModel', dest='skipResolutionModel', default=False, action="store_true", help="Skip detector resolution model using DCB")
  parser.add_option('--skipMC', dest='skipMC', default=False, action="store_true", help="Skip MC reading (NB: The MC are needed for the resolution model!)")
  parser.add_option("--useInterpolation", dest='useInterpolation', default=False, action="store_true", help="Use signal interpolation instead of analytical form")
  parser.add_option("--useDCB", dest='useDCB', default=False, action="store_true", help="Use DCB in signal interpolation approach")
  parser.add_option("--useDiagonalProcForShape", dest='useDiagonalProcForShape', default=False, action="store_true", help="Use shape of diagonal process, keeping normalisation (requires diagonal mapping produced by getDiagProc script)")
  parser.add_option('--skipVertexScenarioSplit', dest='skipVertexScenarioSplit', default=True, action="store_true", help="Skip vertex scenario split")
  parser.add_option('--skipZeroes', dest='skipZeroes', default=False, action="store_true", help="Skip proc x cat is numEntries = 0., or sumEntries < 0.")
  # For systematics
  parser.add_option('--skipSystematics', dest='skipSystematics', default=False, action="store_true", help="Skip shape systematics in signal model")
  parser.add_option('--useDiagonalProcForSyst', dest='useDiagonalProcForSyst', default=False, action="store_true", help="Use diagonal process for systematics (requires diagonal mapping produced by getDiagProc script)")
  parser.add_option("--scales", dest='scales', default='', help="Photon shape systematics: scales")
  parser.add_option("--scalesCorr", dest='scalesCorr', default='', help='Photon shape systematics: scalesCorr')
  parser.add_option("--scalesGlobal", dest='scalesGlobal', default='', help='Photon shape systematics: scalesGlobal')
  parser.add_option("--smears", dest='smears', default='', help='Photon shape systematics: smears')
  # Parameter values
  parser.add_option('--replacementThreshold', dest='replacementThreshold', default=100, type='int', help="Nevent threshold to trigger replacement dataset")
  parser.add_option('--beamspotWidthData', dest='beamspotWidthData', default=3.5, type='float', help="Width of beamspot in data [cm]")
  parser.add_option('--beamspotWidthMC', dest='beamspotWidthMC', default=3.7, type='float', help="Width of beamspot in MC [cm]")
  parser.add_option('--MHPolyOrder', dest='MHPolyOrder', default=0, type='int', help="Order of polynomial for MH dependence")
  parser.add_option('--nBins', dest='nBins', default=400, type='int', help="Number of bins for fit")
  # Minimizer options
  parser.add_option('--minimizerMethod', dest='minimizerMethod', default='TNC', help="(Scipy) Minimizer method")
  parser.add_option('--minimizerTolerance', dest='minimizerTolerance', default=1e-8, type='float', help="(Scipy) Minimizer toleranve")
  return parser.parse_args()
(opt,args) = get_options()

if opt.skipResolutionModel:
  opt.skipMC = True

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

w_indicator = 'W' if 'p' in opt.width else 'kMpl'
lowW = '001' if opt.proc=='rsg' else '0p014'
lowW_str = w_indicator + lowW
nomW_str = w_indicator + opt.width
MHLow = opt.minMass
MHHigh = opt.maxMass
print("Width", opt.width)
print("MHLow", opt.minMass)
print("MHHigh", opt.maxMass)

masses = opt.massPoints.split(",")
MHNominal = masses[len(masses)//2] # Use middle mass point as nominal mass. TOFIX

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SETUP: signal fit
print(" --> Running fit for (proc,cat) = (%s,%s)"%(opt.proc,opt.cat))
if( len(opt.massPoints.split(",")) == 1 )&( opt.MHPolyOrder > 0 ):
  print(" --> [WARNING] Attempting to fit polynomials of O(MH^%g) for single mass point. Setting order to 0"%opt.MHPolyOrder)
  opt.MHPolyOrder=0

# Add stopwatch function

# Load replacement map
if opt.analysis not in globalReplacementMap:
  print(" --> [ERROR] replacement map does not exist for analysis (%s). Please add to tools/replacementMap.py"%opt.analysis)
  leave()
else: rMap = globalReplacementMap[opt.analysis]

# Load XSBR map
if opt.analysis not in globalXSBRMap:
  print(" --> [ERROR] XS * BR map does not exist for analysis (%s). Please add to tools/XSBRMap.py"%opt.analysis)
  leave()
else: xsbrMap = globalXSBRMap[opt.analysis]

if not opt.skipMC:
  # Load RooRealVars from workspace
  nominalWSFileName = glob.glob("%s/output*M%s_%s*%s.root"%(opt.inputWSDir,MHNominal,nomW_str,opt.proc))[0]
  f0 = ROOT.TFile(nominalWSFileName,"read")
  inputWS0 = f0.Get(inputWSName__)
  xvar = inputWS0.var(opt.xvar)
  xvar.setRange(int(MHLow), int(MHHigh))
  xvarFit = xvar.Clone()
  print(xvar)
  dZ = inputWS0.var("dZ")
  true_mass = inputWS0.var("true_mass")
  true_mass.setConstant(True)
  reduced_mass = inputWS0.var("reduced_mass")
  reduced_mass.setRange(-0.1, 0.1)
  aset = ROOT.RooArgSet(xvar,dZ,true_mass,reduced_mass)
  f0.Close()
else:
  xvar = ROOT.RooRealVar("CMS_hgg_mass","CMS_hgg_mass",float(MHNominal),float(MHLow),float(MHHigh))
  xvar.setBins((int(MHHigh)-int(MHLow))*2)
  xvar.setBins(10000, "cache")
  xvarFit = xvar.Clone()
  true_mass = ROOT.RooRealVar("true_mass","true_mass",float(MHNominal),float(MHLow),float(MHHigh))
  true_mass.setBins(int(MHHigh)-int(MHLow))
  true_mass.setConstant(True)
  reduced_mass = ROOT.RooRealVar("reduced_mass","reduced_mass",0., -0.1, 0.1)
  reduced_mass.setBins(100)
  dZ = ROOT.RooRealVar("dZ","dZ",0.,-20.,20.)
  dZ.setBins(40)
  aset = ROOT.RooArgSet(xvar,dZ,true_mass,reduced_mass)

# Create MH var
MH = ROOT.RooRealVar("MH","m_{H}", int(MHLow), int(MHHigh))
MH.setUnit("GeV")
MH.setConstant(True)

G0 = ROOT.RooRealVar("G0","Gamma_{0}", 0, 10000)
G0.setConstant(True)

if opt.skipZeroes:
  # Extract nominal mass dataset and see if entries == 0
  WSFileName = glob.glob("%s/output*M%s_%s*%s.root"%(opt.inputWSDir,MHNominal,nomW_str,opt.proc))[0]
  f = ROOT.TFile(WSFileName,"read")
  inputWS = f.Get(inputWSName__)
  d = reduceDataset(inputWS.data("%s_%s_%s_%s_%s"%(procToData(opt.proc.split("_")[0]),MHNominal,opt.width,sqrts__,opt.cat)),aset)
  if( d.numEntries() == 0. )|( d.sumEntries <= 0. ):
    print(" --> (%s,%s) has zero events. Will not construct signal model"%(opt.proc,opt.cat))
    exit()
  inputWS.Delete()
  f.Close()

# Define proc x cat with which to extract shape: if skipVertexScenarioSplit label all events as "RV"
procRVFit, catRVFit = opt.proc, opt.cat
if opt.skipVertexScenarioSplit:
  print(" --> Skipping vertex scenario split")
else:
  procWVFit, catWVFit = opt.proc, opt.cat

# Options for using diagonal process from getDiagProc output json
if opt.useDiagonalProcForShape:
  if not os.path.exists("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(swd__,opt.ext)):
    print(" --> [ERROR] Diagonal process json from getDiagProc does not exist. Using nominal proc x cat for shape")
  else:
    with open("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(swd__,opt.ext),"r") as jf: dproc = json.load(jf)
    procRVFit = dproc[opt.cat]
    print(" --> Using diagonal proc (%s,%s) for shape"%(procRVFit,opt.cat))
    if not opt.skipVertexScenarioSplit: procWVFit = dproc[opt.cat]

# Process for syst
procSyst = opt.proc
if opt.useDiagonalProcForSyst:
  if not os.path.exists("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(swd__,opt.ext)):
    print(" --> [ERROR] Diagonal process json from getDiagProc does not exist. Using nominal proc x cat for systematics")
  else:
    with open("%s/outdir_%s/getDiagProc/json/diagonal_process.json"%(swd__,opt.ext),"r") as jf: dproc = json.load(jf)
    procSyst = dproc[opt.cat]
    print(" --> Using diagonal proc (%s,%s) for systematics"%(procSyst,opt.cat))

# Define process with which to extract normalisation: nominal
procNorm, catNorm = opt.proc, opt.cat

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTRACT DATASETS TO FIT (for each mass point)
nominalDatasets = od()
# For RV (or if skipping vertex scenario split)
datasetRVForFit = od()
datasetRVForFit['low_w'] = od()
datasetRVForFit['nom_w'] = od()

if not opt.skipVertexScenarioSplit:
  datasetWVForFit = od()
  datasetWVForFit['low_w'] = od()
  datasetWVForFit['nom_w'] = od()

if not opt.skipMC:
  for mp in opt.massPoints.split(","):
    # Load low width samples for true lineshape description
    WSFileName = glob.glob("%s/output*M%s_%s*%s.root"%(opt.inputWSDir,mp,lowW_str,procRVFit))[0]
    f = ROOT.TFile(WSFileName,"read")
    inputWS = f.Get(inputWSName__)
    Roodata = inputWS.data("%s_%s_%s_%s_%s"%(procToData(procRVFit.split("_")[0]),mp,lowW,sqrts__,catRVFit))
    d = reduceDataset(Roodata,aset)
    if opt.skipVertexScenarioSplit: datasetRVForFit['low_w'][mp] = d
    else: datasetRVForFit['low_w'][mp] = splitRVWV(d,aset,mode="RV")
    inputWS.Delete()
    f.Close()

    WSFileName = glob.glob("%s/output*M%s_%s*%s.root"%(opt.inputWSDir,mp,nomW_str,procRVFit))[0]
    f = ROOT.TFile(WSFileName,"read")
    inputWS = f.Get(inputWSName__)
    d = reduceDataset(inputWS.data("%s_%s_%s_%s_%s"%(procToData(procRVFit.split("_")[0]),mp,opt.width,sqrts__,catRVFit)),aset)
    nominalDatasets[mp] = d.Clone()
    if opt.skipVertexScenarioSplit: datasetRVForFit['nom_w'][mp] = d
    else: datasetRVForFit['nom_w'][mp] = splitRVWV(d,aset,mode="RV")
    inputWS.Delete()
    f.Close()

  # Check if nominal yield > threshold (or if +ve sum of weights). If not then use replacement proc x cat
  if( datasetRVForFit['nom_w'][MHNominal].numEntries() < opt.replacementThreshold  )|( datasetRVForFit['nom_w'][MHNominal].sumEntries() < 0. ):
    nominal_numEntries = datasetRVForFit['nom_w'][MHNominal].numEntries()
    procReplacementFit, catReplacementFit = rMap['procRVMap'][opt.cat], rMap['catRVMap'][opt.cat]
    for mp in opt.massPoints.split(","):
      WSFileName = glob.glob("%s/output*M%s_%s*%s.root"%(opt.inputWSDir,mp,nomW_str,procReplacementFit))[0]
      f = ROOT.TFile(WSFileName,"read")
      inputWS = f.Get(inputWSName__)
      d = reduceDataset(inputWS.data("%s_%s_%s_%s_%s"%(procToData(procReplacementFit.split("_")[0]),mp,opt.width,sqrts__,catReplacementFit)),aset)
      if opt.skipVertexScenarioSplit: datasetRVForFit['nom_w'][mp] = d
      else: datasetRVForFit['nom_w'][mp] = splitRVWV(d,aset,mode="RV")
      inputWS.Delete()
      f.Close()

    # Check if replacement dataset has too few entries: if so throw error
    if( datasetRVForFit['nom_w'][MHNominal].numEntries() < opt.replacementThreshold )|( datasetRVForFit['nom_w'][MHNominal].sumEntries() < 0. ):
      print(" --> [ERROR] replacement dataset (%s,%s) has too few entries (%g < %g)"%(procReplacementFit,catReplacementFit,datasetRVForFit[MHNominal].numEntries(),opt.replacementThreshold))
      sys.exit(1)

    else:
      procRVFit, catRVFit = procReplacementFit, catReplacementFit
      if opt.skipVertexScenarioSplit:
        print(" --> Too few entries in nominal dataset (%g < %g). Using replacement (proc,cat) = (%s,%s) for extracting shape"%(nominal_numEntries,opt.replacementThreshold,procRVFit,catRVFit))
        for mp in opt.massPoints.split(","):
          print("     * MH = %s GeV: numEntries = %g, sumEntries = %.6f"%(mp,datasetRVForFit['nom_w'][mp].numEntries(),datasetRVForFit['nom_w'][mp].sumEntries()))
      else:
        print(" --> RV: Too few entries in nominal dataset (%g < %g). Using replacement (proc,cat) = (%s,%s) for extracting shape"%(nominal_numEntries,opt.replacementThreshold,procRVFit,catRVFit))
        for mp in opt.massPoints.split(","):
          print("     * MH = %s: numEntries = %g, sumEntries = %.6f"%(mp,datasetRVForFit['nom_w'][mp].numEntries(),datasetRVForFit['nom_w'][mp].sumEntries()))

  else:
    if opt.skipVertexScenarioSplit:
      print(" --> Using (proc,cat) = (%s,%s) for extracting shape"%(procRVFit,catRVFit))
      for mp in opt.massPoints.split(","):
        print("     * MH = %s: numEntries = %g, sumEntries = %.6f"%(mp,datasetRVForFit['nom_w'][mp].numEntries(),datasetRVForFit['nom_w'][mp].sumEntries()))
    else:
      print(" --> RV: Using (proc,cat) = (%s,%s) for extracting shape"%(procRVFit,catRVFit))
      for mp in opt.massPoints.split(","):
        print("     * MH = %s: numEntries = %g, sumEntries = %.6f"%(mp,datasetRVForFit['nom_w'][mp].numEntries(),datasetRVForFit['nom_w'][mp].sumEntries()))

  # Repeat for WV scenario
  if not opt.skipVertexScenarioSplit:
    for mp in opt.massPoints.split(","):
      # Load low width samples for true lineshape description
      WSFileName = glob.glob("%s/output*M%s_%s*%s.root"%(opt.inputWSDir,mp,lowW_str,procWVFit))[0]
      f = ROOT.TFile(WSFileName,"read")
      inputWS = f.Get(inputWSName__)
      d = reduceDataset(inputWS.data("%s_%s_%s_%s_%s"%(procToData(procWVFit.split("_")[0]),mp,lowW,sqrts__,catWVFit)),aset)
      datasetWVForFit['low_w'][mp] = splitRVWV(d,aset,mode="WV")
      inputWS.Delete()
      f.Close()

      WSFileName = glob.glob("%s/output*M%s_%s*%s.root"%(opt.inputWSDir,mp,nomW_str,procWVFit))[0]
      f = ROOT.TFile(WSFileName,"read")
      inputWS = f.Get(inputWSName__)
      d = reduceDataset(inputWS.data("%s_%s_%s_%s_%s"%(procToData(procWVFit.split("_")[0]),mp,opt.width,sqrts__,catWVFit)),aset)
      datasetWVForFit['nom_w'][mp] = splitRVWV(d,aset,mode="WV")
      inputWS.Delete()
      f.Close()

    # Check nominal mass dataset
    if( datasetWVForFit['nom_w'][MHNominal].numEntries() < opt.replacementThreshold  )|( datasetWVForFit['nom_w'][MHNominal].sumEntries() < 0. ):
      nominal_numEntries = datasetWVForFit['nom_w'][MHNominal].numEntries()
      procReplacementFit, catReplacementFit = rMap['procWV'], rMap['catWV']
      for mp in opt.massPoints.split(","):
        WSFileName = glob.glob("%s/output*M%s_%s*%s.root"%(opt.inputWSDir,mp,nomW_str,procReplacementFit))[0]
        f = ROOT.TFile(WSFileName,"read")
        inputWS = f.Get(inputWSName__)
        d = reduceDataset(inputWS.data("%s_%s_%s_%s_%s"%(procToData(procReplacementFit.split("_")[0]),mp,opt.width,sqrts__,catReplacementFit)),aset)
        datasetWVForFit['nom_w'][mp] = splitRVWV(d,aset,mode="WV")
        inputWS.Delete()
        f.Close()
      # Check if replacement dataset has too few entries: if so throw error
      if( datasetWVForFit['nom_w'][MHNominal].numEntries() < opt.replacementThreshold )|( datasetWVForFit['nom_w'][MHNominal].sumEntries() < 0. ):
        print(" --> [ERROR] replacement dataset (%s,%s) has too few entries (%g < %g)"%(procReplacementFit,catReplacementFit,datasetWVForFit['nom_w'][MHNominal].numEntries,opt.replacementThreshold))
        sys.exit(1)
      else:
        procWVFit, catWVFit = procReplacementFit, catReplacementFit
        print(" --> WV: Too few entries in nominal dataset (%g < %g). Using replacement (proc,cat) = (%s,%s) for extracting shape"%(nominal_numEntries,opt.replacementThreshold,procWVFit,catWVFit))
        for mp in opt.massPoints.split(","):
          print("     * MH = %s: numEntries = %g, sumEntries = %.6f"%(mp,datasetWVForFit['nom_w'][mp].numEntries(),datasetWVForFit['nom_w'][mp].sumEntries()))
    else:
      print(" --> WV: Using (proc,cat) = (%s,%s) for extracting shape"%(procWVFit,catRVFit))
      for mp in opt.massPoints.split(","):
        print("     * MH = %s: numEntries = %g, sumEntries = %.6f"%(mp,datasetWVForFit['nom_w'][mp].numEntries(),datasetWVForFit['nom_w'][mp].sumEntries()))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BEAMSPOT REWEIGHT
if not opt.skipBeamspotReweigh:
  # Datasets for fit
  for mp,d in datasetRVForFit['nom_w'].items():
    drw = beamspotReweigh(datasetRVForFit['nom_w'][mp],opt.beamspotWidthData,opt.beamspotWidthMC,xvar,dZ,_x=opt.xvar)
    datasetRVForFit['nom_w'][mp] = drw
  if not opt.skipVertexScenarioSplit:
    for mp,d in datasetWVForFit['nom_w'].items():
      drw = beamspotReweigh(datasetWVForFit['nom_w'][mp],opt.beamspotWidthData,opt.beamspotWidthMC,xvar,dZ,_x=opt.xvar)
      datasetWVForFit['nom_w'][mp] = drw
    print(" --> Beamspot reweigh: RV(sumEntries) = %.6f, WV(sumEntries) = %.6f"%(datasetRVForFit['nom_w'][mp].sumEntries(),datasetWVForFit['nom_w'][mp].sumEntries()))
  else:
    print(" --> Beamspot reweigh: sumEntries = %.6f"%datasetRVForFit['nom_w'][mp].sumEntries())

  # Nominal datasets for saving to output Workspace: preserve norm for eff * acc calculation
  for mp,d in nominalDatasets.items():
    drw = beamspotReweigh(d,opt.beamspotWidthData,opt.beamspotWidthMC,xvar,dZ,_x=opt.xvar,preserveNorm=True)
    nominalDatasets[mp] = drw

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# If using nGaussian fit then extract nGaussians from fTest json file
if opt.useInterpolation:
  if not opt.useDCB:
    with open("%s/outdir_%s/fTest/json/nGauss_%s.json"%(swd__,opt.ext,catRVFit)) as jf: ngauss = json.load(jf)
    nRV = int(ngauss["%s__%s"%(procRVFit,catRVFit)]['nRV'])
    if opt.skipVertexScenarioSplit: print(" --> Fitting function: convolution of nGaussians (%g)"%nRV)
    else:
      with open("%s/outdir_%s/fTest/json/nGauss_%s.json"%(swd__,opt.ext,catWVFit)) as jf: ngauss = json.load(jf)
      nWV = int(ngauss["%s__%s"%(procWVFit,catWVFit)]['nWV'])
      print(" --> Fitting function: convolution of nGaussians (RV=%g,WV=%g)"%(nRV,nWV))
  else:
    print(" --> Fitting function: DCB + 1 Gaussian")

  if opt.doVoigtian:
    print(" --> Will add natural Higgs width as parameter in Pdf (Gaussians -> Voigtians)")
else:
  print(" --> Fitting function: Analytical form => DCB * BW")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIT: simultaneous signal fit (ssf)
ssfMap = od()
name = "Total" if opt.skipVertexScenarioSplit else "RV"
if 'p' in opt.width:
  width = opt.width.replace('p','.')
  width = f"({float(width)/100})"
else:
  width = "%s.%s"%(opt.width[0],opt.width[1:])
ssfRV = SimultaneousFit(name,opt.proc,opt.cat,datasetRVForFit,xvar.Clone(),true_mass.Clone(),reduced_mass.Clone(),MH,MHLow,MHHigh,
                        width,
                        opt.massPoints,opt.nBins,opt.MHPolyOrder,opt.minimizerMethod,opt.minimizerTolerance)
if opt.useInterpolation:
  if opt.useDCB: ssfRV.buildDCBplusGaussian()
  else: ssfRV.buildNGaussians(nRV)
  ssfRV.runFit()
  ssfRV.buildSplines()
else:
  ssfRV.buildTrueLineshape(decay='hgg', xsec='sm')
  if not opt.skipResolutionModel:
    ssfRV.buildResoModel()
    ssfRV.buildAnalytical()

ssfMap[name] = ssfRV

if not opt.skipVertexScenarioSplit:
  name = "WV"
  ssfWV = SimultaneousFit(name,opt.proc,opt.cat,datasetWVForFit,xvar.Clone(),true_mass.Clone(),reduced_mass.Clone(),MH,MHLow,MHHigh,width,opt.massPoints,opt.nBins,opt.MHPolyOrder,opt.minimizerMethod,opt.minimizerTolerance)
  if opt.useInterpolation:
    if opt.useDCB: ssfWV.buildDCBplusGaussian()
    else: ssfWV.buildNGaussians(nRV)
    ssfWV.runFit()
    ssfWV.buildSplines()
  else:
    ssfWV.buildTrueLineshape(decay='hgg', xsec='sm')
    if not opt.skipResolutionModel:
      ssfRV.buildResoModel()
      ssfRV.buildAnalytical()

  ssfMap[name] = ssfWV

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INTERFERENCE MODEL: construction
print("\n --> Constructing interference model")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FINAL MODEL: construction
if opt.skipResolutionModel:
  print("\n --> We cannot construct the final model without the resolution model...")
else:
  print("\n --> Constructing final model")
  fm = FinalModel(ssfMap,opt.proc,opt.cat,opt.ext,opt.year,sqrts__,nominalDatasets,xvar,MH,MHNominal,MHLow,MHHigh,opt.massPoints,G0,nomW_str,xsbrMap,procSyst,opt.scales,opt.scalesCorr,opt.scalesGlobal,opt.smears,opt.doVoigtian,opt.useDCB,opt.skipVertexScenarioSplit,opt.skipSystematics)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # SAVE: to output workspace
  foutDir = "%s/outdir_%s/signalFit/output"%(swd__,opt.ext)
  foutName = "%s/outdir_%s/signalFit/output/CMS-HGG_sigfit_%s_%s_%s_%s.root"%(swd__,opt.ext,opt.ext,opt.proc,opt.year,opt.cat)
  print("\n --> Saving output workspace to file: %s"%foutName)
  if not os.path.isdir(foutDir): os.system("mkdir %s"%foutDir)
  fout = ROOT.TFile(foutName,"RECREATE")
  outWS = ROOT.RooWorkspace("%s_%s"%(outputWSName__,sqrts__),"%s_%s"%(outputWSName__,sqrts__))
  fm.save(outWS)
  outWS.Write()
  fout.Close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOTTING
if opt.doPlots:
  print("\n --> Making plots...")
  if not os.path.isdir("%s/outdir_%s/signalFit/Plots"%(swd__,opt.ext)): os.system("mkdir %s/outdir_%s/signalFit/Plots"%(swd__,opt.ext))
  if not opt.useInterpolation:
    if not opt.skipResolutionModel:
      if not os.path.isdir("%s/outdir_%s/signalFit/Plots/resolutionDCB"%(swd__,opt.ext)): os.system("mkdir %s/outdir_%s/signalFit/Plots/resolutionDCB"%(swd__,opt.ext))
      plotIndividualDCB(ssfRV,_outdir="%s/outdir_%s/signalFit/Plots/resolutionDCB"%(swd__,opt.ext), _skipMC=opt.skipMC)
      plotIndividualDCB(ssfRV,_outdir="%s/outdir_%s/signalFit/Plots/resolutionDCB"%(swd__,opt.ext), _from_formulas=True, _skipMC=opt.skipMC)
      plotDCBParameters(ssfRV,_outdir="%s/outdir_%s/signalFit/Plots/resolutionDCB"%(swd__,opt.ext))
    if not os.path.isdir("%s/outdir_%s/signalFit/Plots/trueLineshapeBW"%(swd__,opt.ext)): os.system("mkdir %s/outdir_%s/signalFit/Plots/trueLineshapeBW"%(swd__,opt.ext))
    truemass_range = 0.001 if (opt.width == "001" or opt.width == "0p014") else 0.2
    truemass_nbins = 150 if (opt.width == "001" or opt.width == "0p014") else 100
    plotTrueLineshape(ssfRV,_outdir="%s/outdir_%s/signalFit/Plots/trueLineshapeBW"%(swd__,opt.ext),_range=truemass_range,_nbins=truemass_nbins, _skipMC=opt.skipMC)
    # Uncomment the following to plot the comparison to internally produced Pythia samples
    # plotPythiaComparison(ssfRV,_outdir="%s/outdir_%s/signalFit/Plots/trueLineshapeBW"%(swd__,opt.ext),_range=truemass_range,_nbins=truemass_nbins, _skipMC=opt.skipMC, proc=opt.proc)
    if not opt.skipResolutionModel:
      if not os.path.isdir("%s/outdir_%s/signalFit/Plots/analyticalModel"%(swd__,opt.ext)): os.system("mkdir %s/outdir_%s/signalFit/Plots/analyticalModel"%(swd__,opt.ext))
      recomass_range = 0.1 if (opt.width == "001" or opt.width == "0p014") else 0.2
      recomass_bwidth = 1 if (opt.width == "001" or opt.width == "0p014") else 4
      plotAnalyticalModel(ssfRV,_outdir="%s/outdir_%s/signalFit/Plots/analyticalModel"%(swd__,opt.ext),_range=recomass_range,_binwidth=recomass_bwidth,_skipMC=opt.skipMC)
      plotSplines(fm,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.ext),_ext='_ea',_nominalMass=MHNominal,splinesToPlot=['ea'])
      plotSplines(fm,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.ext),_ext='_xs',_nominalMass=MHNominal,splinesToPlot=['xs_interference'])
  else:
    if opt.skipVertexScenarioSplit:
      plotPdfComponents(ssfRV,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.ext),_extension="total_",_proc=procRVFit,_cat=catRVFit, _mass=float(MHNominal))
    if not opt.skipVertexScenarioSplit:
      plotPdfComponents(ssfRV,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.ext),_extension="RV_",_proc=procRVFit,_cat=catRVFit, _mass=float(MHNominal))
      plotPdfComponents(ssfWV,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.ext),_extension="WV_",_proc=procWVFit,_cat=catRVFit, _mass=float(MHNominal))
    # Plot interpolation
    plotInterpolation(fm,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.ext))
    plotSplines(fm,_outdir="%s/outdir_%s/signalFit/Plots"%(swd__,opt.ext),_nominalMass=MHNominal,splinesToPlot=['xs','br','ea'])
