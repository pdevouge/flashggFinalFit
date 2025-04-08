# Script to convert data trees to RooWorkspace (compatible for finalFits)
# Assumes tree names of the format: 
# * Data_<sqrts>_category

import os, sys
import re
from optparse import OptionParser

def get_options():
    parser = OptionParser()
    parser.add_option('--inputConfig', dest='inputConfig', default="", help='Input config: specify list of variables/analysis categories')
    parser.add_option('--inputTreeFile', dest='inputTreeFile', default=None, help='Input tree file')
    parser.add_option('--outputWSDir', dest='outputWSDir', default=None, help='Output dir (default is same as input dir)')
    parser.add_option('--applyMassCut', dest='applyMassCut', default=False, action="store_true", help='Apply cut on mass')
    parser.add_option('--massCutRange', dest='massCutRange', default='100,180', help='mass cut range, e.g. 130,300')
    return parser.parse_args()
(opt, args) = get_options()

from collections import OrderedDict as od
from importlib import import_module

import ROOT
import pandas
import numpy as np
import uproot

# Set global constants
sqrts__ = "13p6TeV"

# Extract mass cut range values
massLow, massHigh = map(float, opt.massCutRange.split(","))
massMid = 0.5 * (massLow + massHigh)

# Name of RooWorkspace directory/object (using mass range in name for clarity)
inputWSName__ = f"tagsDumper/xgg_highmass_{sqrts__}_m{int(massLow)}-{int(massHigh)}"

print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ XGG TREES 2 WS (DATA) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
def leave():
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ XGG TREES 2 WS (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    exit(0)

# Function to add vars to workspace
def add_vars_to_workspace(_ws=None, _dataVars=None):
    # Add intLumi var
    intLumi = ROOT.RooRealVar("intLumi", "intLumi", 1000., 0., 999999999.)
    intLumi.setConstant(True)
    getattr(_ws, 'import')(intLumi)

    _vars = od()
    for var in _dataVars:
        if var == "mass":
            _vars[var] = ROOT.RooRealVar(var, var, massMid, massLow, massHigh)
            _vars[var].setBins(10)
        elif var == "weight":
            _vars[var] = ROOT.RooRealVar(var, var, 0.)
        else:
            _vars[var] = ROOT.RooRealVar(var, var, 1., -999999, 999999)
            _vars[var].setBins(1)
        getattr(_ws, 'import')(_vars[var], ROOT.RooFit.Silence())
    return _vars.keys()

# Function to make RooArgSet
def make_argset(_ws=None, _varNames=None):
    _aset = ROOT.RooArgSet()
    for v in _varNames: _aset.add(_ws.var(v))
    return _aset

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract options from config file:
options = od()
if opt.inputConfig != '':
    if os.path.exists(opt.inputConfig):
    
        # Import config options
        _cfg = import_module(re.sub(".py", "", opt.inputConfig)).trees2wsCfg
        
        #Extract options
        inputTreeDir = _cfg['inputTreeDir']
        dataVars = _cfg['dataVars']
        cats = _cfg['cats']
    else:
        print(f"[ERROR] Config file {opt.inputConfig} does not exist.")
        leave()
else:
    print("[ERROR] Please specify --inputConfig. Leaving.")
    leave()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get tree names using uproot
f = uproot.open(opt.inputTreeFile)
if inputTreeDir == '': listOfTreeNames == f.keys()
else: listOfTreeNames = f[inputTreeDir].keys()

# If auto, determine category names
if cats == 'auto':
    cats = []
    for tn in listOfTreeNames:
        if "sigma" in tn: continue
        c = tn.split("_%s_" % sqrts__)[-1].split(";")[0]
        cats.append(c)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Open input ROOT file
f = ROOT.TFile(opt.inputTreeFile)

# Output workspace directory and file
if opt.outputWSDir is not None: outputWSDir = opt.outputWSDir+"/ws"
else: outputWSDir = "/".join(opt.inputTreeFile.split("/")[:-1])+"/ws"
if not os.path.exists(outputWSDir): os.system("mkdir %s"%outputWSDir)
outputWSFile = outputWSDir+"/"+opt.inputTreeFile.split("/")[-1]
print(" --> Creating output workspace: (%s)"%outputWSFile)

# Create output file and workspace
fout = ROOT.TFile(outputWSFile, "RECREATE")
foutdir = fout.mkdir(inputWSName__.split("/")[0])
foutdir.cd()
ws = ROOT.RooWorkspace(inputWSName__.split("/")[1], inputWSName__.split("/")[1])

# Add vars to workspace and build argset
varNames = add_vars_to_workspace(ws, dataVars)
aset = make_argset(ws, varNames)

# Loop over categories and fill datasets
for cat in cats:
  print(" --> Extracting events from category: %s"%cat)
  if inputTreeDir == '': treeName = "Data_%s_%s"%(sqrts__,cat)
  else: treeName = "%s/Data_%s_%s"%(inputTreeDir,sqrts__,cat)
  print("    * tree: %s"%treeName)
  t = f.Get(treeName)

  # Define dataset for cat
  dname = "Data_%s_%s"%(sqrts__,cat)
  d = ROOT.RooDataSet(dname,dname,aset,'weight')

  nEv = 0
  # Loop over events in tree and add to dataset with weight 1
  for ev in t:
    if opt.applyMassCut:
      if(getattr(ev,"mass") < float(opt.massCutRange.split(",")[0])) | (getattr(ev,"mass") > float(opt.massCutRange.split(",")[1])): continue
    for var in dataVars: 
      if var == "weight": continue
      ws.var(var).setVal(getattr(ev,var))
    d.add(aset,1.)
    nEv += 1
    
  print(f"    -> Added {nEv} events to dataset '{dname}'")
  
  # Add dataset to worksapce
  getattr(ws, 'import')(d)

# Save the workspace
ws.Write()
fout.Close()

print(" --> DONE. Workspace written successfully.")
