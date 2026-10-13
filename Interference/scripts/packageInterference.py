# Script to package individual interference models for a single category in one ROOT file.
# Merging per-era extensions the same way Signal/scripts/packageSignal.py does.

import os, sys
import glob
import ROOT
from optparse import OptionParser
from commonObjects import *

def get_options():
  parser = OptionParser()
  parser.add_option("--cat", dest='cat', default='', help="RECO category to package")
  parser.add_option("--proc", dest='proc', default='', help="Signal process")
  parser.add_option("--exts", dest='exts', default='', help="Comma separated list of per-era extensions to merge")
  parser.add_option("--outputExt", dest='outputExt', default='packaged', help="Output extension")
  parser.add_option("--mergeYears", dest='mergeYears', default=False, action="store_true", help="Merge specified categories across years")
  parser.add_option("--year", dest="year", default="2016", help="If not merging, then specify year for output file name")
  return parser.parse_args()
(opt,args) = get_options()

def rooiter(x):
  it = x.iterator()
  ret = it.Next()
  while ret:
    yield ret
    ret = it.Next()

# Extract all files to be merged
fNames = {}
for ext in opt.exts.split(","):
  fNames[ext] = glob.glob("results/outdir_%s/computeIntf/output/CMS-HGG_intfm_%s_*_%s.root"%(ext,ext,opt.cat))

print(" --> Packaging output workspace")
packagedWS = ROOT.RooWorkspace("%s"%(intfWSName__),"%s"%(intfWSName__))
packagedWS.imp = getattr(packagedWS,"import")

for ext, fNames_by_ext in fNames.items():
  if not fNames_by_ext:
    print(" --> [ERROR] No interference model file found for ext=%s, cat=%s. "
          "Did the interference step run for that era?"%(ext,opt.cat))
    sys.exit(1)
  for fName in fNames_by_ext:
    fin = ROOT.TFile(fName)
    wsin = fin.Get(intfWSName__)
    if not wsin: continue
    allVars, allFunctions, allPdfs = {}, {}, {}
    for _var in rooiter(wsin.allVars()): allVars[_var.GetName()] = _var
    for _func in rooiter(wsin.allFunctions()): allFunctions[_func.GetName()] = _func
    for _pdf in rooiter(wsin.allPdfs()): allPdfs[_pdf.GetName()] = _pdf

    for _varName, _var in allVars.items(): packagedWS.imp(_var,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
    for _funcName, _func in allFunctions.items(): packagedWS.imp(_func,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
    for _pdfName, _pdf in allPdfs.items(): packagedWS.imp(_pdf,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
    fin.Close()

# Save to file
if not os.path.isdir("results/outdir_%s/computeIntf/output"%opt.outputExt):
  os.system("mkdir -p results/outdir_%s/computeIntf/output"%opt.outputExt)

if opt.mergeYears:
  foutName = "results/outdir_%s/computeIntf/output/CMS-HGG_intfm_%s_%s_%s.root"%(opt.outputExt,opt.outputExt,opt.proc,opt.cat)
else:
  foutName = "results/outdir_%s/computeIntf/output/CMS-HGG_intfm_%s_%s_%s_%s.root"%(opt.outputExt,opt.outputExt,opt.proc,opt.cat,opt.year)

print(" --> Writing to: ./%s"%foutName)
f = ROOT.TFile(foutName,"RECREATE")
packagedWS.Write()
f.Close()
