import ROOT
import os
import sys
import json
import re
import numpy as np
import pandas as pd
import pickle
from collections import OrderedDict as od
from commonObjects import *
from commonTools import *

class InterferenceModel:
  # Constructor
  def __init__(self,_proc,_cat,_ext,_year,_sqrts,_xvar,_MH,_massPoints,_width):
    self.proc = _proc
    self.cat = _cat
    self.ext = _ext
    self.year = _year
    self.sqrts = _sqrts
    self.name = "%s_%s_%s_%s"%(self.proc,self.year,self.cat,self.sqrts)
    self.xvar = _xvar
    self.MH = _MH
    self.massPoints = _massPoints
    self.width = _width
    self.intLumi = ROOT.RooRealVar("IntLumi","IntLumi",1.,0.,999999999.) # in pb^-1
    self.dPhase = ROOT.RooRealVar("delta","delta",1,0.,np.pi)
    self.Pdfs = od()
    self.Functions = od()
    self.Splines = od()
    self.retrieveSigInfos()
    self.buildInterference()
    self.MH.setVal(750)
    self.MH.setConstant(True)
    self.xvar.setVal(750)
    self.xvar.setConstant(True)


  def retrieveSigInfos(self):
    fsigName = "%s/outdir_%s/signalFit/output/CMS-HGG_sigfit_%s_%s_%s_%s.root"%(swd__,self.ext,self.ext,self.proc,self.year,self.cat)
    fin = ROOT.TFile(fsigName)
    wsin = fin.Get("%s_%s"%(outputWSName__,sqrts__))
    self.Splines['ea'] = wsin.function("fea_%s"%(self.name)).Clone()
    self.Pdfs[self.name] = wsin.pdf("%s_%s"%(outputWSObjectTitle__,self.name)).Clone()


  def make_interference_real(self):
    # (mgg**2-Mx**2)/sqrt((mgg**2-Mx**2)**2 + Mx**2*Gx**2)
    dependents = ROOT.RooArgList()
    Gx = f"sqrt(2) * {self.width}^2 * MH " if self.proc=='rsg' else f"{self.width}*MH"
    formula = f"(CMS_hgg_mass^2-MH^2) / sqrt((CMS_hgg_mass^2-MH^2)^2 + MH^2*({Gx})^2)"

    dependents.add(self.xvar)
    dependents.add(self.MH)

    self.Functions['I_re'] = ROOT.RooFormulaVar("interference_re","",formula, dependents)


  def make_interference_imaginary(self):
    # Mx*Gx/sqrt((mgg**2-Mx**2)**2 + Mx**2*Gx**2)
    dependents = ROOT.RooArgList()
    Gx = f"sqrt(2) * {self.width}^2 * MH " if self.proc=='rsg' else f"{self.width}*MH"
    formula = f"MH*{Gx} / sqrt((CMS_hgg_mass^2-MH^2)^2 + MH^2*({Gx})^2)"

    dependents.add(self.xvar)
    dependents.add(self.MH)

    self.Functions['I_im'] = ROOT.RooFormulaVar("interference_im","",formula, dependents)

  def make_Mb(self):
    script_dir = os.path.abspath( os.path.dirname( __file__ ) )
    ggbox_xs = pd.read_csv('%s/csv/mcfm_xsec_ggbox.csv'%(script_dir)).set_index('mh').eval('xsec/width')

    self.Splines['ggbox_xsec'] = ROOT.RooSpline1D("ggbox_xs_%s"%(self.name),"ggbox_xs_%s"%(self.name), self.xvar, len(ggbox_xs), ggbox_xs.index.to_numpy(), ggbox_xs.to_numpy())

    self.Functions['Mbkg'] = ROOT.RooFormulaVar("Mbkg","Mbkg","@0*@1",ROOT.RooArgList(self.Splines['ea'],self.Splines['ggbox_xsec']))

  def make_Ms(self):
    self.Functions['Msig'] = ROOT.RooFormulaVar("Msig","Msig","@0*@1",ROOT.RooArgList(self.Splines['ea'],self.Pdfs[self.name]))

  def buildInterference(self):
    self.make_interference_imaginary()
    self.make_interference_real()
    self.make_Mb()
    self.make_Ms()

    dependents = ROOT.RooArgList()
    dependents.add(self.Functions['Mbkg'])
    dependents.add(self.Functions['Msig'])
    dependents.add(self.dPhase)
    dependents.add(self.Functions['I_re'])
    dependents.add(self.Functions['I_im'])

    formula = "2*sqrt(@0*@1)*(@3*cos(@2)-@4*sin(@2))"

    self.Pdfs['interference'] = ROOT.RooFormulaVar("interference","",formula, dependents)


