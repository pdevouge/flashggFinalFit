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
  def __init__(self,_proc,_cat,_ext,_year,_sqrts,_xvar,_massPoints,_width):
    self.proc = _proc
    self.cat = _cat
    self.ext = _ext
    self.year = _year
    self.sqrts = _sqrts
    self.name = "%s_%s_%s_%s"%(self.proc,self.year,self.cat,self.sqrts)
    self.xvar = _xvar
    self.massPoints = _massPoints
    self.width = _width
    self.Vars = od()
    self.Vars['dPhi'] = ROOT.RooRealVar("delta","delta",0.,0.,2*np.pi)
    self.Pdfs = od()
    self.Functions = od()
    self.Splines = od()
    self.retrieveSigInfos()
    self.buildInterference()
    self.xvar.setVal(750)
    self.xvar.setConstant(True)


  def retrieveSigInfos(self):
    fsigName = "%s/outdir_%s/signalFit/output/CMS-HGG_sigfit_%s_%s_%s_%s.root"%(swd__,self.ext,self.ext,self.proc,self.year,self.cat)
    fin = ROOT.TFile(fsigName)
    wsin = fin.Get("%s_%s"%(outputWSName__,sqrts__))
    wsin.var('MH').setVal(750)
    self.Vars['MH'] = wsin.var('MH')
    self.Splines['ea'] = wsin.function("fea_%s"%(self.name))
    self.Pdfs['reso_dcb_%s'%self.name] = wsin.pdf('reso_dcb_%s'%self.name)
    self.Pdfs['sig_pdf'] = wsin.pdf("%s_%s"%(outputWSObjectTitle__,self.name))
    self.Functions['sig_norm'] = wsin.function("%s_%s_norm"%(outputWSObjectTitle__,self.name))


  def make_interference_real(self):
    # (mgg**2-Mx**2)/sqrt((mgg**2-Mx**2)**2 + Mx**2*Gx**2)
    dependents = ROOT.RooArgList()
    Gx = f"sqrt(2) * {self.width}^2 * MH " if self.proc=='rsg' else f"{self.width}*MH"
    formula = f"(CMS_hgg_mass^2-MH^2) / sqrt((CMS_hgg_mass^2-MH^2)^2 + MH^2*({Gx})^2)"

    dependents.add(self.xvar)
    dependents.add(self.Vars['MH'])

    self.Functions['I_re'] = ROOT.RooFormulaVar("interference_re","",formula, dependents)


  def make_interference_imaginary(self):
    # Mx*Gx/sqrt((mgg**2-Mx**2)**2 + Mx**2*Gx**2)
    dependents = ROOT.RooArgList()
    Gx = f"sqrt(2) * {self.width}^2 * MH " if self.proc=='rsg' else f"{self.width}*MH"
    formula = f"MH*{Gx} / sqrt((CMS_hgg_mass^2-MH^2)^2 + MH^2*({Gx})^2)"

    dependents.add(self.xvar)
    dependents.add(self.Vars['MH'])

    self.Functions['I_im'] = ROOT.RooFormulaVar("interference_im","",formula, dependents)

  def make_Mb(self):
    script_dir = os.path.abspath( os.path.dirname( __file__ ) )
    ggbox_xs = pd.read_csv('%s/csv/mcfm_xsec_ggbox.csv'%(script_dir)).set_index('mh').eval('xsec/width')
    self.Splines['ggbox_xsec'] = ROOT.RooSpline1D("ggbox_xs_%s"%(self.name),"ggbox_xs_%s"%(self.name), self.xvar, len(ggbox_xs), ggbox_xs.index.to_numpy(), ggbox_xs.to_numpy())

    self.Pdfs['ggbox_pdf'] = ROOT.RooGenericPdf("Mbkg","Mbkg","@0*@1",ROOT.RooArgList(self.Splines['ea'],self.Splines['ggbox_xsec']))

    self.Functions['Mbkg'] = ROOT.RooFFTConvPdf("Mbkg_%s"%self.name, "Mbkg_%s"%self.name, self.xvar, self.Pdfs['ggbox_pdf'], self.Pdfs['reso_dcb_%s'%self.name])

  def make_Ms(self):
    self.Functions['Msig'] = ROOT.RooProduct("Msig_%s"%self.name,"Msig_%s"%self.name,ROOT.RooArgList(self.Functions['sig_norm'],self.Pdfs['sig_pdf']))

  def buildInterference(self):
    self.make_interference_imaginary()
    self.make_interference_real()
    self.make_Mb()
    self.make_Ms()

    dependents = ROOT.RooArgList()
    dependents.add(self.Functions['Mbkg'])
    dependents.add(self.Functions['Msig'])
    dependents.add(self.Vars['dPhi'])
    dependents.add(self.Functions['I_re'])
    dependents.add(self.Functions['I_im'])

    formula = "2*sqrt(@0*@1)*(@3*cos(@2)-@4*sin(@2))"

    self.Functions['interference'] = ROOT.RooFormulaVar("interference","",formula, dependents)

    self.Pdfs['SBI'] = ROOT.RooGenericPdf("sbi_%s"%self.name,"@0+@1+@2",ROOT.RooArgList(self.Functions['Msig'],self.Functions['Mbkg'],self.Functions['interference']))

  def save(self,wsout):
    wsout.imp = getattr(wsout,"import")
    self.xvar.setBins(10000, "cache")
    wsout.imp(self.xvar, ROOT.RooFit.RecycleConflictNodes())
    for sp in self.Splines.keys():
      wsout.imp(self.Splines[sp],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Pdfs['SBI'],ROOT.RooFit.RecycleConflictNodes())
