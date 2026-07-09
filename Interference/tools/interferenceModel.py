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
  def __init__(self,_proc,_cat,_ext,_year,_sqrts,_xvar,_MHLow,_MHHigh,_massPoints,_width):
    self.proc = _proc
    self.cat = _cat
    self.ext = _ext
    self.year = _year
    self.sqrts = _sqrts
    self.name = "%s_%s_%s_%s"%(self.proc,self.year,self.cat,self.sqrts)
    self.xvar = _xvar
    self.MHLow = _MHLow
    self.MHHigh = _MHHigh
    self.massPoints = _massPoints
    self.width = _width
    self.Vars = od()
    self.Vars['dPhi'] = ROOT.RooRealVar("delta","delta",0.,0.,2*np.pi)
    self.Pdfs = od()
    self.Functions = od()
    self.Splines = od()

    fsigName = "%s/outdir_%s/signalFit/output/CMS-HGG_sigfit_%s_%s_%s_%s.root"%(swd__,self.ext,self.ext,self.proc,self.year,self.cat)
    fin = ROOT.TFile(fsigName)
    wsin = fin.Get("%s_%s"%(outputWSName__,sqrts__))
    self.retrieveSigInfos(wsin)
    self.buildInterference(wsin)


  def retrieveSigInfos(self, wsin):

    self.Vars['MH'] = wsin.var('MH')
    self.Vars['truem'] = wsin.var('true_mass')
    self.Splines['ea'] = wsin.function("effs_Total").Clone()
    self.Pdfs['reso_dcb_%s'%self.name] = wsin.pdf('reso_dcb_%s'%self.name).Clone()

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

    self.Pdfs['Mbkg_pdf'] = ROOT.RooGenericPdf("Mbkg_pdf","Mbkg_pdf","@0*@1",ROOT.RooArgList(self.Splines['ea'],self.Splines['ggbox_xsec']))

    self.Functions['Mbkg_func'] = ROOT.RooFormulaVar("Mbkg_func","Mbkg_func","@0*@1",ROOT.RooArgList(self.Splines['ea'],self.Splines['ggbox_xsec']))

    self.Pdfs['Mbkg'] = ROOT.RooFFTConvPdf("Mbkg_%s"%self.name, "Mbkg_%s"%self.name, self.xvar, self.Pdfs['Mbkg_pdf'], self.Pdfs['reso_dcb_%s'%self.name])

    # Normalization
    mp = self.massPoints.split(',')
    minMass, maxMass = int(mp[0]), int(mp[-1])
    mh = np.arange(minMass, maxMass + 1, dtype=np.float64)
    pdf_y = np.empty(len(mh), dtype=np.float64)
    for i, m in enumerate(mh):
        self.Vars['MH'].setVal(m)
        pdf_y[i] = self.Pdfs['Mbkg_pdf'].getNormIntegral(ROOT.RooArgSet(self.xvar)).getVal()

    MbkgPdfName = self.Pdfs['Mbkg'].GetName()
    self.Functions['Mbkg_norm'] = ROOT.RooSpline1D("%s_norm" % MbkgPdfName,"%s_norm" % MbkgPdfName,self.Vars['MH'],len(mh),mh,pdf_y)

  def make_Ms(self, wsin):
    wsin.Print('v')
    self.Functions['Msig_func'] = wsin.function("Msig").Clone()
    # self.Pdfs['Msig'].redirectServers(ROOT.RooArgSet(self.xvar))

  def buildInterference(self, wsin):
    self.make_interference_imaginary()
    self.make_interference_real()
    self.make_Mb()
    self.make_Ms(wsin)

    dependents = ROOT.RooArgList()
    # dependents.add(self.Functions['Mbkg_func'])
    dependents.add(self.Functions['Msig_func'])
    # dependents.add(self.Vars['dPhi'])
    # dependents.add(self.Functions['I_re'])
    # dependents.add(self.Functions['I_im'])

    # Full SBI before resolution smearing
    sbi_formula = "@0" #"@0 + @1 + 2*sqrt(@0*@1)*(@3*cos(@2)-@4*sin(@2))"
    self.Functions['SBI_func'] = ROOT.RooFormulaVar(
        "sbi_func_%s" % self.name, "",
        sbi_formula, dependents
    )
    self.Pdfs['SBI_truth'] = ROOT.RooGenericPdf(
        "sbi_truth_%s" % self.name, "",
        sbi_formula, dependents
    )

    # Single convolution with resolution
    self.Pdfs['SBI'] = ROOT.RooFFTConvPdf(
        "sbi_%s" % self.name, "sbi_%s" % self.name,
        self.xvar,
        self.Pdfs['SBI_truth'],
        self.Pdfs['reso_dcb_%s' % self.name]
    )

    # Normalization
    mp = self.massPoints.split(',')
    # minMass, maxMass = int(mp[0]), int(mp[-1])
    # mh = np.arange(minMass, maxMass + 1, dtype=np.float64)
    mh = np.array([float(m) for m in mp])
    pdf_y = np.empty(len(mh), dtype=np.float64)
    for i, m in enumerate(mh):
        self.Vars['MH'].setVal(m)
        pdf_y[i] = self.Pdfs['SBI_truth'].getNormIntegral(ROOT.RooArgSet(self.xvar)).getVal()

    sbiPdfName = self.Pdfs['SBI'].GetName()
    self.Functions['SBI_norm'] = ROOT.RooSpline1D("%s_norm" % sbiPdfName,"%s_norm" % sbiPdfName,self.Vars['MH'],len(mh),mh,pdf_y)

  def save(self,wsout):
    wsout.imp = getattr(wsout,"import")
    self.xvar.setBins(10000, "cache")
    wsout.imp(self.xvar, ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Pdfs['SBI'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Pdfs['Mbkg'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Functions['SBI_norm'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Functions['Mbkg_norm'],ROOT.RooFit.RecycleConflictNodes())
