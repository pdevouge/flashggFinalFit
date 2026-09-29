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
    self.Integrals = od()

    fsigName = "%s/results/outdir_%s/signalFit/output/CMS-HGG_sigfit_%s_%s_%s_%s.root"%(swd__,self.ext,self.ext,self.proc,self.year,self.cat)
    fin = ROOT.TFile(fsigName)
    wsin = fin.Get("%s_%s"%(outputWSName__,sqrts__))
    self.retrieveSigInfos(wsin)
    self.buildInterference(wsin)


  def retrieveSigInfos(self, wsin):

    self.Vars['MH'] = wsin.var('MH')
    self.Vars['truem'] = wsin.var('true_mass')
    # self.Splines['ea'] = wsin.function("effs_Total")
    # self.Splines['ea'].recursiveRedirectServers(ROOT.RooArgSet(self.xvar))
    self.Pdfs['reso_dcb_%s'%self.name] = wsin.pdf('reso_dcb_%s'%self.name)
    self.Pdfs['reso_dcb_%s'%self.name].recursiveRedirectServers(ROOT.RooArgSet(self.xvar))

  def make_interference_real(self):
    # (mgg**2-Mx**2)/sqrt((mgg**2-Mx**2)**2 + Mx**2*Gx**2)
    dependents = ROOT.RooArgList()
    Gx = f"sqrt(2) * {self.width}^2 * MH " if self.proc=='rsg' else f"{self.width}*MH"
    formula = f"(CMS_hgg_mass^2-MH^2) / sqrt((CMS_hgg_mass^2-MH^2)^2 + MH^2*({Gx})^2)"

    dependents.add(self.xvar)
    dependents.add(self.Vars['MH'])

    self.Functions['I_re'] = ROOT.RooFormulaVar("interference_re_%s"%self.name,"",formula, dependents)


  def make_interference_imaginary(self):
    # Mx*Gx/sqrt((mgg**2-Mx**2)**2 + Mx**2*Gx**2)
    dependents = ROOT.RooArgList()
    Gx = f"sqrt(2) * {self.width}^2 * MH " if self.proc=='rsg' else f"{self.width}*MH"
    formula = f"MH*{Gx} / sqrt((CMS_hgg_mass^2-MH^2)^2 + MH^2*({Gx})^2)"

    dependents.add(self.xvar)
    dependents.add(self.Vars['MH'])

    self.Functions['I_im'] = ROOT.RooFormulaVar("interference_im_%s"%self.name,"",formula, dependents)

  def make_Mb(self):
    script_dir = os.path.abspath( os.path.dirname( __file__ ) )
    # ggbox_xs = pd.read_csv('%s/csv/mcfm_xsec_ggbox.csv'%(script_dir)).set_index('mh').eval('xsec/width')
    # csv is dsigma/dm in fb/GeV; the signal side is pb: divide by 1000.
    ggbox_xs = pd.read_csv('%s/csv/sherpa_xsec_ggbox.csv'%(script_dir), comment='#').set_index('x')['y'].sort_index()/1000.
    self.Splines['ggbox_xsec'] = ROOT.RooSpline1D("ggbox_xs_%s"%(self.name),"ggbox_xs_%s"%(self.name), self.xvar, len(ggbox_xs), ggbox_xs.index.to_numpy(dtype=float), ggbox_xs.to_numpy(dtype=float))
    effCsv = '%s/results/outdir_%s/computeGGBoxEff/output/ggbox_eff_%s_%s_%s_%s.csv'%(iwd__,self.ext,self.ext,self.proc,self.year,self.cat)
    if not os.path.exists(effCsv):
      raise Exception("No ggbox efficiency table at %s -- run computeGGBoxEff.py with first")
    ggbox_eff = pd.read_csv(effCsv).set_index('mNom')['eff'].dropna().sort_index()
    self.Splines['ggbox_eff'] = ROOT.RooSpline1D("ggbox_eff_%s"%(self.name),"ggbox_eff_%s"%(self.name), self.xvar, len(ggbox_eff), ggbox_eff.index.to_numpy(dtype=float), ggbox_eff.to_numpy(dtype=float))

    self.Pdfs['Mbkg_pdf'] = ROOT.RooGenericPdf("Mbkg_pdf_%s"%self.name,"Mbkg_pdf_%s"%self.name,"@0*@1",ROOT.RooArgList(self.Splines['ggbox_eff'],self.Splines['ggbox_xsec']))

    self.Functions['Mbkg_func'] = ROOT.RooFormulaVar("Mbkg_func_%s"%self.name,"Mbkg_func_%s"%self.name,"@0*@1",ROOT.RooArgList(self.Splines['ggbox_eff'],self.Splines['ggbox_xsec']))

    self.Pdfs['Mbkg'] = ROOT.RooFFTConvPdf("Mbkg_%s"%self.name, "Mbkg_%s"%self.name, self.xvar, self.Pdfs['Mbkg_pdf'], self.Pdfs['reso_dcb_%s'%self.name])

    # Normalization: raw integral over the fit range, i.e. the range Combine normalises the pdf on.
    # NB getNormIntegral does NOT return that -- it gave 0 at every mass, so the process vanished.
    mp = self.massPoints.split(',')
    minMass, maxMass = int(mp[0]), int(mp[-1])
    mh = np.arange(minMass, maxMass + 1, dtype=np.float64)
    self.Integrals['Mbkg'] = self.Pdfs['Mbkg_pdf'].createIntegral(ROOT.RooArgSet(self.xvar))
    pdf_y = np.empty(len(mh), dtype=np.float64)
    for i, m in enumerate(mh):
        self.Vars['MH'].setVal(m)
        pdf_y[i] = self.Integrals['Mbkg'].getVal()

    MbkgPdfName = self.Pdfs['Mbkg'].GetName()
    self.Functions['Mbkg_norm'] = ROOT.RooSpline1D("%s_norm" % MbkgPdfName,"%s_norm" % MbkgPdfName,self.Vars['MH'],len(mh),mh,pdf_y)

  def make_Ms(self, wsin):

    # Truth-level (UNsmeared) lineshape x rate
    self.Functions['Msig_func'] = wsin.function("Msig_%s"%self.name)
    self.Functions['Msig_func'].recursiveRedirectServers(ROOT.RooArgSet(self.xvar))

    # Same lineshape but CONVOLVED with the resolution
    self.Pdfs['Msig'] = wsin.pdf("%s_%s"%(outputWSObjectTitle__,self.name))
    self.Pdfs['Msig'].recursiveRedirectServers(ROOT.RooArgSet(self.xvar))
    # Its yield, equal to the integral of Msig_func over the fit range.
    self.Functions['Msig_norm'] = wsin.function("%s_%s_norm"%(outputWSObjectTitle__,self.name))
    self.Functions['Msig_norm'].recursiveRedirectServers(ROOT.RooArgSet(self.xvar))
    # (1 + sum rateConst_i*nuisance_i) from finalModel: kept live below so the SBI yield tracks it.
    self.Functions['rate'] = wsin.function("rate_%s"%self.name)


  def buildInterference(self, wsin):
    self.make_interference_imaginary()
    self.make_interference_real()
    self.make_Mb()
    self.make_Ms(wsin)

    dependents = ROOT.RooArgList()
    dependents.add(self.Functions['Mbkg_func'])
    dependents.add(self.Functions['Msig_func'])  # truth-level, NOT self.Pdfs['Msig'] -- avoids double smearing
    dependents.add(self.Vars['dPhi'])
    dependents.add(self.Functions['I_re'])
    dependents.add(self.Functions['I_im'])

    # Full SBI before resolution smearing
    sbi_formula = "@0 + @1 + 2*sqrt(@0*@1)*(@3*cos(@2)-@4*sin(@2))"
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

    # Normalization: integral SBI = N_S + N_B + N_I, term by term.
    #   N_S  = the signal process' OWN _norm node, N_B = the ggbox process' OWN _norm node.
    # Reusing those two nodes makes Combine model cancel exactly by construction.
    #   N_I  = 2*( cos(delta)*Cre - sin(delta)*Cim ),   Cre/Cim = integral sqrt(B*S)*I_re,I_im
    GK = "RooAdaptiveGaussKronrodIntegrator1D"
    B = self.Functions['Mbkg_func']     # ggbox  |A_B|^2
    S = self.Functions['Msig_func']     # signal |A_S|^2

    # integrand sqrt(B*S)*I_re. 
    # NB: Max because B*S is >= 0 analytically, but the ggbox splines can undershoot to a tiny negative value between knots.
    self.Functions['Cre_integrand'] = ROOT.RooFormulaVar("sbi_cre_%s"%self.name,"",
        "sqrt(TMath::Max(@0*@1,0.))*@2", ROOT.RooArgList(B, S, self.Functions['I_re']))
    # Use improved integration method (see finalModel.py in Signal)
    self.Functions['Cre_integrand'].specialIntegratorConfig(ROOT.kTRUE).method1D().setLabel(GK)
    self.Integrals['Cre'] = self.Functions['Cre_integrand'].createIntegral(ROOT.RooArgSet(self.xvar))

    # integrand sqrt(B*S)*I_im. Same as before
    self.Functions['Cim_integrand'] = ROOT.RooFormulaVar("sbi_cim_%s"%self.name,"",
        "sqrt(TMath::Max(@0*@1,0.))*@2", ROOT.RooArgList(B, S, self.Functions['I_im']))
    self.Functions['Cim_integrand'].specialIntegratorConfig(ROOT.kTRUE).method1D().setLabel(GK)
    self.Integrals['Cim'] = self.Functions['Cim_integrand'].createIntegral(ROOT.RooArgSet(self.xvar))

    # Splines over MH
    mp = self.massPoints.split(',')
    mh = np.arange(int(mp[0]), int(mp[-1]) + 1, dtype=np.float64)
    cre = np.empty(len(mh), dtype=np.float64)
    cim = np.empty(len(mh), dtype=np.float64)
    for i, m in enumerate(mh):
        self.Vars['MH'].setVal(m)
        cre[i] = self.Integrals['Cre'].getVal()
        cim[i] = self.Integrals['Cim'].getVal()

    nameRe = "sbi_int_re_%s"%self.name
    nameIm = "sbi_int_im_%s"%self.name
    self.Splines['Cre'] = ROOT.RooSpline1D(nameRe, nameRe, self.Vars['MH'], len(mh), mh, cre)
    self.Splines['Cim'] = ROOT.RooSpline1D(nameIm, nameIm, self.Vars['MH'], len(mh), mh, cim)

    sbiPdfName = self.Pdfs['SBI'].GetName()
    self.Functions['SBI_norm'] = ROOT.RooFormulaVar("%s_norm"%sbiPdfName, "%s_norm"%sbiPdfName, "@0 + @1 + 2*(@2*cos(@3) - @4*sin(@3))",
                                                    ROOT.RooArgList(self.Functions['Msig_norm'], self.Functions['Mbkg_norm'],self.Splines['Cre'], self.Vars['dPhi'], self.Splines['Cim']))

  def save(self,wsout):
    wsout.imp = getattr(wsout,"import")
    self.xvar.setBins(10000, "cache")
    wsout.imp(self.xvar, ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Pdfs['SBI'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Pdfs['Mbkg'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Pdfs['Msig'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Functions['SBI_norm'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Functions['Mbkg_norm'],ROOT.RooFit.RecycleConflictNodes())
    wsout.imp(self.Functions['Msig_norm'],ROOT.RooFit.RecycleConflictNodes())
