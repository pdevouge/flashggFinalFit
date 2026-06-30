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


  def retrieveSigInfos(self):
    fsigName = "%s/outdir_%s/signalFit/output/CMS-HGG_sigfit_%s_%s_%s_%s.root"%(swd__,self.ext,self.ext,self.proc,self.year,self.cat)
    fin = ROOT.TFile(fsigName)
    wsin = fin.Get("%s_%s"%(outputWSName__,sqrts__))
    wsin.var('MH').setVal(750)
    self.Vars['MH'] = wsin.var('MH')
    self.Vars['truem'] = wsin.var('true_mass')
    self.Splines['ea'] = wsin.function("effs_Total")
    self.Splines['ea'].redirectServers(ROOT.RooArgSet(self.xvar))
    self.Pdfs['reso_dcb_%s'%self.name] = wsin.pdf('reso_dcb_%s'%self.name)
    self.Splines['ea'].redirectServers(ROOT.RooArgSet(self.xvar))

    # self.Pdfs['sig_pdf'] = wsin.pdf("%s_%s"%(outputWSObjectTitle__,self.name))
    # self.Functions['sig_norm'] = wsin.function("%s_%s_norm"%(outputWSObjectTitle__,self.name))


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

    MbkgPdfName = self.Pdfs['Mbkg'].GetName()
    # Can be replaced by self.Functions['Mbkg_func'].createIntegral
    self.Functions['Mbkg_norm'] = self.Pdfs['Mbkg_pdf'].createIntegral(
        ROOT.RooArgSet(self.xvar)
    )
    self.Functions['Mbkg_norm'].SetName("%s_norm" % MbkgPdfName)

  def make_Ms(self):
    self.ZZ_buildTrueLineshape()

    self.Pdfs['Msig'] = ROOT.RooFFTConvPdf("Msig_%s"%self.name, "Msig_%s"%self.name, self.xvar, self.Pdfs['Msig_pdf'], self.Pdfs['reso_dcb_%s'%self.name])

    MsigPdfName = self.Pdfs['Msig'].GetName()
    # createIntegral returns a live RooAbsReal that recomputes when MH changes
    self.Functions['Msig_norm'] = self.Pdfs['Msig_pdf'].createIntegral(
        ROOT.RooArgSet(self.xvar)
    )
    self.Functions['Msig_norm'].SetName("%s_norm" % MsigPdfName)

  def buildInterference(self):
    self.make_interference_imaginary()
    self.make_interference_real()
    self.make_Mb()
    self.make_Ms()

    dependents = ROOT.RooArgList()
    dependents.add(self.Functions['Mbkg_func'])
    dependents.add(self.Functions['Msig_func'])
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

    sbiPdfName = self.Pdfs['SBI'].GetName()
    # createIntegral returns a live RooAbsReal that recomputes when MH changes
    self.Functions['SBI_norm'] = self.Pdfs['SBI_truth'].createIntegral(
        ROOT.RooArgSet(self.xvar)
    )
    self.Functions['SBI_norm'].SetName("%s_norm" % sbiPdfName)

  def buildSignalSplines(self):
    script_dir = os.path.abspath( os.path.dirname( __file__ ) )
    gtot_o_m3 = pd.read_csv('%s/tools/csv/gtot_o_m3_%s.csv'%(swd__,self.proc)).set_index('m_x')
    xsec_m = pd.read_csv('%s/tools/csv/xsec_%s.csv'%(swd__,self.proc)).set_index('m_x')['xsec']

    var_map = {
      'MH': self.Vars['MH'],
      'm': self.xvar,
      'truem': self.Vars['truem']
    }
    for vname, var in var_map.items():
      for mp in self.massPoints.split(','):
        self.Splines[f'gtot_o_m3_{mp}_{vname}'] = ROOT.RooSpline1D("gtot_o_m3_%s_%s_%s"%(mp,vname,self.name),"gtot_o_m3_%s_%s_%s"%(mp,vname,self.name), var, len(gtot_o_m3), gtot_o_m3.index.to_numpy(), gtot_o_m3[f'{mp}.0'].to_numpy())
      self.Splines[f'xsec_{vname}'] = ROOT.RooSpline1D("xsec_%s_%s"%(vname,self.name),"xsec_%s_%s"%(vname,self.name), var, len(xsec_m), xsec_m.index.to_numpy(), xsec_m.to_numpy())

    if self.proc == 'rsg':
      gtot = pd.read_csv('%s/tools/csv/gtot_rsg.csv'%(swd__)).set_index('m_x')
      self.Splines[f'gtot_MH'] = ROOT.RooSpline1D("gtot_MH_%s"%(self.name),"gtot_MH_%s"%(self.name), var, len(gtot), gtot.index.to_numpy(), gtot.to_numpy())

  def buildTrueLineshape(self):
    dependents = ROOT.RooArgList()
    self.buildSignalSplines()

    # -- xsec --
    dependents.add(self.Splines[f'xsec_MH'])
    dependents.add(self.Splines[f'xsec_m'])
    kf = "(%s / %s)"%(self.Splines[f'xsec_m'].GetName(), self.Splines[f'xsec_MH'].GetName())

    # -- Gtot --
    if self.proc == 'rsg':
      gx_corr = self.Splines[f'gtot_MH'].GetName()
      dependents.add(self.Splines[f'gtot_MH'])
      power_bw = 1
    else:
      gx_corr = 1
      power_bw = 3
    Gx = f"sqrt(2) * {self.width}^2 * MH " if self.proc=='rsg' else f"{self.width}*MH"

    # Here we construct Gtot(m,MX). The df gtot_o_m3 accounts for the mX dependence, by having multiple columns, one for each mX.
    # To avoid the tedious task of constructing one signal model per mass point, we express Gtot as follow:
    # ∑overX:(MH-mX)<eps*(Gtot(m,mX))  -> If MH~=mX, then Gtot si correctly computed
    eps = 1e-6
    terms = []
    for mp in self.massPoints.split(','):
        Gtot_expr = "%s * %s * %s/%s * (CMS_hgg_mass/MH)^%s"%(Gx,gx_corr,self.Splines[f'gtot_o_m3_{mp}_m'].GetName(),self.Splines[f'gtot_o_m3_{mp}_MH'].GetName(),power_bw)
        term = f"(abs(MH-{mp})<{eps}) * ({Gtot_expr})"
        terms.append(term)
        dependents.add(self.Splines[f'gtot_o_m3_{mp}_m'])
        dependents.add(self.Splines[f'gtot_o_m3_{mp}_MH'])

    Gtot = " + ".join(terms)

    dependents.add(self.xvar)
    dependents.add(self.Vars['MH'])

    # -- Gf(m)/Gf(mx) --
    # here Gf = Gtot (for the Pythia model, effective branching ratio set to 1)
    terms = []
    for mp in self.massPoints.split(','):
      Gtot_ratio_expr = "%s/%s"%(self.Splines[f'gtot_o_m3_{mp}_m'].GetName(),self.Splines[f'gtot_o_m3_{mp}_MH'].GetName())
      term = f"(abs(MH-{mp})<{eps}) * ({Gtot_ratio_expr})"
      terms.append(term)

    ratio = " + ".join(terms)

    # -- EffxAcc --
    eff = self.Splines['ea'] .GetName()
    dependents.add(self.Splines['ea'] )

    # -- BW --
    #  - for RSGrav: sigma(m) \propto \kappa()m^2 = kappa0 x m^2 : power_xsec = 2
    #  - power = 2 (from m^2) + 3 (from G_tot) + power_xsec
    if self.proc == 'rsg':
      power = 2 + 3 + 2
    else:
      power = 2 + 3
    formula = f"(CMS_hgg_mass/MH)^{power} / ((CMS_hgg_mass^2 - MH^2)^2 + CMS_hgg_mass^2*({Gtot})^2) * {kf} * ({ratio}) * {eff}"

    self.Pdfs['Msig_pdf'] = ROOT.RooGenericPdf("Msig_pdf","",formula, dependents)
    self.Functions['Msig'] = ROOT.RooFormulaVar("Msig","",formula, dependents)

  def ZZ_buildSignalSplines(self, decay='hgg'):
    script_dir = os.path.abspath( os.path.dirname( __file__ ) )
    mcfm = pd.read_csv('%s/tools/csv/mcfm_xsec_v2.csv'%swd__).set_index('m_x')
    limit_pb = pd.read_csv('%s/tools/csv/limitExp_spin0_138fb_Oct2025.csv'%swd__).set_index('mh')['up_pb']
    gx_lhc_gev = pd.read_csv('%s/tools/csv/lhchwg_hsm_width_v2.csv'%swd__).set_index('mh')
    xsec_lhc_pb = pd.read_csv('%s/tools/csv/lhchwg_hsm_xsec.csv'%swd__).set_index('mh')['xsec']
    xsec_mcfm_lo: pd.Series = mcfm['xsec_lo']
    xsec_mcfm_lo[xsec_mcfm_lo.index <= 200] = xsec_lhc_pb[xsec_lhc_pb.index <= 200]
    # xsec_mcfm_lo = mcfm.query('xsec_lo > 1e-10')['xsec_lo']
    xsec_mcfm_pl = mcfm['xsec']
    xsec_mcfm_py = mcfm['xsec_pythia']
    kf_mcfm = xsec_lhc_pb[200] / xsec_mcfm_lo[200]

    self.Splines['xsec_ul'] = ROOT.RooSpline1D("xsec_ul_%s"%(self.name),"xsec_ul_%s"%(self.name), self.Vars['MH'], len(limit_pb), limit_pb.index.to_numpy(), limit_pb.to_numpy())
    xsec_sm_ = lambda x : np.interp(x, xsec_mcfm_lo.index, xsec_mcfm_lo.values * kf_mcfm)
    xsec_n3lo_ggF = 48.6
    kf_n3lo =  xsec_n3lo_ggF / xsec_sm_(125)
    kf_pythia = xsec_n3lo_ggF / xsec_mcfm_py.loc[120]

    var_map = {
      'MH': self.Vars['MH'],
      'm': self.xvar,
      'truem': self.Vars['truem']
    }
    for vname, var in var_map.items():
      self.Splines[f'xsec_sm_{vname}'] = ROOT.RooSpline1D("xsec_sm_%s_%s"%(vname,self.name),"xsec_sm_%s_%s"%(vname,self.name), var, len(xsec_mcfm_lo), xsec_mcfm_lo.index.to_numpy(), xsec_mcfm_lo.to_numpy() * kf_mcfm) # * kf_n3lo)
      self.Splines[f'xsec_pl_{vname}'] = ROOT.RooSpline1D("xsec_pl_%s_%s"%(vname,self.name),"xsec_pl_%s_%s"%(vname,self.name), var, len(xsec_mcfm_pl), xsec_mcfm_pl.index.to_numpy(), xsec_mcfm_pl.to_numpy())
      self.Splines[f'xsec_py_{vname}'] = ROOT.RooSpline1D("xsec_py_%s_%s"%(vname,self.name),"xsec_py_%s_%s"%(vname,self.name), var, len(xsec_mcfm_py), xsec_mcfm_py.index.to_numpy(), xsec_mcfm_py.to_numpy() * kf_pythia)
      self.Splines[f'ghgg_sm_{vname}'] = ROOT.RooSpline1D("ghgg_sm_%s_%s"%(vname,self.name),"ghgg_sm_%s_%s"%(vname,self.name), var, len(gx_lhc_gev), gx_lhc_gev.index.to_numpy(), gx_lhc_gev[decay].to_numpy())

  def ZZ_buildTrueLineshape(self, decay='hgg', xsec='sm'):
    dependents = ROOT.RooArgList()
    self.ZZ_buildSignalSplines(decay)
    br_x = "(%s / %s)"%(self.Splines['xsec_ul'].GetName(),self.Splines['xsec_sm_MH'].GetName())
    dependents.add(self.Splines['xsec_ul'])
    dependents.add(self.Splines['xsec_sm_MH'])

    if xsec == 'sm':
      xsec_type = 'xsec_sm'
    elif xsec == 'pythia':
      xsec_type = 'xsec_py'
    else:
      xsec_type = 'xsec_pl'

    if xsec != 'sm': dependents.add(self.Splines[f'{xsec_type}_MH'])
    dependents.add(self.Splines[f'{xsec_type}_m'])

    kf = "(%s / %s)"%(self.Splines['xsec_sm_MH'].GetName(), self.Splines[f'{xsec_type}_MH'].GetName())
    xsec = self.Splines[f'{xsec_type}_m'].GetName()

    dependents.add(self.Splines['ghgg_sm_MH'])
    dependents.add(self.Splines['ghgg_sm_m'])

    m_mx = "(%s / %s)"%(self.xvar.GetName(),self.Vars['MH'].GetName())
    dependents.add(self.xvar)
    dependents.add(self.Vars['MH'])

    # -- EffxAcc --
    eff = self.Splines['ea'] .GetName()
    dependents.add(self.Splines['ea'] )

    formula = f"{br_x} * {xsec} * {kf} * ghgg_sm_m_{self.name} / ghgg_sm_MH_{self.name} * {eff}\
                * 2/pi * {self.width} * ({m_mx})^2 / ((({m_mx})^2 - 1)^2 + {self.width}^2) * 1 / MH"

    self.Pdfs['Msig_pdf'] = ROOT.RooGenericPdf("Msig_pdf","",formula, dependents)
    self.Functions['Msig_func'] = ROOT.RooFormulaVar("Msig_func","",formula, dependents)

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
