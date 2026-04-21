# Class for performing simultaneous signal fit
import os
import ROOT
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import scipy.stats
from collections import OrderedDict as od
from array import array
import ctypes

# Parameter lookup table for initialisation
# So far defined up to MHPolyOrder=2
pLUT = od()
pLUT['Reso_DCB'] = od()
pLUT['Reso_DCB']['dm'] = [0, -0.01, 0.01]
pLUT['Reso_DCB']['sigma'] = [0.01, 0.001, 0.1]
# pLUT['Reso_DCB']['n1'] = [4., 3., 10.] # nL
# pLUT['Reso_DCB']['n2'] = [10., 4., 20.] # nR
# pLUT['Reso_DCB']['a1'] = [1.5, 0.5, 2.0] # alphaL
# pLUT['Reso_DCB']['a2'] = [1.5, 0.5, 2.0] # alphaR
pLUT['Reso_DCB']['n1'] = [4., 3., 10.] # nL
pLUT['Reso_DCB']['n2'] = [10., 4., 20.] # nR
pLUT['Reso_DCB']['a1'] = [1., 0.8, 2.0] # alphaL
pLUT['Reso_DCB']['a2'] = [1.5, 1.2, 1.8] # alphaR
pLUT['Reso_func'] = od()
pLUT['Reso_func']['dm'] = ['[0]+[1]*x']
pLUT['Reso_func']['sigma'] = ['[0]+[1]*x']
pLUT['Reso_func']['n1'] = ['[0]+[1]*x']
pLUT['Reso_func']['n2'] = ['[0]+[1]*x']
pLUT['Reso_func']['a1'] = ['[0]+[1]*x']
pLUT['Reso_func']['a2'] = ['[0]+[1]*x']
pLUT['DCB'] = od()
pLUT['DCB']['dm_p0'] = [0.1,-2.5,2.5]
pLUT['DCB']['dm_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['dm_p2'] = [0.0,-0.001,0.001]
pLUT['DCB']['sigma_p0'] = [2.,1.,20.]
pLUT['DCB']['sigma_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['sigma_p2'] = [0.0,-0.001,0.001]
pLUT['DCB']['n1_p0'] = [20.,1.00001,500]
pLUT['DCB']['n1_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['n1_p2'] = [0.0,-0.001,0.001]
pLUT['DCB']['n2_p0'] = [20.,1.00001,500]
pLUT['DCB']['n2_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['n2_p2'] = [0.0,-0.001,0.001]
pLUT['DCB']['a1_p0'] = [1.,1.,10.]
pLUT['DCB']['a1_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['a1_p2'] = [0.0,-0.001,0.001]
pLUT['DCB']['a2_p0'] = [1.,1.,20.]
pLUT['DCB']['a2_p1'] = [0.0,-0.1,0.1]
pLUT['DCB']['a2_p2'] = [0.0,-0.001,0.001]
pLUT['Gaussian_wdcb'] = od()
pLUT['Gaussian_wdcb']['dm_p0'] = [0.1,-1.5,1.5]
pLUT['Gaussian_wdcb']['dm_p1'] = [0.01,-0.01,0.01]
pLUT['Gaussian_wdcb']['dm_p2'] = [0.01,-0.01,0.01]
pLUT['Gaussian_wdcb']['sigma_p0'] = [1.5,1.0,4.]
pLUT['Gaussian_wdcb']['sigma_p1'] = [0.0,-0.1,0.1]
pLUT['Gaussian_wdcb']['sigma_p2'] = [0.0,-0.001,0.001]
pLUT['Frac'] = od()
pLUT['Frac']['p0'] = [0.25,0.01,0.99]
pLUT['Frac']['p1'] = [0.,-0.05,0.05]
pLUT['Frac']['p2'] = [0.,-0.0001,0.0001]
pLUT['Gaussian'] = od()
pLUT['Gaussian']['dm_p0'] = [0.1,-5.,5.]
pLUT['Gaussian']['dm_p1'] = [0.0,-0.01,0.01]
pLUT['Gaussian']['dm_p2'] = [0.0,-0.01,0.01]
pLUT['Gaussian']['sigma_p0'] = ['func',0.5,10.0]
pLUT['Gaussian']['sigma_p1'] = [0.0,-0.01,0.01]
pLUT['Gaussian']['sigma_p2'] = [0.0,-0.01,0.01]
pLUT['FracGaussian'] = od()
pLUT['FracGaussian']['p0'] = ['func',0.01,0.99]
pLUT['FracGaussian']['p1'] = [0.01,-0.005,0.005]
pLUT['FracGaussian']['p2'] = [0.00001,-0.00001,0.00001]

# Function to convert sumw2 variance to poisson interval
def poisson_interval(x,eSumW2,level=0.68):
  neff = x**2/(eSumW2**2)
  scale = abs(x)/neff
  l = scipy.stats.gamma.interval(level, neff, scale=scale,)[0]
  u = scipy.stats.gamma.interval(level, neff+1, scale=scale,)[1]
  # protect against no effective entries
  l[neff==0] = 0.
  # protect against no variance
  l[eSumW2==0.] = 0.
  u[eSumW2==0.] = np.inf
  # convert to upper and lower errors
  eLo, eHi = abs(l-x),abs(u-x)
  return eLo, eHi

# Function to calc chi2 for binned fit given pdf, RooDataHist and xvar as inputs
#def calcChi2(x,pdf,d,errorType="Sumw2",_verbose=False,fitRange=[100,180]):
#def calcChi2(x,pdf,d,errorType="Poisson",_verbose=False,fitRange=[110,140]):
def calcChi2(x,pdf,d,errorType="Poisson",_verbose=False,fitRange=[250, 500]):

  k = 0. # number of non empty bins (for calc degrees of freedom)
  normFactor = d.sumEntries()

  # Using numpy and poisson error
  bins, nPdf, nData, eDataSumW2 = [], [],[],[]
  for i in range(d.numEntries()):
    p = d.get(i)
    x.setVal(p.getRealValue(x.GetName()))
    if( x.getVal() < fitRange[0] )|( x.getVal() > fitRange[1] ): continue
    ndata = d.weight()
    if ndata*ndata == 0: continue
    npdf = pdf.getVal(ROOT.RooArgSet(x))*normFactor*d.binVolume()
    eLo, eHi = ctypes.c_double(), ctypes.c_double()
    #eLo, eHi = ROOT.Double(), ROOT.Double()
    d.weightError(eLo,eHi,ROOT.RooAbsData.SumW2)
    bins.append(i)
    nPdf.append(npdf)
    nData.append(ndata)
    eDataSumW2.append(eHi) if npdf>ndata else eDataSumW2.append(eLo)
    k += 1

  # Convert to numpy array
  nPdf = np.asarray(nPdf)
  nData = np.asarray(nData)
  eDataSumW2 = np.asarray([e.value for e in eDataSumW2], dtype=float)
  #eDataSumW2 = np.asarray(eDataSumW2)

  if errorType == 'Poisson':
    # Change error to poisson intervals: take max interval as error
    eLo,eHi = poisson_interval(nData,eDataSumW2,level=0.68)
    #eDataPoisson = 0.5*(eHi+eLo)
    eDataPoisson = np.maximum(eHi,eLo)
    #eDataPoisson = (nPdf>nData)*eHi + (nPdf<=nData)*eLo
    e = eDataPoisson
    # Calculate chi2 terms
    terms = (nPdf-nData)**2/(eDataPoisson**2)
  elif errorType == "Expected":
    # Change error to sqrt pdf entries
    eExpected = np.sqrt(nPdf)
    e = eExpected
    # Calculate chi2 terms
    terms = (nPdf-nData)**2/(eExpected**2)
  else:
    # Use SumW2 terms to calculate chi2
    e = eDataSumW2
    terms = (nPdf-nData)**2/(eDataSumW2**2)

  # If verbose: print to screen
  if _verbose:
    for i in range(len(terms)):
      print(" --> [DEBUG] Bin %g : nPdf = %.6f, nData = %.6f, e(%s) = %.6f --> chi2 term = %.6f"%(bins[i],nPdf[i],nData[i],errorType,e[i],terms[i]))

  # Sum terms
  result = terms.sum()

  return result,k

# Function to add chi2 for multiple mass points
def nChi2Addition(X,ssf,verbose=False):
  # X: vector of param values (updated with minimise function)
  # Loop over parameters and set RooVars
  for i in range(len(X)): ssf.FitParameters[i].setVal(X[i])
  # Loop over datasets: adding chi2 for each mass point
  chi2sum = 0
  K = 0 # number of non empty bins
  C = len(X)-1 # number of fit params (-1 for MH)
  for mp,d in ssf.DataHists['reco_mass'].items():
    ssf.MH.setVal(int(mp))
    chi2, k  = calcChi2(ssf.xvar,ssf.Pdfs['final'],d,_verbose=verbose)
    chi2sum += chi2
    K += k
  # N degrees of freedom
  ndof = K-C
  ssf.setNdof(ndof)
  return chi2sum

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class SimultaneousFit:
  # Constructor
  def __init__(self,_name,_proc,_cat,_effAcc,_datasetForFit,_xvar,_true_mass,_reduced_mass,_MH,_MHLow,_MHHigh,_width,_massPoints,_nBins,_MHPolyOrder,_minimizerMethod,_minimizerTolerance,verbose=True):
    self.name = _name
    self.proc = _proc
    self.cat = _cat
    self.effAcc = _effAcc
    self.datasetForFit = _datasetForFit
    self.xvar = _xvar
    self.true_mass = _true_mass
    self.reduced_mass = _reduced_mass
    self.MH = _MH
    self.MHLow = _MHLow
    self.MHHigh = _MHHigh
    self.width = _width
    self.massPoints = _massPoints
    self.nBins = _nBins
    self.MHPolyOrder = _MHPolyOrder
    self.minimizerMethod = _minimizerMethod
    self.minimizerTolerance = _minimizerTolerance
    self.verbose = verbose
    # Prepare vars
    self.MH.setConstant(False)
    self.MH.setVal(400) # was 125
    self.MH.setBins(int(self.MHHigh) - int(self.MHLow))
    self.dMH = ROOT.RooFormulaVar("dMH","dMH","@0-400.0",ROOT.RooArgList(self.MH))
    # self.xvar.setVal(125)
    self.xvar.setBins(self.nBins)
    # Dicts to store all fit vars, polynomials, pdfs and splines
    self.nGaussians = 1
    self.Vars = od()
    self.Varlists = od()
    self.Polynomials = od()
    self.Pdfs = od()
    self.ResoFuncs = od()
    self.Coeffs = od()
    self.Splines = od()
    # Prepare RooDataHists for fit
    self.DataHists = od()
    self.DataHists['reco_mass'] = od()
    self.DataHists['reduced_mass'] = od()
    self.DataHists['gen_mass'] = od()
    self.prepareDataHists()
    # Fit containers
    self.FitParameters = None
    self.Ndof = None
    self.Chi2 = None
    self.FitResult = None

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function for setting N degrees of freedom
  def setNdof(self,_ndof): self.Ndof = _ndof

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to extract param bounds
  def extractXBounds(self):
    XBounds = []
    for i in range(len(self.FitParameters)): XBounds.append((self.FitParameters[i].getMin(),self.FitParameters[i].getMax()))
    return np.asarray(XBounds)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to extract initial param value vector
  def extractX0(self):
    X0 = []
    for i in range(len(self.FitParameters)): X0.append(self.FitParameters[i].getVal())
    return np.asarray(X0)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to normalise datasets and convert to RooDataHists for calc chi2
  def prepareDataHists(self):
    # Construct low width data
    for k,d in self.datasetForFit['low_w'].items():
      sumw = d.sumEntries()
      drw_scaled   = d.emptyClone()
      self.Vars['scaled_weight'] = ROOT.RooRealVar("weight","weight",-10000,10000)
      for i in range(0,d.numEntries()):
        val = d.get(i).getRealValue(self.reduced_mass.GetName())
        # Artificially enforce the cuts on reduced mass
        # If not done, ROOFit will create under/overflow bins, worsening the fit
        if (val < self.reduced_mass.getMin()) or (val > self.reduced_mass.getMax()):
          continue
        self.reduced_mass.setVal(d.get(i).getRealValue(self.reduced_mass.GetName()))
        self.Vars['scaled_weight'].setVal((1/sumw)*d.weight()*8.1*1000) #TODO: change to correct lumi
        drw_scaled.add(ROOT.RooArgSet(self.xvar,self.true_mass,self.reduced_mass,self.Vars['scaled_weight']),self.Vars['scaled_weight'].getVal())
      # Convert to RooDataHist
      self.DataHists['reduced_mass'][k] = ROOT.RooDataHist("%s_hist_reduced"%d.GetName(),"%s_hist_reduced"%d.GetName(),ROOT.RooArgSet(self.reduced_mass),drw_scaled)

    # Loop over datasets and normalise to 1
    for k,d in self.datasetForFit['nom_w'].items():
      sumw = d.sumEntries()
      drw = d.emptyClone()
      self.Vars['weight'] = ROOT.RooRealVar("weight","weight",-10000,10000)
      for i in range(0,d.numEntries()):
        self.xvar.setVal(d.get(i).getRealValue(self.xvar.GetName()))
        self.true_mass.setVal(d.get(i).getRealValue(self.true_mass.GetName()))
        self.Vars['weight'].setVal((1/sumw)*d.weight())
        drw.add(ROOT.RooArgSet(self.xvar,self.true_mass,self.reduced_mass,self.Vars['weight']),self.Vars['weight'].getVal())
      # Convert to RooDataHist
      self.DataHists['reco_mass'][k] = ROOT.RooDataHist("%s_hist_reco"%d.GetName(),"%s_hist_reco"%d.GetName(),ROOT.RooArgSet(self.xvar),drw)
      self.DataHists['gen_mass'][k] = ROOT.RooDataHist("%s_hist_true"%d.GetName(),"%s_hist_true"%d.GetName(),ROOT.RooArgSet(self.true_mass),drw)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to construct the resolution model from the available mass points
  def fitReducedMasses(self):
    for mass in self.massPoints.split(','):
      k = "res_param_%s"%mass
      self.Varlists[k] = ROOT.RooArgList("%s"%k)

      for f in ['dm','sigma','a1','n1','a2','n2']:
        self.Vars['%s_%s'%(k,f)] = ROOT.RooRealVar("%s_%s"%(k,f),"%s_%s"%(k,f),pLUT['Reso_DCB']["%s"%f][0],pLUT['Reso_DCB']["%s"%f][1],pLUT['Reso_DCB']["%s"%f][2])
        self.Varlists[k].add( self.Vars['%s_%s'%(k,f)] )

      # Fix 'a1' and 'a2'
      self.Vars['%s_n1'%(k)].setVal(5.3) #3.5 #5.3
      self.Vars['%s_n1'%(k)].setConstant(True)

      self.Vars['%s_n2'%(k)].setVal(7) #5 #7
      self.Vars['%s_n2'%(k)].setConstant(True)

      # Build DCB for individual mass
      self.Pdfs['dcb_reso_%s'%mass] = ROOT.RooDoubleCBFast("dcb_reso_%s"%mass,"dcb_reso_%s"%mass,self.reduced_mass,
                                                           self.Vars['%s_dm'%k],
                                                           self.Vars['%s_sigma'%k],
                                                           self.Vars['%s_a1'%k],
                                                           self.Vars['%s_n1'%k],
                                                           self.Vars['%s_a2'%k],
                                                           self.Vars['%s_n2'%k])
      # Fit single DCB to data
      result = self.Pdfs['dcb_reso_%s'%mass].fitTo(self.DataHists['reduced_mass'][mass], ROOT.RooFit.Save(),  ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.SumW2Error(True))
      result.Print()

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def buildResoModel(self,_verbose='Q'):

    self.fitReducedMasses()

    for f in ['dm', 'sigma', 'a1', 'n1', 'a2', 'n2']:
      self.Varlists['reso_func_%s_params'%f] = ROOT.RooArgList('reso_func_%s_params'%f)

      # Create graph for single DCB parameters
      g = ROOT.TGraphErrors(len(self.massPoints.split(',')))
      for i, mass in enumerate(self.massPoints.split(',')):
        k = f"res_param_{mass}"
        var = self.Vars['%s_%s' % (k, f)]
        val = var.getVal()
        err = var.getError()  # Get fit uncertainty
        g.SetPoint(i, float(mass), val)
        g.SetPointError(i, 0, err)
      # Fit function to parameter(mass)
      func = ROOT.TF1("fit_func_%s"%f, pLUT['Reso_func']["%s"%f][0], float(self.MHLow), float(self.MHHigh))
      g.Fit(func, _verbose)
      self.ResoFuncs["%s_function"%f] = func

      # Extract parameters from TF1
      n_params = func.GetNpar()
      formula = pLUT['Reso_func'][f][0].replace('x', 'MH')
      for i in range(n_params):
        pval = func.GetParameter(i)
        self.Vars['reso_func_%s_p%s'%(f,i)] = ROOT.RooRealVar(f"{f}_p{i}", f"{f}_p{i}", pval)
        self.Varlists['reso_func_%s_params'%f].add( self.Vars['reso_func_%s_p%s'%(f,i)] )
        formula = formula.replace('[%d]'%i, '@%d'%i)
      self.Varlists['reso_func_%s_params'%f].add(self.MH)

      # Create formula from function
      roo_formula = ROOT.RooFormulaVar(
          f"{f}_formula", formula, self.Varlists['reso_func_%s_params'%f])

      # Save the RooFormulaVar for use in DCB later
      self.ResoFuncs["%s_formula"%f] = roo_formula
      self.ResoFuncs["%s_formula"%f].Print()

    # Fit new DCB with some of the parameters fixed
    fixed = ['sigma', 'a1', 'a2', 'n1', 'n2']
    for mass in self.massPoints.split(','):
      k = f"res_param_{mass}"
      self.MH.setVal(float(mass))
      self.MH.setConstant(True)

      self.Pdfs['dcb_reso_%s_from_func'%mass] = ROOT.RooDoubleCBFast("dcb_reso_%s"%mass,"dcb_reso_%s"%mass,self.reduced_mass,
                                                          self.ResoFuncs['dm_formula'] if 'dm' in fixed else self.Vars[f'{k}_dm'],
                                                          self.ResoFuncs['sigma_formula'] if 'sigma' in fixed else self.Vars[f'{k}_sigma'],
                                                          self.ResoFuncs['a1_formula'] if 'a1' in fixed else self.Vars[f'{k}_a1'],
                                                          self.ResoFuncs['n1_formula'] if 'n1' in fixed else self.Vars[f'{k}_n1'],
                                                          self.ResoFuncs['a2_formula'] if 'a2' in fixed else self.Vars[f'{k}_a2'],
                                                          self.ResoFuncs['n2_formula'] if 'n2' in fixed else self.Vars[f'{k}_n2'])

      result = self.Pdfs['dcb_reso_%s_from_func'%mass].fitTo(self.DataHists['reduced_mass'][mass], ROOT.RooFit.Save(), ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1), )
      result.Print()

    # Get updated values for unfixed parameters
    for f in ['dm', 'sigma', 'a1', 'n1', 'a2', 'n2']:
      if f in fixed: continue

      g = ROOT.TGraphErrors(len(self.massPoints.split(',')))
      for i, mass in enumerate(self.massPoints.split(',')):
        k = f"res_param_{mass}"
        var = self.Vars['%s_%s' % (k, f)]
        val = var.getVal()
        err = var.getError()  # Get fit uncertainty
        g.SetPoint(i, float(mass), val)
        g.SetPointError(i, 0, err)
      # Fit function to parameter(mass)
      func = ROOT.TF1("fit_func_%s"%f, pLUT['Reso_func']["%s"%f][0], float(self.MHLow), float(self.MHHigh))
      g.Fit(func, _verbose)
      self.ResoFuncs["%s_function"%f] = func

      # Extract parameters from TF1
      n_params = func.GetNpar()
      for i in range(n_params):
        pval = func.GetParameter(i)
        self.Vars['reso_func_%s_p%s'%(f,i)].setVal(pval)
    # NB: no need to store the new formula, since it is already stored in the list of vars, we can just update it

    self.Pdfs['final_dcb_reso_from_func'] = ROOT.RooDoubleCBFast("dcb_reso_model","dcb_reso_model",self.reduced_mass,
                                                        self.ResoFuncs['dm_formula'],
                                                        self.ResoFuncs['sigma_formula'],
                                                        self.ResoFuncs['a1_formula'],
                                                        self.ResoFuncs['n1_formula'],
                                                        self.ResoFuncs['a2_formula'],
                                                        self.ResoFuncs['n2_formula'])

    # Create the resolution function
    self.ResoFuncs['dm_scaled'] = ROOT.RooFormulaVar("dm_scaled", "@0 * @1", ROOT.RooArgList(self.ResoFuncs['dm_formula'], self.MH))
    self.ResoFuncs['sigma_scaled'] = ROOT.RooFormulaVar("sigma_scaled", "@0 * @1", ROOT.RooArgList(self.ResoFuncs['sigma_formula'], self.MH))
    self.Pdfs['final_reso_model'] = ROOT.RooDoubleCBFast("dcb_reso_model","dcb_reso_model",self.xvar,
                                                          self.ResoFuncs['dm_scaled'],
                                                          self.ResoFuncs['sigma_scaled'],
                                                          self.ResoFuncs['a1_formula'],
                                                          self.ResoFuncs['n1_formula'],
                                                          self.ResoFuncs['a2_formula'],
                                                          self.ResoFuncs['n2_formula'])

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def buildSignalSplines(self):
    script_dir = os.path.abspath( os.path.dirname( __file__ ) )
    gtot_o_m3 = pd.read_csv('%s/csv/gtot_o_m3_%s.csv'%(script_dir,self.proc)).set_index('m_x')
    xsec_m = pd.read_csv('%s/csv/xsec_%s.csv'%(script_dir,self.proc)).set_index('m_x')['xsec']
    effs_m = pd.Series(self.effAcc)
    # effs_m = pd.read_csv('%s/csv/effs_%s.csv'%(script_dir,self.proc)).set_index('m')['effs']
    self.Splines['effs'] = ROOT.RooSpline1D("effs_%s"%(self.name),"effs_%s"%(self.name), self.xvar, len(effs_m), effs_m.index.astype(float).to_numpy(), effs_m.to_numpy())

    var_map = {
      'MH': self.MH,
      'm': self.xvar,
      'truem': self.true_mass
    }
    for vname, var in var_map.items():
      for mp in self.massPoints.split(','):
        self.Splines[f'gtot_o_m3_{mp}_{vname}'] = ROOT.RooSpline1D("gtot_o_m3_%s_%s_%s"%(mp,vname,self.name),"gtot_o_m3_%s_%s_%s"%(mp,vname,self.name), var, len(gtot_o_m3), gtot_o_m3.index.to_numpy(), gtot_o_m3[f'{mp}.0'].to_numpy())
      self.Splines[f'xsec_{vname}'] = ROOT.RooSpline1D("xsec_%s_%s"%(vname,self.name),"xsec_%s_%s"%(vname,self.name), var, len(xsec_m), xsec_m.index.to_numpy(), xsec_m.to_numpy())

    if self.proc == 'rsg':
      gtot = pd.read_csv('%s/csv/gtot_rsg.csv'%(script_dir)).set_index('m_x')
      self.Splines[f'gtot_MH'] = ROOT.RooSpline1D("gtot_MH_%s"%(self.name),"gtot_MH_%s"%(self.name), var, len(gtot), gtot.index.to_numpy(), gtot.to_numpy())

  # Construct Pythia model: m^2*Gtot(m,mX)/()^2+m^2*Gtot(m,mX)
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
    dependents.add(self.MH)

    # -- Gf(m)/Gf(mx) --
    # here Gf = Gtot (for the Pythia model, effective branching ratio set to 1)
    terms = []
    for mp in self.massPoints.split(','):
      Gtot_ratio_expr = "%s/%s"%(self.Splines[f'gtot_o_m3_{mp}_m'].GetName(),self.Splines[f'gtot_o_m3_{mp}_MH'].GetName())
      term = f"(abs(MH-{mp})<{eps}) * ({Gtot_ratio_expr})"
      terms.append(term)

    ratio = " + ".join(terms)

    # -- EffxAcc --
    eff = self.Splines['effs'].GetName()
    dependents.add(self.Splines['effs'])

    # -- BW --
    #  - for RSGrav: sigma(m) \propto \kappa()m^2 = kappa0 x m^2 : power_xsec = 2
    #  - power = 2 (from m^2) + 3 (from G_tot) + power_xsec
    if self.proc == 'rsg':
      power = 2 + 3 + 2
    else:
      power = 2 + 3
    formula = f"(CMS_hgg_mass/MH)^{power} / ((CMS_hgg_mass^2 - MH^2)^2 + CMS_hgg_mass^2*({Gtot})^2) * {kf} * ({ratio}) * {eff}"

    self.Pdfs['rel_bw'] = ROOT.RooGenericPdf("rel_bw","",formula, dependents)

    # DEBUG
    # self.xvar.setVal(490)
    # self.MH.setVal(500)
    # print(f'xvar: {self.xvar.getVal()}')
    # print(f'MH: {self.MH.getVal()}')
    # if self.proc == 'rsg': print(f'gxcorr: {self.Splines[f"gtot_MH"].getVal()}')
    # print(f'xsec_test_MH: {self.Splines[f"xsec_MH"].getVal()}')
    # print(f'xsec_test_m: {self.Splines[f"xsec_m"].getVal()}')
    # print(f'gtot_o_m3_m: {self.Splines[f"gtot_o_m3_500_m"].getVal()}')
    # print(f'gtot_o_m3_MH: {self.Splines[f"gtot_o_m3_500_MH"].getVal()}')
    # print(f'effs: {self.Splines["effs"].getVal()}')
    # print(f'rel_bw: {self.Pdfs["rel_bw"].getVal()}')

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def buildAnalytical(self):

    self.xvar.setBins(10000, "cache")
    self.Pdfs['final'] = ROOT.RooFFTConvPdf("final_model", "rel_bw * dcb", self.xvar, self.Pdfs['rel_bw'], self.Pdfs['final_reso_model'])

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def buildDCBplusGaussian(self,_recursive=True):

    # DCB
    # Define polynominal functions (in dMH)
    for f in ['dm','sigma','n1','n2','a1','a2']:
      k = "%s_dcb"%f
      self.Varlists[k] = ROOT.RooArgList("%s_coeffs"%k)
      # Create coeff for polynominal of order MHPolyOrder: y = a+bx+cx^2+...
      for po in range(0,self.MHPolyOrder+1):
        self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar("%s_p%g"%(k,po),"%s_p%g"%(k,po),pLUT['DCB']["%s_p%s"%(f,po)][0],pLUT['DCB']["%s_p%s"%(f,po)][1],pLUT['DCB']["%s_p%s"%(f,po)][2])
        self.Varlists[k].add( self.Vars['%s_p%g'%(k,po)] )
      # Define polynominal
      self.Polynomials[k] = ROOT.RooPolyVar(k,k,self.dMH,self.Varlists[k])
    # Mean function
    self.Polynomials['mean_dcb'] = ROOT.RooFormulaVar("mean_dcb","mean_dcb","(@0+@1)",ROOT.RooArgList(self.MH,self.Polynomials['dm_dcb']))
    # Build DCB
    print(self.Polynomials['mean_dcb'],self.Polynomials['sigma_dcb'],self.Polynomials['a1_dcb'],self.Polynomials['n1_dcb'],self.Polynomials['a2_dcb'],self.Polynomials['n2_dcb'])
    self.Pdfs['dcb'] = ROOT.RooDoubleCBFast("dcb","dcb",self.xvar,self.Polynomials['mean_dcb'],self.Polynomials['sigma_dcb'],self.Polynomials['a1_dcb'],self.Polynomials['n1_dcb'],self.Polynomials['a2_dcb'],self.Polynomials['n2_dcb'])

    # Gaussian
    # Define polynomial function for sigma (in dMH). Gaussian defined to have same mean as DCB
    f = "sigma"
    k = "%s_gaus"%f
    self.Varlists[k] = ROOT.RooArgList("%s_coeffs"%k)
    # Create coeff for polynominal of order MHPolyOrder: y = a+bx+cx^2+...
    for po in range(0,self.MHPolyOrder+1):
      self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar("%s_p%g"%(k,po),"%s_p%g"%(k,po),pLUT['Gaussian_wdcb']["%s_p%s"%(f,po)][0],pLUT['Gaussian_wdcb']["%s_p%s"%(f,po)][1],pLUT['Gaussian_wdcb']["%s_p%s"%(f,po)][2])
      self.Varlists[k].add( self.Vars['%s_p%g'%(k,po)] )
    # Define polynomial
    self.Polynomials[k] = ROOT.RooPolyVar(k,k,self.dMH,self.Varlists[k])
    # Build Gaussian
    self.Pdfs['gaus'] = ROOT.RooGaussian("gaus","gaus",self.xvar,self.Polynomials['mean_dcb'],self.Polynomials['sigma_gaus'])


    # Relative fraction: also polynomial of order MHPolyOrder
    self.Varlists['frac'] = ROOT.RooArgList("frac_coeffs")
    for po in range(0,self.MHPolyOrder+1):
      self.Vars['frac_p%g'%po] = ROOT.RooRealVar("frac_p%g"%po,"frac_p%g"%po,pLUT['Frac']['p%g'%po][0],pLUT['Frac']['p%g'%po][1],pLUT['Frac']['p%g'%po][2])
      self.Varlists['frac'].add( self.Vars['frac_p%g'%po] )
    # Define Polynomial
    self.Polynomials['frac'] = ROOT.RooPolyVar('frac','frac',self.dMH,self.Varlists['frac'])
    # Constrain fraction to not be above 1 or below 0
    self.Polynomials['frac_constrained'] = ROOT.RooFormulaVar("frac_constrained","frac_constrained","(@0>0)*(@0<1)*@0+(@0>1.0)*0.9999",ROOT.RooArgList(self.Polynomials['frac']))
    self.Coeffs['frac_constrained'] = self.Polynomials['frac_constrained' ]

    # Define total PDF
    _pdfs, _coeffs = ROOT.RooArgList(), ROOT.RooArgList()
    for pdf in ['dcb','gaus']: _pdfs.add(self.Pdfs[pdf])
    _coeffs.add(self.Coeffs['frac_constrained'])
    self.Pdfs['final'] = ROOT.RooAddPdf("%s_%s"%(self.proc,self.cat),"%s_%s"%(self.proc,self.cat),_pdfs,_coeffs,_recursive)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def buildNGaussians(self,nGaussians,_recursive=True):

    # Set number of gaussians
    self.nGaussians = nGaussians

    # Loop over NGaussians
    for g in range(0,nGaussians):
      # Define polynominal functions for mean and sigma (in MH)
      for f in ['dm','sigma']:
        k = "%s_g%g"%(f,g)
        self.Varlists[k] = ROOT.RooArgList("%s_coeffs"%k)
        # Create coeff for polynominal of order MHPolyOrder: y = a+bx+cx^2+...
        for po in range(0,self.MHPolyOrder+1):
          # p0 value of sigma is function of g (creates gaussians of increasing width)
          if(f == "sigma")&(po==0):
            self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar("%s_p%g"%(k,po),"%s_p%g"%(k,po),(g+1)*1.0,pLUT['Gaussian']["%s_p%s"%(f,po)][1],pLUT['Gaussian']["%s_p%s"%(f,po)][2])
          else:
            self.Vars['%s_p%g'%(k,po)] = ROOT.RooRealVar("%s_p%g"%(k,po),"%s_p%g"%(k,po),pLUT['Gaussian']["%s_p%s"%(f,po)][0],pLUT['Gaussian']["%s_p%s"%(f,po)][1],pLUT['Gaussian']["%s_p%s"%(f,po)][2])
          self.Varlists[k].add( self.Vars['%s_p%g'%(k,po)] )
        # Define polynominal
        self.Polynomials[k] = ROOT.RooPolyVar(k,k,self.dMH,self.Varlists[k])
      # Mean function
      self.Polynomials['mean_g%g'%g] = ROOT.RooFormulaVar("mean_g%g"%g,"mean_g%g"%g,"(@0+@1)",ROOT.RooArgList(self.MH,self.Polynomials['dm_g%g'%g]))
      # Build Gaussian
      self.Pdfs['gaus_g%g'%g] = ROOT.RooGaussian("gaus_g%g"%g,"gaus_g%g"%g,self.xvar,self.Polynomials['mean_g%g'%g],self.Polynomials['sigma_g%g'%g])

      # Relative fractions: also polynomials of order MHPolyOrder (define up to n=nGaussians-1)
      if g < nGaussians-1:
        self.Varlists['frac_g%g'%g] = ROOT.RooArgList("frac_g%g_coeffs"%g)
        for po in range(0,self.MHPolyOrder+1):
          if po == 0:
            self.Vars['frac_g%g_p%g'%(g,po)] = ROOT.RooRealVar("frac_g%g_p%g"%(g,po),"frac_g%g_p%g"%(g,po),0.5-0.05*g,pLUT['FracGaussian']['p%g'%po][1],pLUT['FracGaussian']['p%g'%po][2])
          else:
            self.Vars['frac_g%g_p%g'%(g,po)] = ROOT.RooRealVar("frac_g%g_p%g"%(g,po),"frac_g%g_p%g"%(g,po),pLUT['FracGaussian']['p%g'%po][0],pLUT['FracGaussian']['p%g'%po][1],pLUT['FracGaussian']['p%g'%po][2])
          self.Varlists['frac_g%g'%g].add( self.Vars['frac_g%g_p%g'%(g,po)] )
        # Define Polynomial
        self.Polynomials['frac_g%g'%g] = ROOT.RooPolyVar("frac_g%g"%g,"frac_g%g"%g,self.dMH,self.Varlists['frac_g%g'%g])
        # Constrain fraction to not be above 1 or below 0
        self.Polynomials['frac_g%g_constrained'%g] = ROOT.RooFormulaVar('frac_g%g_constrained'%g,'frac_g%g_constrained'%g,"(@0>0)*(@0<1)*@0+ (@0>1.0)*0.9999",ROOT.RooArgList(self.Polynomials['frac_g%g'%g]))
        self.Coeffs['frac_g%g_constrained'%g] = self.Polynomials['frac_g%g_constrained'%g]
    # End of loop over n Gaussians

    # Define total PDF
    _pdfs, _coeffs = ROOT.RooArgList(), ROOT.RooArgList()
    for g in range(0,nGaussians):
      _pdfs.add(self.Pdfs['gaus_g%g'%g])
      if g < nGaussians-1: _coeffs.add(self.Coeffs['frac_g%g_constrained'%g])
    self.Pdfs['final'] = ROOT.RooAddPdf("%s_%s"%(self.proc,self.cat),"%s_%s"%(self.proc,self.cat),_pdfs,_coeffs,_recursive)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def runFit(self):
    # Extract fit variables: remove xvar from fit parameters
    fv = self.Pdfs['final'].getVariables().Clone()
    fv.remove(self.xvar)
    self.FitParameters = ROOT.RooArgList(fv)

    # Create initial vector of parameters and calculate initial Chi2
    if self.verbose: print("\n --> (%s) Initialising fit parameters"%self.name)
    x0 = self.extractX0()
    xbounds = self.extractXBounds()
    self.Chi2 = self.getChi2()
    # Print parameter pre-fit values
    if self.verbose: self.printFitParameters(title="Pre-fit")

    # Run fit
    if self.verbose: print(" --> (%s) Running fit"%self.name)
    self.FitResult = minimize(nChi2Addition,x0,args=self,bounds=xbounds)#,method=self.minimizerMethod)
    self.Chi2 = self.getChi2()
    #self.Chi2 = nChi2Addition(self.FitResult['x'],self)
    # Print parameter post-fit values
    if self.verbose: self.printFitParameters(title="Post-fit")

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function build RooSpline1D from Polynomials to model MH dependence on fit params
  def buildSplines(self):
    # Loop over polynomials
    for k, poly in self.Polynomials.items():
      _x, _y = [], []
      _mh = float(self.MHLow)
      while(_mh<float(self.MHHigh) + .1):
        self.MH.setVal(_mh)
        _x.append(_mh)
        _y.append(poly.getVal())
        _mh += 25
      # Convert to arrays
      arr_x, arr_y = array('f',_x), array('f',_y)
      # Create spline and save to dict
      self.Splines[k] = ROOT.RooSpline1D(poly.GetName(),poly.GetName(),self.MH,len(_x),arr_x,arr_y)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to print fit values
  def printFitParameters(self,title="Fit"):
    print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(" --> (%s) %s parameter values:"%(self.name,title))
    # Skip MH
    for i in range(1,len(self.FitParameters)): print("    * %-20s = %.6f"%(self.FitParameters[i].GetName(),self.FitParameters[i].getVal()))
    print("    ~~~~~~~~~~~~~~~~")
    print("    * chi2 = %.6f, n(dof) = %g --> chi2/n(dof) = %.3f"%(self.getChi2(),int(self.Ndof),self.getChi2()/int(self.Ndof)))
    print("    ~~~~~~~~~~~~~~~~")
    print("    * [VERBOSE] chi2 = %.6f"%(self.getChi2(verbose=False)))
    print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to set vars from json file
  def setVars(self,_json):
    with open(_json,"r") as jf: _vals = json.load(jf)
    for k,v in _vals.items(): self.Vars[k].setVal(v)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to re-calculate chi2 after setting vars
  def getChi2(self,verbose=False):
    x = self.extractX0()
    self.Chi2 = nChi2Addition(x,self,verbose=verbose)
    return self.Chi2

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Function to re-calculate chi2/ndof after setting vars
  def getReducedChi2(self):
    x = self.extractX0()
    self.Chi2 = nChi2Addition(x,self)
    return self.Chi2/int(self.Ndof)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

