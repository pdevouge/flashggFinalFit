# Functions for plotting the interference step.
#  * the ggbox efficiency tables that feed Mbkg   (called by computeGGBoxEff.py)
#  * the interference model itself                (called by computeInterference.py)
import ROOT
import numpy as np
import pandas as pd
import os
from collections import OrderedDict as od
from commonObjects import *

# ROOT objects drawn on a pad are not owned by python; without a reference here they are
# collected before SaveAs and the pads come out blank.
_keep = []

def _graph(_x, _y, _ex=None, _ey=None, _color=ROOT.kBlack, _style=1, _width=2, _marker=None):
  x, y = np.ascontiguousarray(_x, dtype=np.float64), np.ascontiguousarray(_y, dtype=np.float64)
  if _ey is None:
    g = ROOT.TGraph(len(x), x, y)
  else:
    ex = np.zeros(len(x)) if _ex is None else np.ascontiguousarray(_ex, dtype=np.float64)
    g = ROOT.TGraphErrors(len(x), x, y, np.ascontiguousarray(ex, dtype=np.float64),
                          np.ascontiguousarray(_ey, dtype=np.float64))
  g.SetLineColor(_color); g.SetLineWidth(_width); g.SetLineStyle(_style)
  g.SetMarkerColor(_color)
  if _marker is not None: g.SetMarkerStyle(_marker)
  _keep.append(g)
  return g

def _legend(_x1=0.6, _y1=0.7, _x2=0.88, _y2=0.88):
  leg = ROOT.TLegend(_x1, _y1, _x2, _y2)
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.04)
  _keep.append(leg)
  return leg

def _save(_canv, _name, _outdir, _extension):
  ext = "_%s"%_extension if _extension != '' else ''
  for suf in ('png','pdf'):
    _canv.SaveAs("%s/%s%s.%s"%(_outdir,_name,ext,suf))
  print("   --> %s/%s%s.pdf"%(_outdir,_name,ext))

def _loadEffCsv(_csv):
  if not os.path.exists(_csv): raise Exception("No efficiency table at %s"%_csv)
  df = pd.read_csv(_csv, comment='#')
  # computeGGBoxEff leaves eff empty where no MC lands in the window
  return df[df['eff'].notna()].copy()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plotGGBoxEfficiency(_csv,_outdir='./',_extension=''):
  """Efficiency vs mass, with the spline fit interferenceModel actually interpolates."""
  df = _loadEffCsv(_csv)
  canv = ROOT.TCanvas(); canv.SetLeftMargin(0.15); _keep.append(canv)

  g = _graph(df['mNom'], df['eff'], _ey=df['eff_err'], _color=ROOT.kBlack, _marker=20)
  g.Draw("APE")
  g.SetTitle("GGBox efficiency;m_{#gamma#gamma} [GeV];#epsilon")
  leg = _legend(); leg.AddEntry(g, "efficiency", "lep")

  if 'eff_fit' in df.columns:
    gf = _graph(df['mNom'], df['eff_fit'], _color=ROOT.kRed+1)
    gf.Draw("L SAME")
    leg.AddEntry(gf, "spline fit", "l")
  leg.Draw()

  canv.SetGridx(); canv.SetGridy()
  _save(canv, "ggbox_eff", _outdir, _extension)

def plotGGBoxNumDen(_csv,_outdir='./',_extension=''):
  """
  Numerator and denominator as densities, linear on top and log below.
  """
  df = _loadEffCsv(_csv)
  dm = (df['mHigh'] - df['mLow']).to_numpy(dtype=float)
  x = df['mNom'].to_numpy(dtype=float)

  canv = ROOT.TCanvas("c_numden","c_numden",700,800); _keep.append(canv)
  pLin = ROOT.TPad("pLin","",0,0.5,1,1.0); pLin.SetBottomMargin(0.02); pLin.SetLeftMargin(0.15)
  pLog = ROOT.TPad("pLog","",0,0.0,1,0.5); pLog.SetTopMargin(0.02); pLog.SetBottomMargin(0.15)
  pLog.SetLeftMargin(0.15); pLog.SetLogy()
  pLin.Draw(); pLog.Draw(); _keep.extend([pLin,pLog])

  for pad, isLog in ((pLin,False),(pLog,True)):
    pad.cd()
    gn = _graph(x, df['num']/dm, _ey=df['num_err']/dm, _color=ROOT.kAzure+1, _marker=20)
    gd = _graph(x, df['den']/dm, _ey=df['den_err']/dm, _color=ROOT.kOrange+7, _marker=21)
    gn.Draw("APE")
    gn.SetTitle(";m_{#gamma#gamma} [GeV];d#sigma/dm (a.u.)" if isLog else ";;d#sigma/dm (a.u.)")
    if isLog:
      pos = np.concatenate([(df['num']/dm).values, (df['den']/dm).values])
      pos = pos[pos > 0]
      if len(pos): gn.SetMinimum(0.5*pos.min()); gn.SetMaximum(5*pos.max())
    gd.Draw("PE SAME")
    if not isLog:
      leg = _legend(); leg.AddEntry(gn,"numerator (reco, in category)","lep")
      leg.AddEntry(gd,"denominator (gen)","lep"); leg.Draw()
    pad.SetGridx(); pad.SetGridy()

  _save(canv, "ggbox_num_den", _outdir, _extension)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plotInterferenceModel(ifm,_mass,_lumi,_outdir='./',_extension='',_delta=0.,_nPoints=1200):
  """
  Three panels at one mass: the templates, the interference alone, and I_re/I_im.
  Each density is what Combine sees for that process: <pdf>_norm * pdf(x) * lumi.
  """
  ifm.Vars['MH'].setVal(float(_mass))
  ifm.Vars['dPhi'].setVal(float(_delta))
  lo, hi = float(ifm.MHLow), float(ifm.MHHigh)
  xs = np.linspace(lo, hi, _nPoints)
  aset = ROOT.RooArgSet(ifm.xvar)

  def dens(pdf, norm):
    n = norm.getVal(); out = np.empty(len(xs))
    for i, x in enumerate(xs):
      ifm.xvar.setVal(float(x)); out[i] = pdf.getVal(aset)*n*_lumi
    return out

  S = dens(ifm.Pdfs['Msig'], ifm.Functions['Msig_norm'])
  B = dens(ifm.Pdfs['Mbkg'], ifm.Functions['Mbkg_norm'])
  SBI = dens(ifm.Pdfs['SBI'], ifm.Functions['SBI_norm'])
  I = SBI - S - B

  ire, iim = np.empty(len(xs)), np.empty(len(xs))
  for i, x in enumerate(xs):
    ifm.xvar.setVal(float(x))
    ire[i] = ifm.Functions['I_re'].getVal(); iim[i] = ifm.Functions['I_im'].getVal()

  canv = ROOT.TCanvas("c_model","c_model",900,1000); _keep.append(canv)
  p1 = ROOT.TPad("p1","",0,0.55,1,1.00); p1.SetBottomMargin(0.02); p1.SetLogy()
  p2 = ROOT.TPad("p2","",0,0.30,1,0.55); p2.SetTopMargin(0.02); p2.SetBottomMargin(0.02)
  p3 = ROOT.TPad("p3","",0,0.00,1,0.30); p3.SetTopMargin(0.02); p3.SetBottomMargin(0.28)
  for p in (p1,p2,p3): p.SetLeftMargin(0.13); p.Draw()
  _keep.extend([p1,p2,p3])

  p1.cd()
  gS = _graph(xs,S,_color=ROOT.kGreen+2); gB = _graph(xs,B,_color=ROOT.kRed+1)
  gSBI = _graph(xs,SBI,_color=ROOT.kBlue+1,_style=2)
  pos = np.concatenate([S[S>0],B[B>0],SBI[SBI>0]])
  gS.Draw("AL")
  gS.SetTitle("M_{X} = %g GeV, #delta = %.3f;;Events / GeV"%(_mass,_delta))
  gS.GetXaxis().SetLimits(lo,hi)
  gS.SetMinimum(max(pos.min(), pos.max()*1e-7)); gS.SetMaximum(pos.max()*10)
  gB.Draw("L SAME"); gSBI.Draw("L SAME")
  leg = _legend(0.62,0.62,0.89,0.88)
  leg.AddEntry(gS,"Signal S","l"); leg.AddEntry(gB,"ggbox B","l")
  leg.AddEntry(gSBI,"S+B+I","l"); leg.Draw()

  p2.cd()
  gI = _graph(xs,I,_color=ROOT.kMagenta+1)
  gI.Draw("AL"); gI.SetTitle(";;Interference [Events / GeV]")
  gI.GetXaxis().SetLimits(lo,hi)
  span = np.abs(I).max()*1.2 or 1.
  gI.SetMinimum(-span); gI.SetMaximum(span)
  gI.GetYaxis().SetTitleSize(0.09); gI.GetYaxis().SetLabelSize(0.075)
  gI.GetYaxis().SetTitleOffset(0.6)
  z = ROOT.TLine(lo,0.,hi,0.); z.SetLineStyle(3); z.Draw(); _keep.append(z)

  p3.cd()
  gRe = _graph(xs,ire,_color=ROOT.kOrange+7); gIm = _graph(xs,iim,_color=ROOT.kAzure+1,_style=2)
  gRe.Draw("AL"); gRe.SetTitle(";m_{#gamma#gamma} [GeV];I_{re}, I_{im}")
  gRe.GetXaxis().SetLimits(lo,hi); gRe.SetMinimum(-1.15); gRe.SetMaximum(1.15)
  gRe.GetXaxis().SetTitleSize(0.10); gRe.GetXaxis().SetLabelSize(0.085)
  gRe.GetYaxis().SetTitleSize(0.09); gRe.GetYaxis().SetLabelSize(0.075)
  gRe.GetYaxis().SetTitleOffset(0.6); gRe.GetXaxis().SetTitleOffset(1.1)
  gIm.Draw("L SAME")
  z2 = ROOT.TLine(lo,0.,hi,0.); z2.SetLineStyle(3); z2.Draw(); _keep.append(z2)
  leg2 = _legend(0.72,0.72,0.89,0.95)
  leg2.AddEntry(gRe,"I_{re}","l"); leg2.AddEntry(gIm,"I_{im}","l"); leg2.Draw()

  _save(canv, "interference_model_M%g_delta%.3f"%(_mass,_delta), _outdir, _extension)
