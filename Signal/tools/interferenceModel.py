import ROOT
import os
import sys
import json
import re
import numpy as np
import pandas
import pickle
from collections import OrderedDict as od
from commonObjects import *
from commonTools import *
from signalTools import *

class InterferenceModel:
  # Constructor
  def __init__(self,_ssfMap,_proc,_cat,_ext,_year,_sqrts,_datasets,_xvar,_MH,_MHNominal,_MHLow,_MHHigh,_massPoints,_width):
    self.ssfMap = _ssfMap
    self.proc = _proc
    self.cat = _cat
    self.ext = _ext
    self.year = _year
    self.sqrts = _sqrts
    self.name = "%s_%s_%s_%s"%(self.proc,self.year,self.cat,self.sqrts)
    self.datasets = _datasets
    self.xvar = _xvar
    self.aset = ROOT.RooArgSet(self.xvar)
    self.MH = _MH
    self.MHNominal = _MHNominal
    self.MHLow = _MHLow
    self.MHHigh = _MHHigh
    self.massPoints = _massPoints
    self.width = _width
    self.intLumi = ROOT.RooRealVar("IntLumi","IntLumi",1.,0.,999999999.) # in pb^-1