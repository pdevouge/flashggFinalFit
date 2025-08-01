# Config file: options for signal fitting

_year = '2022preEE'

signalScriptCfg = {

  # Setup
  'inputWSDir':'../../low_width_RSG_all_masses/root/ws_rsg/',
  'procs':'rsg', # if auto: inferred automatically from filenames
  'cats':'inclusive', # if auto: inferred automatically from (0) workspace
  'ext':'highmass_%s_500-1000'%_year,
  'analysis':'highMassAnalysis', # To specify which replacement dataset mapping (defined in ./python/replacementMap.py)
  'year':'%s'%_year, # Use 'combined' if merging all years: not recommended
  'width':'001',
  'massPoints':'500,550,600,650,700,750,800,900,1000', #'130,150,200,250,300,350,400,450,500,550,600,650,700,750,800,900,1000,1250,1500,1750,2000,2250,2500,3000',

  #Photon shape systematics
  'scales':'', # separate nuisance per year
  'scalesCorr':'', # correlated across years
  'scalesGlobal':'', # affect all processes equally, correlated across years
  'smears':'', # separate nuisance per year

  # Job submission options
  'batch':'local', # ['condor','SGE','IC','local']
  'queue':'espresso',

}
