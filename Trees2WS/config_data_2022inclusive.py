
trees2wsCfg = {

  # Name of RooDirectory storing input tree
  'inputTreeDir': 'DiphotonTree',

  # Variables to be added to dataframe (MC): not used here
  'mainVars': ["mass", "weight"],

  # Variables for data RooDataSets
  'dataVars': ["mass", "weight"],

  # STXS var: not used
  'stxsVar': '',

  # Systematic vars: not needed right now
  'systematicsVars': ["mass", "weight"],

  # No theory weights
  'theoryWeightContainers': {},

  # No systematics
  'systematics': [],

  # Explicit category name
  'cats': ["2022inclusive"]
}