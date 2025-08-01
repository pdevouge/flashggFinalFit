
backgroundScriptCfg = {

  # Setup
  'inputWS': '../../2022preEE/workspaces/data/allData_PreEE.root',  # Input RooWorkspace

  # 'cats': 'auto',         # Let it automatically detect category names from the workspace
  'cats': 'inclusive',
  'catOffset': 0,         # No offset needed since this is a single file
  'ext': 'highmass_500-1000_new', # Will name output directories like: outputs_fTest_2022inclusive/
  'year': '2022preEE',         # Shown in plots; adjust if merging multiple years

  # Job submission options
  'batch': 'local',       # You can change to 'condor' if you'd prefer distributed batch mode
  'queue': 'espresso'     # Ignored when batch is 'local'

}
