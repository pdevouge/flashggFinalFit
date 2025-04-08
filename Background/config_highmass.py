
backgroundScriptCfg = {

  # Setup
  'inputWS': 'workspaces/data/data_2022inclusive.root',  # Input RooWorkspace

  # 'cats': 'auto',         # Let it automatically detect category names from the workspace
  'cats': '2022inclusive',
  'catOffset': 0,         # No offset needed since this is a single file
  'ext': 'm130-300', # Will name output directories like: outputs_fTest_2022inclusive/
  'year': '2022',         # Shown in plots; adjust if merging multiple years

  # Job submission options
  'batch': 'local',       # You can change to 'condor' if you'd prefer distributed batch mode
  'queue': 'espresso'     # Ignored when batch is 'local'

}
