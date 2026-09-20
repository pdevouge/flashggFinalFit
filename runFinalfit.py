# Script for running the entire pipeline for flashggFinalFit
import os, sys, yaml, glob, shutil, subprocess
import numpy as np
from optparse import OptionParser
from collections import OrderedDict as od

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING FINALFIT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

def get_options():
  parser = OptionParser()
  parser.add_option('--runOnly', dest='run_only', default='', help="Run only given steps (trees, signal, background, datacard, combine, text2ws, limits, collect)")
  parser.add_option('--doSystematics', dest='do_syst', action='store_true', help="Run with systematics")
  parser.add_option('--skipIntf', dest='skip_intf', action='store_true', help="Skip interference making")
  return parser.parse_args()
(opt,args) = get_options()

base_dir = os.getcwd()

def leave():
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING FINALFIT (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  exit(0)

def run(cmd):
  ret = subprocess.call(cmd, shell=True)
  if ret < 0:
    print(f"[ERROR] Command was interrupted by signal {-ret}. Stopping runFinalfit.py.")
    sys.exit(1)
  elif ret != 0:
    print(f"[ERROR] Command failed (exit code {ret}). Stopping runFinalfit.py.")
    sys.exit(ret)
  return ret

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract options from config file
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXTRACTING CONFIG FROM YAML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
if opt.do_syst: print("-----> 'doSystematics': FinalFit will run with systematics included.")
else: print("You have chosen to run FinalFit without systematics. If this was an error please use '--doSystematics'")
if (args_count := len(sys.argv)) > 1:
    config_file = sys.argv[1]
    if os.path.exists( config_file ):
        with open(config_file) as f:
            cfg = yaml.safe_load(f)
else:
  print("[ERROR] no config file was specified. Leaving...")
  leave()

year = cfg["common"]["year"]
ext = cfg["common"]["extension"]
input_dir = os.path.join(os.getcwd(),cfg["tree"]["input_dir"])
cat = cfg["common"]["cats"]
MLow, MHigh = cfg["common"]["binning"].split(",")
MBins = cfg["common"]["nbins"]
MNom = cfg["signal"]["mass_points"].split(",")[len(cfg["signal"]["mass_points"].split(","))//2]

# Create config files for each step
# For signal
py_config = f"""_year = '{year}'

signalScriptCfg = {{

  # Setup
  'inputWSDir': '{input_dir}/{cfg["signal"]["dir"]}/ws_{cfg["signal"]["procs"]}/',
  'procs': '{cfg["signal"]["procs"]}', # if auto: inferred automatically from filenames
  'cats': '{cat}', # if auto: inferred automatically from workspace
  'ext': '{ext}',
  'analysis': 'highMassAnalysis', # To specify which replacement dataset mapping (defined in ./python/replacementMap.py)
  'year': '%s' % _year, # Use 'combined' if merging all years: not recommended
  'width': '{cfg["signal"]["width"]}',
  'massPoints': '{cfg["signal"]["mass_points"]}',

  # Photon shape systematics
  'scales': 'Scale', # separate nuisance per year
  'scalesCorr': '', # correlated across years
  'scalesGlobal': '', # affect all processes equally, correlated across years
  'smears': 'Smearing', # separate nuisance per year

  # Job submission options
  'batch': 'local', # ['condor','SGE','IC','local']
  'queue': 'espresso',

}}
"""

with open("Signal/config_high_mass.py", "w") as f:
    f.write(py_config)

with open("Interference/config_high_mass.py", "w") as f:
    f.write(py_config)


# For background
py_config = f"""

backgroundScriptCfg = {{

  # Setup
  'inputWS': '{input_dir}/data/ws/{cfg["background"]["file"]}',  # Input RooWorkspace

  'cats': '{cfg["common"]["cats"]}',
  'catOffset': 0,         # No offset needed since this is a single file
  'ext': '{ext}', # Will name output directories like: outputs_fTest_2022inclusive/
  'year': '{year}',         # Shown in plots; adjust if merging multiple years

  # Job submission options
  'batch': 'local',       # You can change to 'condor' if you'd prefer distributed batch mode
  'queue': 'espresso'     # Ignored when batch is 'local'

}}
"""

with open("Background/config_high_mass.py", "w") as f:
    f.write(py_config)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
os.chdir("Trees2WS")
if len(opt.run_only) == 0 or "trees" in opt.run_only:
    # Run Tree2Workspace
    syst_opt = '--doSystematics' if opt.do_syst else ''
    dir = cfg["signal"]["dir"]

    cmd = f"python3 RunWSScripts.py --inputDir {input_dir}/{dir}/ --inputConfig config_high_mass.py \
        --year {year} --mode trees2ws --batch local --modeOpts \"--minMass {MLow} --maxMass {MHigh} {syst_opt} \""
    run(cmd)

    cmd = f"python3 RunWSScripts.py --inputDir {input_dir}/data/ --inputConfig config_high_mass.py \
        --year {year} --mode trees2ws_data --batch local --modeOpts \"--applyMassCut --massCutRange {MLow},{MHigh} \""
    run(cmd)

if len(opt.run_only) == 0 or "signal" in opt.run_only:
    # Run Signal fit
    os.chdir("../Signal")

    if opt.do_syst:
        cmd = f"""python3 RunSignalScripts.py --inputConfig config_high_mass.py --mode calcPhotonSyst \
            --modeOpts \" --nBins {MBins}  --minMass {MLow} --maxMass {MHigh}\""""
        run(cmd)

    syst_opt = '' if opt.do_syst else '--skipSystematics'

    cmd = f"""python3 RunSignalScripts.py --inputConfig config_high_mass.py --mode signalFit \
        --modeOpts \" --doPlots {syst_opt} --skipVertexScenarioSplit --skipBeamspotReweigh --nBins {MBins}  --minMass {MLow} --maxMass {MHigh} {cfg['signal']['options']} \""""
    run(cmd)

    # --outputExt tags the packaged workspace's dir AND filename with {ext}
    cmd = f"""python3 RunPackager.py --cats {cat} --exts {ext} --outputExt {ext} --year {year} \
        --massPoints {cfg["signal"]["mass_points"]} --batch local"""
    run(cmd)

if (len(opt.run_only) == 0 or "interference" in opt.run_only) and not (opt.skip_intf):
    # Run Background fit
    os.chdir("../Interference")

    bins = np.concatenate([
        np.linspace(   0,  200, 51)[:-1],
        np.linspace( 200,  500, 31)[:-1],
        np.linspace( 500, 1000, 11)[:-1],
        np.linspace(1000, 2000, 11)[:-1],
        np.linspace(2000, 5000, 16),
    ])
    mPoints = ','.join([str(element) for element in bins])

    cmd = f"""python3 computeGGBoxEff.py --config tools/{year}_cfg.yaml \
        --massList {mPoints} --outCsv tools/csv/ggbox_eff_{year}_{cat}_09_09.csv"""
    print(cmd)
    run(cmd)

    # cmd = f"""python3 RunInterferenceScripts.py --inputConfig config_high_mass.py --mode computeIntf \
    #     --modeOpts \"  --minMass {MLow} --maxMass {MHigh} \""""
    # run(cmd)

if len(opt.run_only) == 0 or "background" in opt.run_only:
    # Run Background fit
    os.chdir("../Background")

    cmd = f"python3 RunBackgroundScripts.py --inputConfig config_high_mass.py --mode fTestParallel"
    run(cmd)

if len(opt.run_only) == 0 or "datacard" in opt.run_only:
    # Run Datacard maker
    os.chdir("../Datacard")
    syst_opt = '--doSystematics' if opt.do_syst else ''
    intf_opt = '--skipIntf' if opt.skip_intf else ''
    dir = cfg["signal"]["dir"]

    cmd = f"""python3 RunYields.py --inputWSDirMap {year}={input_dir}/{dir}/ws_{cfg["signal"]["procs"]} {syst_opt} {intf_opt}\
        --cats {cat} --procs {cfg["signal"]["procs"]} --ext {ext} --sigModelExt {ext} --skipCOWCorr --batch local --mass {MNom} --width {cfg["signal"]["width"]}"""
    print("------>", cmd)
    run(cmd)

    cmd = f"""python3 makeDatacard.py --ext {ext} --years {year} --prune {syst_opt} \
        --skipCOWCorr --doMCStatUncertainty --saveDataFrame --output Datacard_{ext} --mass {MNom}"""
    print("------>", cmd)
    run(cmd)

if len(opt.run_only) == 0 or "combine" in opt.run_only:
    # Each subrange/extension gets its own Combine/<ext>/ area
    combine_dir = os.path.join(base_dir, "Combine", ext)
    os.makedirs(os.path.join(combine_dir, "Models", "signal"), exist_ok=True)
    os.makedirs(os.path.join(combine_dir, "Models", "background"), exist_ok=True)
    if not opt.skip_intf:
        os.makedirs(os.path.join(combine_dir, "Models", "interference"), exist_ok=True)
    os.chdir(combine_dir)

    cat_list = cat.split(",")

    for c in cat_list:
        run(f"cp ../../Background/outdir_{ext}/CMS-HGG_multipdf_{c}_{year}*.root Models/background/")
    if not opt.skip_intf:
        for c in cat_list:
            run(f"cp ../../Interference/outdir_{ext}/computeIntf/output/CMS-HGG_intfm_{c}_{year}*.root Models/interference/")
    run(f"cp ../../Datacard/Datacard_{ext}.txt .")

    packaged_files = []
    for c in cat_list:
        packaged_files += glob.glob(f"../../Signal/outdir_{ext}/CMS-HGG_sigfit_{ext}_{c}_{year}*.root")
    if not packaged_files:
        print(f"[WARNING] No packaged signal file found for {ext}/{cat}_{year} in Signal/outdir_{ext}/. Did the signal step run?")
    for src in packaged_files:
        shutil.copy(src, os.path.join("Models", "signal", os.path.basename(src)))

if len(opt.run_only) == 0 or "text2ws" in opt.run_only:
    # Build the physics-model workspace (mu_inclusive by default) from the datacard
    os.chdir(os.path.join(base_dir, "Combine", ext))

    combine_mode = cfg.get("combine", {}).get("mode", "mu_inclusive")
    common_opts = cfg.get("combine", {}).get(
        "common_opts", f"-m {MNom} higgsMassRange={MLow},{MHigh}"
    )

    cmd = f"""python3 ../RunText2Workspace.py --mode {combine_mode} --batch local \
        --ext _{ext} --common_opts \"{common_opts}\""""
    print("------>", cmd)
    run(cmd)

if len(opt.run_only) == 0 or "limits" in opt.run_only:
    os.chdir(os.path.join(base_dir, "Combine", ext))

    combine_mode = cfg.get("combine", {}).get("mode", "mu_inclusive")
    workspace = f"Datacard_{ext}_{combine_mode}.root"
    title = cfg.get("limits", {}).get("title", f"{cfg['signal']['procs']}, {year}")
    mass_points = cfg.get("limits", {}).get("mass_points", cfg["signal"]["mass_points"])

    cmd = f"""python3 ../RunLimits.py {workspace} --outdir Limits --extension {ext} \
        --mass_points {mass_points} --title \"{title}\""""
    print("------>", cmd)
    run(cmd)

if len(opt.run_only) == 0 or "collect" in opt.run_only:
    # Gather everything (plots, workspaces, datacard, limits).
    results_dir = os.path.join(base_dir, "Results", ext)
    os.makedirs(results_dir, exist_ok=True)

    def collect(src, dst_name):
        src = os.path.join(base_dir, src)
        if os.path.isdir(src):
            shutil.copytree(src, os.path.join(results_dir, dst_name), dirs_exist_ok=True)
        elif os.path.exists(src):
            shutil.copy(src, os.path.join(results_dir, dst_name))
        else:
            print(f"[WARNING] Could not collect '{src}': not found. Skipping.")

    collect(f"Signal/outdir_{ext}", "signal")
    collect(f"Background/outdir_{ext}", "background")
    if not opt.skip_intf:
        collect(f"Interference/outdir_{ext}", "interference")
    collect(f"Datacard/Datacard_{ext}.txt", "Datacard.txt")
    collect(f"Combine/{ext}", "combine")

    print(f" --> Collected results for '{ext}' into {results_dir}")