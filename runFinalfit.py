# Script for running the entire pipeline for flashggFinalFit
import os, sys, yaml, glob, shutil, subprocess
import numpy as np
from optparse import OptionParser
from collections import OrderedDict as od
# setup.sh puts tools/ on PYTHONPATH, but don't depend on it having been sourced
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools"))
from commonObjects import lumiMap

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING FINALFIT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

def get_options():
  parser = OptionParser()
  parser.add_option('--runOnly', dest='run_only', default='', help="Run only given steps (trees, signal, background, datacard, combine, text2ws, limits, collect)")
  parser.add_option('--doSystematics', dest='do_syst', action='store_true', help="Run with systematics")
  parser.add_option('--skipIntf', dest='skip_intf', action='store_true', help="Skip interference making")
  parser.add_option('--limitJobs', dest='limit_jobs', default=None, type='int', help="Number of mass points to run at the same time in the limits stage (overrides limits.jobs in the config)")
  parser.add_option('--onlyYear', dest='only_year', default='', help="Run just this one era out of a multi-era config, without merging. The era must be listed in common.years. Requires --backgroundFile, since a multi-era campaign stages the MERGED background.")
  parser.add_option('--backgroundFile', dest='background_file', default='', help="Override background.file from the config (name of the file under <tree_input_dir>/data/). Needed with --onlyYear to point at that era's own data instead of the merged all-era file.")
  parser.add_option('--dry-run', dest='dry_run', action='store_true', help="Print every command instead of running it")
  return parser.parse_args()
(opt,args) = get_options()

base_dir = os.getcwd()

def fail(msg):
  print(f"[ERROR] {msg}")
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING FINALFIT (FAILED) ~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  sys.exit(1)

def run(cmd):
  print(f"\n$ {cmd}\n")
  if opt.dry_run:
    return 0
  ret = subprocess.call(cmd, shell=True)
  if ret < 0:
    print(f"[ERROR] Command was interrupted by signal {-ret}. Stopping runFinalfit.py.")
    sys.exit(1)
  elif ret != 0:
    print(f"[ERROR] Command failed (exit code {ret}). Stopping runFinalfit.py.")
    sys.exit(ret)
  return ret

# --dry-run must not touch the filesystem.
def write_file(path, text):
  if opt.dry_run:
    print(f"[DRY RUN] Would write {path}")
    return
  with open(path, "w") as f:
    f.write(text)

def copy_file(src, dst):
  if opt.dry_run:
    print(f"[DRY RUN] Would copy {src} -> {dst}")
    return
  shutil.copy(src, dst)

def ensure_dir(path):
  if opt.dry_run:
    if not os.path.isdir(path): print(f"[DRY RUN] Would create {path}/")
    return
  os.makedirs(path, exist_ok=True)

def chdir(path):
  # Under --dry-run the Combine/<ext> area may not exist yet
  if opt.dry_run and not os.path.isdir(path):
    print(f"[DRY RUN] Would cd {path}")
    return
  os.chdir(path)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract options from config file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXTRACTING CONFIG FROM YAML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
if opt.do_syst: print("-----> 'doSystematics': FinalFit will run with systematics included.")
else: print("You have chosen to run FinalFit without systematics. If this was an error please use '--doSystematics'")
if (args_count := len(sys.argv)) > 1:
    config_file = sys.argv[1]
    if os.path.exists( config_file ):
        with open(config_file) as f:
            cfg = yaml.safe_load(f)
    else:
        fail(f"config file '{config_file}' does not exist.")
else:
  fail("no config file was specified.")

# Eras come from `years`, whatever the count: run_campaign.py always writes it,
# single era included, and merge_years is decided on its length alone -- the same
# rule the campaign layer uses. Accepted either as a comma string (what
# run_campaign.py writes, matching cats/mass_points) or as a YAML list (how
# campaigns/*.yaml spells it), so pointing this at either file works.
if "years" in cfg["common"]:
    raw = cfg["common"]["years"]
    years = [y.strip() for y in (raw.split(",") if isinstance(raw, str) else raw)]
else:
    fail("common must set 'years' (e.g. years: \"2022postEE\" or years: [\"2022preEE\", \"2022postEE\"]).")

if not years or not all(years):
    fail(f"common.years is empty or contains a blank era: {cfg['common'].get('years')!r}")

era_qualified = len(years) > 1

if opt.only_year:
    if opt.only_year not in years:
        fail(f"--onlyYear '{opt.only_year}' is not in common.years ({','.join(years)}).")
    years = [opt.only_year]
    merge_years = False

else:
    merged_label = cfg["common"].get("merged_label")
    merge_years = len(years) > 1 and bool(merged_label)

    if len(years) > 1 and not merged_label:
        fail(f"common.years lists {len(years)} eras ({','.join(years)}) but no merged_label. "
            f"Set merged_label to merge them into one datacard, or pass "
            f"--onlyYear <era> to run a single one.")

# A multi-era campaign usually stages the MERGED background (allData_<merged_label>.root),
# so running one era against it would fit that era's signal to all-era data and
# label the plot with that era's luminosity. Fail instead of quietly doing it.
bkg_file = opt.background_file or cfg.get("background", {}).get("file")
if not bkg_file:
    fail("background.file is not set in the config, and no --backgroundFile was given.")
if opt.only_year and era_qualified and not opt.background_file:
    fail(f"--onlyYear {opt.only_year} needs --backgroundFile: this config's "
         f"background.file ({cfg['background']['file']}) is the merged all-era dataset, and "
         f"fitting one era's signal against it would mismatch the luminosity. Pass that era's "
         f"own allData file (staged under {cfg['tree']['input_dir']}/data/).")
if opt.background_file:
    print(f" --> Background file overridden: {cfg.get('background', {}).get('file')} -> {bkg_file}")

# TODO: interference merge years
if merge_years and not opt.skip_intf:
    fail("Multi-year (mergeYears) campaigns don't support interference modelling yet. Pass --skipIntf.")

# Single label used for the one combined background fit and for Combine/limits naming. 
bkg_year = merged_label if merge_years else years[0]

# Integrated luminosity shown on the limit plot. Override with `limits: lumi: "..."` in the config.
def get_lumi_label():
    override = cfg.get("limits", {}).get("lumi")
    if override not in (None, ""):
        return str(override)
    if bkg_year not in lumiMap:
        fail(f"No luminosity for '{bkg_year}' in lumiMap (tools/commonObjects.py). "
             f"Add it there, or set `limits: lumi: \"<value>\"` in {config_file}.")
    return "%.2f" % lumiMap[bkg_year]

lumi_label = get_lumi_label()
print(f" --> Luminosity label for the limit plot: {lumi_label} fb^-1 (lumiMap['{bkg_year}'])")

ext = cfg["common"]["extension"]
config_name = f"config_{ext}.py"
# Configs now live in each section's config/ subdir; scripts run with cwd = the section dir.
config_rel = f"config/{config_name}"
input_dir = os.path.join(os.getcwd(),cfg["tree"]["input_dir"])
cat = cfg["common"]["cats"]
MLow, MHigh = cfg["common"]["binning"].split(",")
MBins = cfg["common"]["nbins"]
MNom = cfg["signal"]["mass_points"].split(",")[len(cfg["signal"]["mass_points"].split(","))//2]

# Create config files for each step
# For signal: written per-era
def write_signal_config(year, ext_tag, ws_dir):
    py_config = f"""_year = '{year}'

signalScriptCfg = {{

  # Setup
  'inputWSDir': '{ws_dir}',
  'procs': '{cfg["signal"]["procs"]}', # if auto: inferred automatically from filenames
  'cats': '{cat}', # if auto: inferred automatically from workspace
  'ext': '{ext_tag}',
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
    write_file(os.path.join(base_dir, "Signal", "config", config_name), py_config)

# TODO: Interference logic related
if not merge_years:
    write_signal_config(years[0], ext, f'{input_dir}/{cfg["signal"]["dir"]}/ws_{cfg["signal"]["procs"]}/')
    copy_file(os.path.join(base_dir, "Signal", "config", config_name),
              os.path.join(base_dir, "Interference", "config", config_name))


# For background
py_config = f"""

backgroundScriptCfg = {{

  # Setup
  'inputWS': '{input_dir}/data/ws/{bkg_file}',  # Input RooWorkspace

  'cats': '{cfg["common"]["cats"]}',
  'catOffset': 0,         # No offset needed since this is a single file
  'ext': '{ext}', # Will name output directories like: outputs_fTest_2022inclusive/
  'year': '{bkg_year}',         # Shown in plots; already the combined-era label when merge_years

  # Job submission options
  'batch': 'local',       # You can change to 'condor' if you'd prefer distributed batch mode
  'queue': 'espresso'     # Ignored when batch is 'local'

}}
"""

write_file(os.path.join(base_dir, "Background", "config", config_name), py_config)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Trees2ws step
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
chdir("Trees2WS")
if len(opt.run_only) == 0 or "trees" in opt.run_only:
    # Run Tree2Workspace
    syst_opt = '--doSystematics' if opt.do_syst else ''
    dir = cfg["signal"]["dir"]

    for year in years:
        # Each era's signal trees live in their own <dir>_<year>
        print(f"\n---- trees2ws for era '{year}' ----")
        year_suffix = f"_{year}" if era_qualified else ""
        cmd = f"python3 RunWSScripts.py --inputDir {input_dir}/{dir}{year_suffix}/ --inputConfig config/config_high_mass.py \
            --year {year} --ext {ext} --mode trees2ws --batch local --modeOpts \"--minMass {MLow} --maxMass {MHigh} {syst_opt} \""
        run(cmd)

    # Convert ONLY 'inputPattern'.root background file. <tree_input_dir>/data/ is shared
    cmd = f"python3 RunWSScripts.py --inputDir {input_dir}/data/ --inputPattern {bkg_file} --inputConfig config/config_high_mass.py \
        --year {bkg_year} --ext {ext} --mode trees2ws_data --batch local --modeOpts \"--applyMassCut --massCutRange {MLow},{MHigh} \""
    run(cmd)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Signal step
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if len(opt.run_only) == 0 or "signal" in opt.run_only:
    # Run Signal fit, once per era.
    chdir("../Signal")

    ext_years = []
    for year in years:
        ext_year = f"{ext}_{year}" if era_qualified else ext
        ext_years.append(ext_year)
        if era_qualified: # more than one year in the list
            print(f"\n---- Signal fit for era '{year}' (ext={ext_year}) ----")
            ws_dir = f'{input_dir}/{cfg["signal"]["dir"]}_{year}/ws_{cfg["signal"]["procs"]}/'
            write_signal_config(year, ext_year, ws_dir)

        if opt.do_syst:
            cmd = f"""python3 RunSignalScripts.py --inputConfig {config_rel} --mode calcPhotonSyst \
                --modeOpts \" --nBins {MBins}  --minMass {MLow} --maxMass {MHigh}\""""
            run(cmd)

        syst_opt = '' if opt.do_syst else '--skipSystematics'

        cmd = f"""python3 RunSignalScripts.py --inputConfig {config_rel} --mode signalFit \
            --modeOpts \" --doPlots {syst_opt} --skipVertexScenarioSplit --skipBeamspotReweigh --nBins {MBins}  --minMass {MLow} --maxMass {MHigh} {cfg['signal']['options']} \""""
        run(cmd)

    if merge_years:
        cmd = f"""python3 RunPackager.py --cats {cat} --exts {','.join(ext_years)} --outputExt {ext} --mergeYears \
            --massPoints {cfg["signal"]["mass_points"]} --batch local"""
    else:
        # ext_years[0] rather than ext: under --onlyYear the single fit still
        # lives in results/outdir_<ext>_<era>, so --exts has to carry the era tag.
        cmd = f"""python3 RunPackager.py --cats {cat} --exts {ext_years[0]} --outputExt {ext} --year {years[0]} \
            --massPoints {cfg["signal"]["mass_points"]} --batch local"""
    run(cmd)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Interference step
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (len(opt.run_only) == 0 or "interference" in opt.run_only) and not (opt.skip_intf):
    # Run Background fit
    # TODO: no mergeYears support here yet 
    chdir("../Interference")

    bins = np.concatenate([
        np.linspace(   0,  200, 51)[:-1],
        np.linspace( 200,  500, 31)[:-1],
        np.linspace( 500, 1000, 11)[:-1],
        np.linspace(1000, 2000, 11)[:-1],
        np.linspace(2000, 5000, 16),
    ])
    mPoints = ','.join([str(element) for element in bins])

    cmd = f"""python3 computeGGBoxEff.py --config tools/{years[0]}_cfg.yaml \
        --massList {mPoints} --outCsv tools/csv/ggbox_eff_{years[0]}_{cat}_09_09.csv"""
    run(cmd)

    # cmd = f"""python3 RunInterferenceScripts.py --inputConfig {config_rel} --mode computeIntf \
    #     --modeOpts \"  --minMass {MLow} --maxMass {MHigh} \""""
    # run(cmd)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Background step
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if len(opt.run_only) == 0 or "background" in opt.run_only:
    # Run Background fit
    chdir("../Background")

    cmd = f"python3 RunBackgroundScripts.py --inputConfig {config_rel} --mode fTestParallel"
    run(cmd)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Datacard step
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if len(opt.run_only) == 0 or "datacard" in opt.run_only:
    # Run Datacard maker
    chdir("../Datacard")
    syst_opt = '--doSystematics' if opt.do_syst else ''
    intf_opt = '--skipIntf' if opt.skip_intf else ''
    dir = cfg["signal"]["dir"]

    if era_qualified:
        ws_dir_map = ",".join(f"{y}={input_dir}/{dir}_{y}/ws_{cfg['signal']['procs']}" for y in years)
    else:
        ws_dir_map = f"{years[0]}={input_dir}/{dir}/ws_{cfg['signal']['procs']}"

    if merge_years:
        # Tag for the background model (eg 2022pre+post -> Tag=2022)
        mergeYears_opt = f"--mergeYears --bkgModelTag {bkg_year}"
    else:
        mergeYears_opt = ""

    cmd = f"""python3 RunYields.py --inputWSDirMap {ws_dir_map} {syst_opt} {intf_opt} {mergeYears_opt}\
        --cats {cat} --procs {cfg["signal"]["procs"]} --ext {ext} --sigModelExt {ext} --skipCOWCorr --batch local --mass {MNom} --width {cfg["signal"]["width"]}"""
    run(cmd)

    cmd = f"""python3 makeDatacard.py --ext {ext} --years {','.join(years)} --prune {syst_opt} \
        --skipCOWCorr --doMCStatUncertainty --saveDataFrame --output results/Datacard_{ext} --mass {MNom}"""
    run(cmd)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine step
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if len(opt.run_only) == 0 or "combine" in opt.run_only:
    # Each subrange/extension gets its own Combine/<ext>/ area
    combine_dir = os.path.join(base_dir, "Combine", "results", ext)
    models_sig = os.path.join(combine_dir, "Models", "signal")
    models_bkg = os.path.join(combine_dir, "Models", "background")
    models_intf = os.path.join(combine_dir, "Models", "interference")
    ensure_dir(models_sig)
    ensure_dir(models_bkg)
    if not opt.skip_intf:
        ensure_dir(models_intf)
    chdir(combine_dir)

    cat_list = cat.split(",")

    # Copy everything into Combine dir.
    for c in cat_list:
        src = os.path.join(base_dir, "Background", f"results/outdir_{ext}",
                           f"CMS-HGG_multipdf_{c}_{bkg_year}.root")
        if not os.path.exists(src):
            fail(f"No background model at {src}. Did the background step run for '{bkg_year}'?")
        run(f"cp {src} {models_bkg}/")
    if not opt.skip_intf:
        for c in cat_list:
            intf_glob = os.path.join(base_dir, "Interference", f"results/outdir_{ext}", "computeIntf",
                                     "output", f"CMS-HGG_intfm_{c}_{years[0]}*.root")
            run(f"cp {intf_glob} {models_intf}/")
    run(f"cp {os.path.join(base_dir, 'Datacard', 'results', f'Datacard_{ext}.txt')} {combine_dir}/")

    packaged_files = []
    for c in cat_list:
        name = f"CMS-HGG_sigfit_{ext}_{c}.root" if merge_years else f"CMS-HGG_sigfit_{ext}_{c}_{years[0]}*.root"
        packaged_files += glob.glob(os.path.join(base_dir, "Signal", f"results/outdir_{ext}", name))
    if not packaged_files:
        print(f"[WARNING] No packaged signal file found for {ext}/{cat} in Signal/results/outdir_{ext}/. Did the signal step run?")
    for src in packaged_files:
        copy_file(src, os.path.join(models_sig, os.path.basename(src)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Text2WS step
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if len(opt.run_only) == 0 or "text2ws" in opt.run_only:
    # Build the physics-model workspace (mu_inclusive by default) from the datacard
    chdir(os.path.join(base_dir, "Combine", "results", ext))

    combine_mode = cfg.get("combine", {}).get("mode", "mu_inclusive")
    common_opts = cfg.get("combine", {}).get(
        "common_opts", f"-m {MNom} higgsMassRange={MLow},{MHigh}"
    )

    cmd = f"""python3 ../../RunText2Workspace.py --mode {combine_mode} --batch local \
        --ext _{ext} --common_opts \"{common_opts}\""""
    run(cmd)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Limits step
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if len(opt.run_only) == 0 or "limits" in opt.run_only:
    chdir(os.path.join(base_dir, "Combine", "results", ext))

    combine_mode = cfg.get("combine", {}).get("mode", "mu_inclusive")
    workspace = f"Datacard_{ext}_{combine_mode}.root"
    title = cfg.get("limits", {}).get("title", f"{cfg['signal']['procs']}, {bkg_year}")
    mass_points = cfg.get("limits", {}).get("mass_points", cfg["signal"]["mass_points"])
    limit_jobs = opt.limit_jobs if opt.limit_jobs else cfg.get("limits", {}).get("jobs", 1)

    cmd = f"""python3 ../../RunLimits.py {workspace} --outdir Limits --extension {ext} \
        --mass_points {mass_points} --parallel {limit_jobs} --title \"{title}\" --lumi \"{lumi_label}\""""
    run(cmd)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Collect results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if len(opt.run_only) == 0 or "collect" in opt.run_only:
    # Gather everything (plots, workspaces, datacard, limits).
    results_dir = os.path.join(base_dir, "Results", ext)
    ensure_dir(results_dir)

    def collect(src, dst_name):
        src = os.path.join(base_dir, src)
        dst = os.path.join(results_dir, dst_name)
        if not os.path.exists(src):
            print(f"[WARNING] Could not collect '{src}': not found. Skipping.")
            return
        if opt.dry_run:
            print(f"[DRY RUN] would collect {src} -> {dst}")
            return
        if os.path.isdir(src):
            shutil.copytree(src, dst, dirs_exist_ok=True)
        else:
            shutil.copy(src, dst)

    collect(f"Signal/results/outdir_{ext}", "signal")
    collect(f"Background/results/outdir_{ext}", "background")
    if not opt.skip_intf:
        collect(f"Interference/results/outdir_{ext}", "interference")
    collect(f"Datacard/results/Datacard_{ext}.txt", "Datacard.txt")
    collect(f"Combine/results/{ext}", "combine")

    print(f" --> Collected results for '{ext}' into {results_dir}")