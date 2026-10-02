# Orchestrates runFinalfit.py (signal/background/datacard/combine/limits)
# for every mass subrange of a campaign in a single command
#
# This assumes signal_X-Y/ and data/ under the campaign's tree_input_dir are
# already populated. First run:
#
#   HiggsDNA environment
#   python3 prepare_campaign_inputs.py campaigns/spin0_09_27.yaml
#
#   new shell, or `deactivate` then cmsenv/source setup.sh
#   python3 run_campaign.py campaigns/spin0_09_27.yaml
#
# Usage:
#   python3 run_campaign.py campaigns/spin0_09_27.yaml
#   python3 run_campaign.py campaigns/spin0_09_27.yaml --only 500-1000
#   python3 run_campaign.py campaigns/spin0_09_27.yaml --dry-run
#
# Each subrange ends up in Results/<campaign>_<subrange>/ (plots, packaged
# signal workspace, datacard, Combine workspace, limits json+plot).
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
from optparse import OptionParser

import yaml

FINALFIT_DIR = os.path.dirname(os.path.abspath(__file__))

# setup.sh puts tools/ on PYTHONPATH, but don't depend on it having been sourced.
sys.path.insert(0, os.path.join(FINALFIT_DIR, "tools"))
from campaignSubmissionTools import submit_campaign

FINALFIT_YAML_TEMPLATE = """\
############################
##   Common Ingredients   ##
############################
common:
  extension: "{extension}"
  binning: "{binning}"
  nbins: "{nbins}"
{year_lines}
  cats: "{cats}"

############################
##    Trees 2 workspace   ##
############################
tree:
  input_dir: "{tree_input_dir}"

############################
##         Signal         ##
############################
signal:
  dir: "{signal_dir}"
  procs: "{procs}"
  mass_points: "{mass_points}"
  width: "{width}"
  options: "{options}"

############################
##       Background       ##
############################
background:
  file: "{background_file}"
"""

FINALFIT_YAML_LIMITS_BLOCK = """
############################
##         Limits         ##
############################
limits:
{limits_lines}
"""

FINALFIT_YAML_INTERFERENCE_BLOCK = """
############################
##      Interference      ##
############################
interference:
  catDict: "{catdict}"
"""


def get_options():
    parser = OptionParser(usage="usage: %prog campaign.yaml [options]")
    parser.add_option("--only", dest="only", default="",
                       help="Comma separated list of subrange names to run (default: all)")
    parser.add_option("--stages", dest="stages", default="",
                       help="Passed through to runFinalfit.py --runOnly (default: run everything)")
    parser.add_option("--jobs", dest="jobs", default=1, type="int",
                       help="Number of subranges to run at the same time (default: 1, i.e. one after another)" 
                       "Use 0 to run all subrange locally at the same time.")
    parser.add_option("--doSystematics", dest="do_syst", action="store_true",
                       help="Run FinalFit with systematics (default: systematics off)")
    parser.add_option("--skipIntf", dest="skip_intf", action="store_true",
                       help="Skip interference modelling, passed through to runFinalfit.py")
    parser.add_option("--limitJobs", dest="limit_jobs", default=1, type="int",
                       help="Number of mass points to run at the same time inside each subrange's")
    parser.add_option("--impactJobs", dest="impact_jobs", default=1, type="int",
                       help="Number of nuisance-parameter fits to run at the same time inside each "
                       "subrange's impacts stage")
    parser.add_option("--lumiscale", dest="lumiscale", default=None, type="float",
                       help="Passed through to runFinalfit.py: also produce a lumi-projected limit "
                       "(frozen rateParam * * <factor> on top of the nominal datacard) into Limits_lumiscale/, "
                       "alongside the nominal Limits/.")
    parser.add_option("--batch", dest="batch", default="local",
                       choices=["local", "condor"],
                       help="local (default) runs the campaign here; condor submits the WHOLE campaign.")
    parser.add_option("--flavour", dest="flavour", default="tomorrow",
                       help="HTCondor +JobFlavour for --batch condor (default: tomorrow = 1 day).")
    parser.add_option("--condor-dir", dest="condor_dir", default="condor_campaign",
                       help="Where to write the condor .sh/.sub and their logs, relative to this campaign's.")
    parser.add_option("--dry-run", dest="dry_run", action="store_true",
                       help="Print every command instead of running it")
    return parser.parse_args()


def run(cmd, dry_run, cwd=None):
    print(f"\n$ {cmd}\n")
    if dry_run:
        return 0
    return subprocess.call(cmd, shell=True, cwd=cwd)


def write_finalfit_config(cfg, subrange, extension):
    common = cfg["common"]
    years = common["years"]
    year_lines = '  years: "%s"' % ",".join(years)
    if len(years) > 1:
        # Multi-year (mergeYears) only: tag for the single combined background fit
        year_lines += '\n  merged_label: "%s"' % common["merged_label"]

    yaml_text = FINALFIT_YAML_TEMPLATE.format(
        extension=extension,
        binning=subrange["binning"],
        nbins=subrange["nbins"],
        year_lines=year_lines,
        cats=subrange.get("cats", common["cats"]),
        tree_input_dir=common["tree_input_dir"],
        signal_dir=subrange["signal_dir"],
        procs=common["procs"],
        mass_points=subrange["mass_points"],
        width=common["width"],
        options=common.get("options", ""),
        background_file=subrange.get("background_file", common["background_file"]),
    )

    # interference step has to categorize the GGBox background so it needs this subrange's catDict too
    cat_dict = subrange.get("catDict", cfg.get("signal_postprocessing", {}).get("catDict"))
    if cat_dict:
        yaml_text += FINALFIT_YAML_INTERFERENCE_BLOCK.format(catdict=cat_dict)

    # Optional: run limits on a different mass-point list than the one the
    # signal model was fit on (e.g. a finer scan for the limit plot). Falls
    # back to signal.mass_points in runFinalfit.py if omitted here.
    limits_lines = []
    limits_mass_points = subrange.get("limits_mass_points") or common.get("limits_mass_points")
    if limits_mass_points:
        limits_lines.append('  mass_points: "%s"' % limits_mass_points)
    # Optional: force the luminosity label on the limit plot. Normally left unset (derived auto year via lumiMap)
    limits_lumi = subrange.get("limits_lumi") or common.get("limits_lumi")
    if limits_lumi:
        limits_lines.append('  lumi: "%s"' % limits_lumi)
    if limits_lines:
        yaml_text += FINALFIT_YAML_LIMITS_BLOCK.format(limits_lines="\n".join(limits_lines))

    out_dir = os.path.join(FINALFIT_DIR, "generated_configs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{extension}.yaml")
    with open(out_path, "w") as f:
        f.write(yaml_text)
    print(f" --> Wrote {out_path}")
    return out_path


def run_subrange(cfg, subrange, opt, capture):
    """Run one subrange end to end. Many can be run at the same time."""
    extension = f"{cfg['campaign']}_{subrange['name'].replace('-', '_')}"
    header = f"\n==================== Subrange '{subrange['name']}' (ext={extension}) ===================="
    if not capture:
        print(header)

    config_path = write_finalfit_config(cfg, subrange, extension)

    ff_cmd = f"python3 runFinalfit.py {config_path}"
    if opt.do_syst:
        ff_cmd += " --doSystematics"
    if opt.skip_intf:
        ff_cmd += " --skipIntf"
    if opt.stages:
        ff_cmd += f" --runOnly {opt.stages}"
    if opt.limit_jobs != 1:
        ff_cmd += f" --limitJobs {opt.limit_jobs}"
    if opt.impact_jobs != 1:
        ff_cmd += f" --impactJobs {opt.impact_jobs}"
    if opt.lumiscale:
        ff_cmd += f" --lumiscale {opt.lumiscale}"

    if not capture:
        return extension, run(ff_cmd, opt.dry_run, cwd=FINALFIT_DIR), None

    log_dir = os.path.join(FINALFIT_DIR, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f"run_{extension}.log")
    print(f" --> Started subrange '{subrange['name']}' (ext={extension}) -> {log_path}")
    with open(log_path, "w") as lf:
        lf.write(f"$ {ff_cmd}\n\n")
        lf.flush()
        ret = subprocess.call(ff_cmd, shell=True, cwd=FINALFIT_DIR, stdout=lf, stderr=subprocess.STDOUT)
    return extension, ret, log_path


def main():
    (opt, args) = get_options()
    if len(args) != 1:
        print("[ERROR] Usage: run_campaign.py campaign.yaml [options]")
        sys.exit(1)

    with open(args[0]) as f:
        cfg = yaml.safe_load(f)

    only = set(x.strip() for x in opt.only.split(",") if x.strip())
    subranges = [s for s in cfg["subranges"] if not only or s["name"] in only]
    if not subranges:
        print(f"[ERROR] No subranges matched --only={opt.only!r}")
        sys.exit(1)

    years = cfg["common"].get("years")
    if not years or isinstance(years, str):
        print("[ERROR] campaign yaml common block must set 'years' as a YAML list, e.g. [\"2022postEE\"] for a single era.")
        sys.exit(1)
    if opt.limit_jobs < 1:
        print(f"[ERROR] --limitJobs must be >= 1 (got {opt.limit_jobs}).")
        sys.exit(1)
    if opt.impact_jobs < 1:
        print(f"[ERROR] --impactJobs must be >= 1 (got {opt.impact_jobs}).")
        sys.exit(1)
    # TODO remove when interference supports merge years
    if len(years) > 1 and not opt.skip_intf:
        print("[ERROR] Multi-year (mergeYears) campaigns don't support interference modelling yet "
              "-- pass --skipIntf.")
        sys.exit(1)

    print(f"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CAMPAIGN: {cfg['campaign']} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    if opt.batch == "condor":
        # One HTCondor job per subrange, each on its own worker.
        if opt.jobs != 1:
            print(" --> --jobs is ignored with --batch condor: every subrange already gets its own "
                  "dedicated worker, one job each.")
        if opt.limit_jobs != 1:
            print(f" --> Each subrange job will request {opt.limit_jobs} CPUs and run its limits "
                  f"stage {opt.limit_jobs} mass points at a time.")
        if opt.impact_jobs != 1:
            print(f" --> Each subrange job will request {opt.impact_jobs} CPUs and run its impacts "
                  f"stage {opt.impact_jobs} nuisance-parameter fits at a time.")
        sys.exit(submit_campaign(cfg, os.path.abspath(args[0]), opt, FINALFIT_DIR))
    print(f" --> Subranges to run: {[s['name'] for s in subranges]}")
    print(" --> Assuming signal_X-Y/ and data/ are already populated "
          "(run prepare_campaign_inputs.py under dnavenv first if not)")

    n_jobs = len(subranges) if opt.jobs == 0 else max(1, opt.jobs)
    if n_jobs == 1 or opt.dry_run or len(subranges) == 1:
        for subrange in subranges:
            extension, ret, _ = run_subrange(cfg, subrange, opt, capture=False)
            if ret != 0 and not opt.dry_run:
                print(f"[ERROR] runFinalfit.py failed for subrange '{subrange['name']}' (exit code {ret}). Stopping.")
                sys.exit(ret)
            print(f" --> Done: Results/{extension}/")
    else:
        print(f" --> Running {len(subranges)} subranges with up to {n_jobs} at a time; "
              f"per-subrange output goes to logs/run_<ext>.log")
        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            results = list(executor.map(lambda sr: run_subrange(cfg, sr, opt, capture=True), subranges))
        failed = []
        for subrange, (extension, ret, log_path) in zip(subranges, results):
            if ret == 0:
                print(f" --> Done: Results/{extension}/   (log: {log_path})")
            else:
                print(f"[ERROR] runFinalfit.py failed for subrange '{subrange['name']}' "
                      f"(exit code {ret}) -- see {log_path}")
                failed.append(subrange["name"])
        if failed:
            print(f"[ERROR] {len(failed)} of {len(subranges)} subranges failed: {failed}")
            sys.exit(1)

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CAMPAIGN DONE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


if __name__ == "__main__":
    main()
