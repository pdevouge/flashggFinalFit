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
from optparse import OptionParser

import yaml

FINALFIT_DIR = os.path.dirname(os.path.abspath(__file__))

FINALFIT_YAML_TEMPLATE = """\
############################
##   Common Ingredients   ##
############################
common:
  extension: "{extension}"
  binning: "{binning}"
  nbins: "{nbins}"
  year: "{year}"
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
  mass_points: "{limits_mass_points}"
"""


def get_options():
    parser = OptionParser(usage="usage: %prog campaign.yaml [options]")
    parser.add_option("--only", dest="only", default="",
                       help="Comma separated list of subrange names to run (default: all)")
    parser.add_option("--stages", dest="stages", default="",
                       help="Passed through to runFinalfit.py --runOnly (default: run everything)")
    parser.add_option("--doSystematics", dest="do_syst", action="store_true",
                       help="Run FinalFit with systematics (default: systematics off)")
    parser.add_option("--skipIntf", dest="skip_intf", action="store_true",
                       help="Skip interference modelling, passed through to runFinalfit.py")
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
    yaml_text = FINALFIT_YAML_TEMPLATE.format(
        extension=extension,
        binning=subrange["binning"],
        nbins=subrange["nbins"],
        year=common["year"],
        cats=subrange.get("cats", common["cats"]),
        tree_input_dir=common["tree_input_dir"],
        signal_dir=subrange["signal_dir"],
        procs=common["procs"],
        mass_points=subrange["mass_points"],
        width=common["width"],
        options=common.get("options", ""),
        background_file=subrange.get("background_file", common["background_file"]),
    )

    # Optional: run limits on a different mass-point list than the one the
    # signal model was fit on (e.g. a finer scan for the limit plot). Falls
    # back to signal.mass_points in runFinalfit.py if omitted here.
    limits_mass_points = subrange.get("limits_mass_points") or common.get("limits_mass_points")
    if limits_mass_points:
        yaml_text += FINALFIT_YAML_LIMITS_BLOCK.format(limits_mass_points=limits_mass_points)

    out_dir = os.path.join(FINALFIT_DIR, "generated_configs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{extension}.yaml")
    with open(out_path, "w") as f:
        f.write(yaml_text)
    print(f" --> Wrote {out_path}")
    return out_path


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

    print(f"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CAMPAIGN: {cfg['campaign']} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f" --> Subranges to run: {[s['name'] for s in subranges]}")
    print(" --> Assuming signal_X-Y/ and data/ are already populated "
          "(run prepare_campaign_inputs.py under dnavenv first if not)")

    for subrange in subranges:
        extension = f"{cfg['campaign']}_{subrange['name'].replace('-', '_')}"
        print(f"\n==================== Subrange '{subrange['name']}' (ext={extension}) ====================")

        config_path = write_finalfit_config(cfg, subrange, extension)

        ff_cmd = f"python3 runFinalfit.py {config_path}"
        if opt.do_syst:
            ff_cmd += " --doSystematics"
        if opt.skip_intf:
            ff_cmd += " --skipIntf"
        if opt.stages:
            ff_cmd += f" --runOnly {opt.stages}"

        ret = run(ff_cmd, opt.dry_run, cwd=FINALFIT_DIR)
        if ret != 0 and not opt.dry_run:
            print(f"[ERROR] runFinalfit.py failed for subrange '{subrange['name']}' (exit code {ret}). Stopping.")
            sys.exit(ret)

        print(f" --> Done: Results/{extension}/")

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CAMPAIGN DONE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


if __name__ == "__main__":
    main()
