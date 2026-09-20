# Prepares the inputs for a campaign
#
# Usage:
#   python3 prepare_campaign_inputs.py campaigns/spin0_09_27.yaml
#   python3 prepare_campaign_inputs.py campaigns/spin0_09_27.yaml --only 500-1000
#   python3 prepare_campaign_inputs.py campaigns/spin0_09_27.yaml --run-background
#   python3 prepare_campaign_inputs.py campaigns/spin0_09_27.yaml --dry-run
#
# To run prepare_output_file.py on condor instead of locally:
#   python3 prepare_campaign_inputs.py campaigns/spin0_09_27.yaml --run-background \
#       --batch condor --memory 25GB --logs ./logs/
#
# (choices are condor, slurm, slurm/psi, local, local/futures))
#
# When launching job on condor, merge_parquet returns instantly after submission,
# meaning that the root step cannot be run since merged parquets do not exist.
# The script will thus fail, and will need to be rerun once the merged are available:
#
#   python3 prepare_campaign_inputs.py campaigns/spin0_09_27.yaml --run-background \
#       --batch condor --memory 25GB --logs ./logs/          # 1) submits merge jobs
#   # ... wait for the merge jobs to finish ...
#   python3 prepare_campaign_inputs.py campaigns/spin0_09_27.yaml --run-background --root-step \
#       --batch condor --memory 25GB --logs ./logs/          # 2) submits root-conversion jobs
#   # ... wait for those to finish too ...
#   python3 prepare_campaign_inputs.py campaigns/spin0_09_27.yaml --stage-only  # 3) place outputs inside data/signal dir
#
# Per-subrange category overrides: a subrange entry may set its own catDict
# (e.g. different pNN cut boundaries for a different mass window) instead of
# inheriting signal_postprocessing.catDict / background_postprocessing.catDict:
#
#   subranges:
#     - name: "500-1000"
#       catDict: "../HDNA_config_jsons/category_rsg_500-1000.json"
#       cats: "rsg_low_cat,rsg_med_cat,rsg_high_cat"     # optional: overrides common.cats for this subrange
#       ...
#     - name: "1000-2000"
#       catDict: "../HDNA_config_jsons/category_rsg_1000-2000.json"
#       background_file: "allData_PostEE_1000-2000.root"      # required: one file name per catDict
#       background_output_dir: "../RSG_bkg_1000-2000"         # optional: has a sensible default
#       ...
#
# The categories come from catDict, so the data has to be categorized with the
# same catDict as the signal it will be fitted against. That is why background
# postprocessing runs once per distinct catDict rather than once per campaign.
# Each of those runs produces one allData_*.root, which is copied into data/.
#
# What that means when writing the yaml:
#
#   - A subrange with its own catDict needs its own background_file. Two runs
#     writing the same file name would overwrite each other in data/, leaving
#     one subrange fitted against the wrong categorization.
#
#   - background_output_dir is optional. By default each subrange writes to
#     <tree_input_dir>/background_raw_<subrange name>/, which keeps the runs
#     apart and leaves the shared HiggsDNA pool untouched. Set it only if you
#     want the output somewhere specific.
#
#   - Subranges sharing a catDict share one background run, so they must agree
#     on both settings. The default directory is per subrange name, so two such
#     subranges will disagree unless you give them the same explicit
#     background_output_dir.
#
# group_subranges_by_catdict() checks all of this before anything runs, and
# exits with a message naming the subranges that clash.
import os, sys, glob, shutil, subprocess
from optparse import OptionParser

import yaml

FINALFIT_DIR = os.path.dirname(os.path.abspath(__file__))


def get_options():
    parser = OptionParser(usage="usage: %prog campaign.yaml [options]")
    parser.add_option("--only", dest="only", default="",
                       help="Comma separated list of subrange names to run (default: all)")
    parser.add_option("--do-predictions", dest="run_pnn", action="store_true",
                           help="Run pNN prediction maker along with merge and root steps.")
    parser.add_option("--run-background", dest="run_background", action="store_true",
                       help="Actually (re-)run the background/data postprocessing command instead of just "
                            "copying whatever allData_*.root already exists (background doesn't depend on the "
                            "mass subrange, so this normally only needs to be done once, not per campaign run)")
    parser.add_option("--batch", dest="batch", default="",
                       help="Batch system to run prepare_output_file.py on (condor, slurm, slurm/psi, local, local/futures). Default: run locally in the foreground.")
    parser.add_option("--memory", dest="memory", default="",
                       help="Memory allocation for batch jobs (e.g. 25GB).")
    parser.add_option("--logs", dest="logs", default="",
                       help="Log directory for batch job submission.")
    parser.add_option("--root-step", dest="root_step", action="store_true",
                       help="On condor/slurm, submit the root-conversion jobs instead of the merge jobs. Use once a previous merge submission has finished.")
    parser.add_option("--stage-only", dest="stage_only", action="store_true",
                       help="Skip prepare_output_file.py entirely and just stage whatever output already exists, "
                            "into data/ and signal_X-Y/. Use after a condor/slurm run's jobs have finished.")
    parser.add_option("--dry-run", dest="dry_run", action="store_true",
                       help="Print every command instead of running it")
    return parser.parse_args()


def is_real_batch(opt):
    return bool(opt.batch) and opt.batch not in ("local", "local/futures")


def batch_opts(opt):
    if not opt.batch:
        return ""
    parts = [f"--batch {opt.batch}"]
    if opt.memory:
        parts.append(f"--memory {opt.memory}")
    if opt.logs:
        parts.append(f"--logs {opt.logs}")
    return " ".join(parts)


def run(cmd, dry_run, cwd=None):
    print(f"\n$ {cmd}\n")
    if dry_run:
        return 0
    return subprocess.call(cmd, shell=True, cwd=cwd)


def stage_signal_root_files(parquet_dir, signal_dir, dry_run):
    """Copy the output_*.root files produced by prepare_output_file.py into
    the flat signal_X-Y/ dir that FinalFit's Trees2WS step expects.

    They land one directory per sample, at
    <parquet_dir>/root/<sample>/output_*.root
    """
    pattern = os.path.join(parquet_dir, "root", "*", "output_*.root")
    found = sorted(glob.glob(pattern))

    if not found:
        print(f"[WARNING] No 'output_*.root' files found under {parquet_dir} "
              f"(tried: {pattern}). Skipping staging -- you may need to move "
              f"them into {signal_dir} by hand.")
        return

    print(f" --> Found {len(found)} root file(s) matching {pattern}")

    if dry_run:
        for src in found:
            print(f"[DRY RUN] mv {src} -> {signal_dir}/")
        return

    os.makedirs(signal_dir, exist_ok=True)
    for src in found:
        dst = os.path.join(signal_dir, os.path.basename(src))
        print(f" --> mv {src} -> {dst}")
        shutil.move(src, dst)


def group_subranges_by_catdict(bp, common, subranges):
    """Group subranges into one entry per distinct effective catDict, each
    carrying the (background_file, background_output_dir) that group's run
    must use, and validate the collisions described at the top of this file.

    Background/data doesn't depend on the mass *cut* (applied later, per
    subrange, in FinalFit's own trees2ws_data step), but its categorization
    does depend on catDict -- hence once per distinct catDict rather than
    once for the whole campaign.
    """
    # First we create cat_dict -> {"background_file", "output_dir", "subranges": [names]}
    # we essentially check that: same catDict → same file and dir
    groups = {}  
    for subrange in subranges:
        cat_dict = subrange.get("catDict", bp.get("catDict")) # Does subrange has catDict or do we use campaign-wide cats?
        if not cat_dict:
            print(f"[ERROR] subrange '{subrange['name']}' has no catDict -- set it on the subrange itself, "
                  f"or on background_postprocessing.catDict as a campaign-wide default.")
            sys.exit(1)
        bkg_file = subrange.get("background_file", common["background_file"])  # Does subrange has bkg filename or do we use campaign-wide?
        out_dir = subrange.get("background_output_dir", bp.get("output_dir"))  # Does subrange has bkg output dir or do we use campaign-wide?
        if out_dir is None:
            # If no output dir is specified, create one
            # to avoid saving the root file inside the shared space,
            # which could overwrite other files.
            out_dir = os.path.join(common["tree_input_dir"], f"background_raw_{subrange['name']}")
        if out_dir:
            # Resolve against flashggFinalFit/ dir, not the cwd the script happens to be invoked from
            out_dir = os.path.abspath(os.path.join(FINALFIT_DIR, out_dir))

        # Write catDict (with parameters for background postprocessing) in groups (to run later)
        # catDict: (bkg filename, output dir, associated subranges)
        if cat_dict in groups:
            g = groups[cat_dict]
            if g["background_file"] != bkg_file or g["output_dir"] != out_dir:
                print(f"[ERROR] subranges {g['subranges']} and '{subrange['name']}' share catDict {cat_dict!r} "
                      f"but disagree on background_file/background_output_dir -- subranges sharing a catDict "
                      f"must also share the same background output.")
                sys.exit(1)
            g["subranges"].append(subrange["name"]) # simply append subrange name if catDict already in groups
        else:
            groups[cat_dict] = {"background_file": bkg_file, "output_dir": out_dir, "subranges": [subrange["name"]]}

    # Second we check the opposite: same file (or dir) → same catDict
    by_file, by_location = {}, {}
    for cat_dict, g in groups.items():
        if g["background_file"] in by_file and by_file[g["background_file"]] != cat_dict:
            print(f"[ERROR] two different catDicts both use background_file {g['background_file']!r} -- give "
                  f"each distinct catDict its own subrange-level background_file so they don't overwrite "
                  f"each other under data/.")
            sys.exit(1)
        by_file[g["background_file"]] = cat_dict

        location = g["output_dir"] if g["output_dir"] else bp["raw_input_dir"]
        if location in by_location and by_location[location] != cat_dict:
            print(f"[ERROR] two different catDicts both write background output to {location!r} -- give each "
                  f"distinct catDict its own subrange-level background_output_dir so one run's merged/root "
                  f"output doesn't collide with the other's.")
            sys.exit(1)
        by_location[location] = cat_dict

    return groups


def run_background_postprocessing(cfg, opt, subranges):
    bp = cfg.get("background_postprocessing")
    if not bp:
        print(" --> No 'background_postprocessing' block in this campaign -- assuming data/ is already populated.")
        return

    common = cfg["common"]
    data_dir = os.path.abspath(os.path.join(FINALFIT_DIR, common["tree_input_dir"], "data")) # data root file dir.

    # For each catDict and associated parameter (filename, output dir, list of subranges)
    for cat_dict, g in group_subranges_by_catdict(bp, common, subranges).items():
        bkg_file = g["background_file"] # bkg filename
        output_dir = g["output_dir"] # bkg output dir
        data_root_source_dir = output_dir if output_dir else bp["raw_input_dir"] 
        print(f" --> Background for subrange(s) {g['subranges']}: catDict={cat_dict}, background_file={bkg_file!r}")

        # Merge and root steps
        if opt.run_background and not opt.stage_only:

            if output_dir and not opt.dry_run:
                os.makedirs(output_dir, exist_ok=True)

            pnn_step = "--do-predictions" if opt.run_pnn else ""
            maps_opts = f"""--process-map {bp['process_map']} --outfiles-map {bp['outfiles_map']} \
                --varDict {bp['varDict']} --cats --catDict {cat_dict}"""

            if is_real_batch(opt) and opt.root_step:
                # Phase 2: read the merged parquet back from where the merge
                # jobs landed (data_root_source_dir), not the raw input pool.
                cmd = f"""prepare_output_file.py --root --merge-data-only \
                    --input {data_root_source_dir} --output {data_root_source_dir} {maps_opts} {batch_opts(opt)}"""
            elif is_real_batch(opt):
                # Phase 1: merge jobs only -- see the two-step note at the top.
                extra_opts = f"--output {output_dir}" if output_dir else ""
                cmd = f"""prepare_output_file.py {pnn_step} --merge --merge-data-only \
                    --input {bp['raw_input_dir']} {extra_opts} {maps_opts} {batch_opts(opt)}"""
            else:
                extra_opts = f"--output {output_dir}" if output_dir else ""
                cmd = f"""prepare_output_file.py {pnn_step} --merge --merge-data-only --root \
                    --input {bp['raw_input_dir']} {extra_opts} {maps_opts} {batch_opts(opt)}"""

            ret = run(cmd, opt.dry_run, cwd=FINALFIT_DIR)
            if ret != 0 and not opt.dry_run:
                print(f"[ERROR] background postprocessing failed (exit code {ret}). Stopping.")
                sys.exit(ret)
            if is_real_batch(opt) and not opt.dry_run:
                if not opt.root_step:
                    print(f" --> Merge jobs submitted to {opt.batch} -- once they finish, re-run with "
                          f"--run-background --root-step.")
                else:
                    print(f" --> Root-conversion jobs submitted to {opt.batch} -- once they finish, re-run with "
                          f"--stage-only to pick up {bkg_file!r}.")
                # Jobs are only queued, not finished: nothing to stage yet, and
                # globbing now would copy a previous run's stale file into data/.
                continue
        else:
            reason = "--stage-only set" if opt.stage_only else "--run-background not set"
            print(f" --> {reason}: not (re-)running background postprocessing, "
                  f"just copying {bkg_file!r} if it already exists.")

        # mv files step -- reached from both branches above: after a completed
        # local run, and when we're only picking up what's already on disk.
        pattern = os.path.join(data_root_source_dir, "root", "Data", bkg_file)

        if opt.dry_run:
            print(f"[DRY RUN] would glob {pattern} and copy matches into {data_dir}/")
            continue

        matches = sorted(glob.glob(pattern))
        if not matches:
            print(f"[WARNING] No file matching {bkg_file!r} found at {pattern}. Run with "
                  f"--run-background, or check background_postprocessing.raw_input_dir/output_dir and "
                  f"background_file in the campaign yaml.")
            continue

        os.makedirs(data_dir, exist_ok=True)
        for src in matches:
            dst = os.path.join(data_dir, os.path.basename(src))
            print(f" --> cp {src} -> {dst}")
            shutil.copy(src, dst)


def build_dataset_selection(pp, common, subrange, dry_run):
    """Build a local dir of empty subdirs, one per dataset we want, to feed
    prepare_output_file.py's --folder-structure.

    Pulling from a shared pool that mixes every mass point and era (e.g.
    .../HDNA_output/2022/RSGrav/), prepare_output_file.py has no per-dataset
    filter: it processes every top-level subdir it's told to list.
    --folder-structure redirects that *listing* step only; --input/--output
    still control where data is read and written.
    """
    template = pp["dataset_name_template"]
    masses = [m.strip() for m in subrange["mass_points"].split(",")]
    dataset_names = [template.format(width=common["width"], mass=m, year=common["year"]) for m in masses]

    selection_dir = os.path.join(FINALFIT_DIR, ".dataset_selection", subrange["parquet_dir"])
    print(f" --> Dataset selection for '{subrange['name']}' ({len(dataset_names)} datasets):")
    for name in dataset_names:
        print(f"       {name}")

    if dry_run:
        print(f"[DRY RUN] would (re)create {selection_dir} with the above as empty subdirs")
        return selection_dir

    if os.path.isdir(selection_dir):
        shutil.rmtree(selection_dir)
    for name in dataset_names:
        os.makedirs(os.path.join(selection_dir, name))
    return selection_dir


def run_postprocessing(cfg, subrange, opt):
    common = cfg["common"]
    pp = cfg["signal_postprocessing"]
    dry_run = opt.dry_run

    # Where this subrange's output lands, whether produced fresh here or
    # already sitting there from a pre-split input (no raw_input_dir).
    output_dir = os.path.abspath(os.path.join(FINALFIT_DIR, common["tree_input_dir"], subrange["parquet_dir"]))

    if pp.get("raw_input_dir"):
        # Shared pool: select just this subrange's datasets via
        # --folder-structure and write to our own dir, rather than
        # processing (or writing into) the whole pool.
        selection_dir = build_dataset_selection(pp, common, subrange, dry_run)
        input_dir = pp["raw_input_dir"]
        extra_opts = f"--folder-structure {selection_dir} --output {output_dir}"
    else:
        input_dir = output_dir
        extra_opts = ""

    if not dry_run:
        os.makedirs(output_dir, exist_ok=True)

    cat_dict = subrange.get("catDict", pp.get("catDict"))
    if not cat_dict:
        print(f"[ERROR] subrange '{subrange['name']}' has no catDict -- set it on the subrange itself, "
              f"or on signal_postprocessing.catDict as a campaign-wide default.")
        sys.exit(1)

    if not opt.stage_only:
        syst_opt = "--syst" if pp.get("syst", True) else ""
        pnn_step = "--do-predictions" if opt.run_pnn else ""
        maps_opts = f"""--cats --catDict {cat_dict} \
            --process-map {pp['process_map']} \
            --outfiles-map {pp['outfiles_map']} \
            --varDict {pp['varDict']} {syst_opt}"""

        if is_real_batch(opt) and opt.root_step:
            # Phase 2: read the merged parquet back from output_dir (where
            # the merge jobs wrote it), not the original input_dir.
            cmd = f"""prepare_output_file.py --root --input {output_dir} {extra_opts} \
                {maps_opts} {batch_opts(opt)}"""
        elif is_real_batch(opt):
            # Phase 1: merge jobs only -- see the two-step note at the top.
            cmd = f"""prepare_output_file.py {pnn_step} --merge --input {input_dir} {extra_opts} \
                {maps_opts} {batch_opts(opt)}"""
        else:
            cmd = f"""prepare_output_file.py {pnn_step} --merge --root --input {input_dir} {extra_opts} \
                {maps_opts} {batch_opts(opt)}"""

        ret = run(cmd, dry_run, cwd=FINALFIT_DIR)
        if ret != 0 and not dry_run:
            print(f"[ERROR] postprocessing failed for subrange '{subrange['name']}' (exit code {ret}). Stopping.")
            sys.exit(ret)
        if is_real_batch(opt) and not dry_run:
            if not opt.root_step:
                print(f" --> Merge jobs submitted to {opt.batch} -- once they finish, re-run with --root-step "
                      f"for subrange '{subrange['name']}'.")
            else:
                print(f" --> Root-conversion jobs submitted to {opt.batch} -- once they finish, re-run with "
                      f"--stage-only to stage subrange '{subrange['name']}'.")
    else:
        print(f" --> --stage-only set: just staging existing output for subrange '{subrange['name']}'.")

    signal_dir = os.path.abspath(os.path.join(FINALFIT_DIR, common["tree_input_dir"], subrange["signal_dir"]))
    stage_signal_root_files(output_dir, signal_dir, dry_run)


def main():
    (opt, args) = get_options()
    if len(args) != 1:
        print("[ERROR] Usage: prepare_campaign_inputs.py campaign.yaml [options]")
        sys.exit(1)

    with open(args[0]) as f:
        cfg = yaml.safe_load(f)

    only = set(x.strip() for x in opt.only.split(",") if x.strip())
    if only: # only certain subranges
        subranges = [s for s in cfg["subranges"] if s["name"] in only]
    else: # all subranges
        subranges = cfg["subranges"]
    if not subranges:
        print(f"[ERROR] No subranges matched --only={opt.only!r}")
        sys.exit(1)

    print(f"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PREPARING INPUTS: {cfg['campaign']} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f" --> Subranges to run: {[s['name'] for s in subranges]}")

    # Once per distinct signal catDict.
    run_background_postprocessing(cfg, opt, subranges)

    for subrange in subranges:
        print(f"\n==================== Subrange '{subrange['name']}' ====================")
        run_postprocessing(cfg, subrange, opt)

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INPUTS READY -- now run run_campaign.py under cmsenv/setup.sh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


if __name__ == "__main__":
    main()
