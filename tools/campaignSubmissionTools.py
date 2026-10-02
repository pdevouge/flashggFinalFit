# Batch submission helpers for run_campaign.py.
#
# --batch condor submits ONE HTCondor job PER SUBRANGE, each re-running
# run_campaign.py --only <subrange> --batch local on its own worker node.
# Every stage for that subrange still runs in order inside its own single
# process (easy to track down with one log per subrange).
# The trade-off is that the fits WITHIN one subrange's job stay serial, with
# one exception: the limits stage honours --limitJobs, which is forwarded to
# the inner run_campaign.py and mirrored into RequestCpus so combineTool's
# process pool actually has cores to run on.
import os
import subprocess

DEFAULT_SCRAM_ARCH = "el9_amd64_gcc12"

# A subrange's chain (signal+background+datacard+combine+limits) runs for
# hours, so the short flavours would be killed mid-run.
# espresso 20min / microcentury 1h / longlunch 2h / workday 8h / tomorrow 1d.
DEFAULT_FLAVOUR = "tomorrow"


def build_inner_command(campaign_yaml, opt, subrange_name):
    """The run_campaign.py invocation that one worker node will run, scoped
    to exactly one subrange.

    --batch is deliberately forced back to "local": passing the batch mode
    through would have the job submit another job, forever.
    """
    cmd = f"python3 -u run_campaign.py {campaign_yaml} --batch local --only {subrange_name}"
    if opt.stages:
        cmd += f" --stages {opt.stages}"
    if opt.do_syst:
        cmd += " --doSystematics"
    if opt.skip_intf:
        cmd += " --skipIntf"
    if opt.limit_jobs != 1:
        cmd += f" --limitJobs {opt.limit_jobs}"
    if opt.impact_jobs != 1:
        cmd += f" --impactJobs {opt.impact_jobs}"
    if opt.lumiscale:
        cmd += f" --lumiscale {opt.lumiscale}"
    return cmd


def build_wrapper_script(inner_cmd, finalfit_dir, cmssw_base, scram_arch=None):
    """Bash wrapper that sets CMSSW up on the worker and runs one subrange."""
    return f"""#!/bin/bash
            set -e
            export SCRAM_ARCH={scram_arch or os.environ.get("SCRAM_ARCH", DEFAULT_SCRAM_ARCH)}
            source /cvmfs/cms.cern.ch/cmsset_default.sh
            cd {cmssw_base}/src
            eval `scramv1 runtime -sh`
            export PYTHONPATH=$PYTHONPATH:{cmssw_base}/src/flashggFinalFit/tools
            cd {finalfit_dir}
            {inner_cmd}
            """


def eos_redirector(path):
    """root:// redirector for an /eos path"""
    p = os.path.realpath(path)
    if p.startswith("/eos/cms"):
        return "root://eoscms.cern.ch/" + p
    if p.startswith("/eos/user") or p.startswith("/eos/home-"):
        return "root://eosuser.cern.ch/" + p
    return None


def build_condor_submit_file(tag, job_dir, flavour, request_cpus=1):
    """
    HTCondor submit description for one subrange's job.

    RequestCpus follows max(--limitJobs, --impactJobs): 1 CPU per mass point / per
    nuisance-parameter fit running in parallel. The limits and impacts stages run
    one after another within the same job, never at once, so this is a max, not a sum.
    """
    destination = eos_redirector(job_dir)
    return f"""universe                = vanilla
            initialdir              = {job_dir}
            executable              = {job_dir}/{tag}.sh
            output                  = {tag}.$(ClusterId).out
            error                   = {tag}.$(ClusterId).err
            log                     = {job_dir}/{tag}.$(ClusterId).log
            output_destination      = {destination}
            MY.XRDCP_CREATE_DIR     = True
            getenv                  = False
            RequestCpus             = {request_cpus}
            +JobFlavour             = "{flavour}"
            MY.SendCredential       = true
            on_exit_hold            = (ExitBySignal == True) || (ExitCode != 0)
            periodic_release        = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)
            requirements            = Machine =!= LastRemoteHost
            queue
            """


def submit_campaign(cfg, campaign_yaml, opt, finalfit_dir):
    """
    Write one wrapper + submit file per subrange and submit each as its
    own HTCondor job, so every subrange gets a dedicated worker.
    """
    campaign = cfg["campaign"]
    cmssw_base = os.path.dirname(os.path.dirname(finalfit_dir))
    job_dir = os.path.abspath(os.path.join(finalfit_dir, opt.condor_dir))

    only = set(x.strip() for x in opt.only.split(",") if x.strip()) if opt.only else None
    subranges = [s["name"] for s in cfg["subranges"] if not only or s["name"] in only]
    if not subranges:
        print(f"[ERROR] No subranges matched --only={opt.only!r}")
        return 1

    if not opt.dry_run:
        os.makedirs(job_dir, exist_ok=True)

    has_failed = 0
    for subrange_name in subranges:
        tag = f"run_campaign_{campaign}_{subrange_name.replace('-', '_')}"
        inner = build_inner_command(campaign_yaml, opt, subrange_name)
        sh = build_wrapper_script(inner, finalfit_dir, cmssw_base)
        sub = build_condor_submit_file(tag, job_dir, opt.flavour, max(1, opt.limit_jobs, opt.impact_jobs))

        if opt.dry_run:
            print(f"[DRY RUN] would write {job_dir}/{tag}.sh and {tag}.sub, then condor_submit it")
            print("\n----- wrapper -----\n" + sh + "----- submit file -----\n" + sub)
            continue

        sh_path = os.path.join(job_dir, tag + ".sh")
        sub_path = os.path.join(job_dir, tag + ".sub")
        with open(sh_path, "w") as f:
            f.write(sh)
        os.chmod(sh_path, 0o755)
        with open(sub_path, "w") as f:
            f.write(sub)

        print(f" --> Wrote {sh_path}")
        print(f" --> Wrote {sub_path}")

        # -spool matches HiggsDNA's remote/htcondor.py submit_jobs(): needed
        # alongside output_destination for the schedd to accept /eos submissions.
        ret = subprocess.call(f"condor_submit -spool {sub_path}", shell=True, cwd=job_dir)
        if ret != 0:
            print(f"[ERROR] condor_submit failed for subrange '{subrange_name}' (exit code {ret}).")
            if not has_failed:
                has_failed = ret
        else:
            print(f" --> Submitted subrange '{subrange_name}'.")

    if opt.dry_run or has_failed:
        return has_failed

    print(f"\n --> {len(subranges)} subrange job(s) submitted, one per worker. Watch them with:")
    print("       condor_q")
    print(f"       tail -f {job_dir}/run_campaign_{campaign}_*.*.out")
    print(f" --> Per-subrange progress also appears under {finalfit_dir}/logs/run_<ext>.log.")
    return 0
