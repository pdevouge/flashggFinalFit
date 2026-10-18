# Final Fits (higgsdnafinalfit)

This is the branch for using final fits with the output of HiggsDNA.

Welcome to the new Final Fits package. Here lies a a series of scripts which are used to run the final stages of CMS Hgg analyses: signal modelling, background modelling, datacard creation, final statistical interpretation and final result plots.

You can follow the latest Final Fits tutorial [here](https://gitlab.cern.ch/jspah/higgsdna_finalfits_tutorial_24/-/tree/master). Slides from an older flashgg-based tutorial series can be found [here](https://indico.cern.ch/event/963619/contributions/4112177/attachments/2151275/3627204/finalfits_tutorial_201126.pdf)

## Download and setup instructions

```
export SCRAM_ARCH=el9_amd64_gcc12
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv

COMBINE_TAG=high_mass_analysis
COMBINEHARVESTER_TAG=high_mass_analysis
FINALFIT_TAG=high_mass_finalfit

# Install Combine with the latest EL9 compatible branch
git clone https://github.com/pdevouge/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit && git fetch origin ${COMBINE_TAG} && git checkout ${COMBINE_TAG}

# Install CombineTools in CombineHarvester
cd ${CMSSW_BASE}/src
git clone -b $COMBINEHARVESTER_TAG https://github.com/pdevouge/CombineHarvester.git

# Compile libraries
cd ${CMSSW_BASE}/src
cmsenv
scram b clean
scram b -j 8

# Install Final Fit package
git clone -b $FINALFIT_TAG https://github.com/pdevouge/flashggFinalFit.git
cd flashggFinalFit/
source setup.sh
```

## Contents
The Final Fits package contains several subfolders, each covering one step of the pipeline:

* `Trees2WS` -- converts HiggsDNA/flashgg ROOT trees into RooWorkspaces (signal and data)
* `Signal` -- fits the signal model and packages it per category/year
* `Interference` -- builds the signal-background-interference (SBI) model for the ggbox continuum (optional)
* `Background` -- fits the data-driven background model (discrete profiling / f-test)
* `Datacard` -- computes signal yields/systematics and writes the Combine text datacard
* `Combine` -- runs `text2workspace.py` and `combineTool.py` to extract limits, and holds plotting/checks scripts
* `Plots` -- standalone scripts to produce result plots (impacts, S+B, summary, etc.)

`Trees2WS`, `Signal`, `Interference` and `Background` each have their own `config/` (per-run `.py` configs) and `results/` (fit outputs) subdirectories, kept separate from the scripts themselves. `Datacard` and `Combine` only have `results/`, since their inputs are passed as command-line options rather than a config file.

Objects and tools common to all subfolders (luminosity map, batch submission helpers, etc.) live in `tools/`. If your input workspaces differ from the flashgg/HiggsDNA output structure, you may need to change the options there.

Each of these folders has its own `README.md` describing its scripts and options in more detail.

On top of these per-step folders, this fork adds an orchestration layer that drives the whole pipeline (`Trees2WS` -> `Signal`/`Interference` -> `Background` -> `Datacard` -> `Combine`) from a single small config, described below.

## How to run finalfit
### One mass subrange/config: `runFinalfit.py`

`runFinalfit.py` runs the entire pipeline above for one analysis (one mass subrange, one category set, one or several years merged via `--mergeYears`) from a single YAML config:

```
python3 runFinalfit.py <config.yaml> [--doSystematics] [--skipIntf] [--onlyYear <era>] \
    [--runOnly trees,signal,interference,background,datacard,combine,text2ws,limits,collect] \
    [--lumiscale <factor>] [--dry-run]
```

It writes the per-step `config/config_<extension>.py` files, `cd`s into each section in turn, and finally collects everything (plots, packaged signal workspace, datacard, Combine workspace, limits json+plot) into `Results/<extension>/`. Use `--dry-run` to print every command without running it, and `--runOnly` to re-run just a subset of steps (e.g. after fixing one stage).

The config schema (`common`, `tree`, `signal`, `background`, optional `limits`/`combine` blocks) is documented in the option help and inline comments at the top of `runFinalfit.py` itself.

### Multiple subranges/years: campaigns

As our high-mass analyse span several mass subranges (and several years merged together), the entry point is a **campaign** YAML under `campaigns/` (see e.g. `campaigns/rsg_multicats_2022_full.yaml`) rather than calling `runFinalfit.py` directly. A campaign describes, once, the raw HiggsDNA postprocessing inputs and the common settings, plus a list of `subranges:` (mass binning, categorization `catDict`, mass points, etc.).

Running a campaign is two steps:

```
# 1) Stage the postprocessed signal/background inputs (needs a HiggsDNA environment)
python3 prepare_campaign_inputs.py campaigns/rsg_multicats_2022_full.yaml [--run-background] [--dry-run]

# 2) Run runFinalfit.py for every subrange (needs cmsenv / source setup.sh)
python3 run_campaign.py campaigns/rsg_multicats_2022_full.yaml [--only 500-1000] [--skipIntf] [--dry-run]
```

`run_campaign.py` writes one `generated_configs/<campaign>_<subrange>.yaml` per subrange and calls `runFinalfit.py` on each; results land in `Results/<campaign>_<subrange>/`. Pass `--batch condor` to submit one HTCondor job per subrange instead of running locally, or use `--jobs N` to run several subranges locally at once.

### Where results end up

Whichever way you ran it, each analysis (one subrange, one extension) ends up collected under `Results/<extension>/`:
```
Results/<extension>/
  signal/        # packaged signal workspaces + fit plots
  background/    # multipdf background workspaces + f-test plots
  interference/  # SBI model (unless --skipIntf)
  Datacard.txt
  combine/       # Combine workspace, Models/, Limits/ (limit json + plot)
```

