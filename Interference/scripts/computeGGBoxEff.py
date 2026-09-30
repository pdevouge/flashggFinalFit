import os
import sys
import glob
import json
import yaml
from optparse import OptionParser
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.interpolate import UnivariateSpline

from commonObjects import *
from plottingTools import *

xsec_bkg = {
  'GGBox_MGG200To500': 7.351e3,
  'GGBox_MGG500To1000': 0.4631e3,
  'GGBox_MGG1000To2000': 0.05151e3,
  'GGBox_MGG2000': 0.003824e3,
  'GGBox_MGG80': 88.75e3 * 0.98
}


# ----------------------------------------------------------------------
def get_options():
    parser = OptionParser()
    parser.add_option('--verbose', dest='verbose', action='store_true')
    parser.add_option("--ext", dest='ext', default='', help="Extension of this run.")
    parser.add_option("--proc", dest='proc', default='', help="Signal process")
    parser.add_option("--cat", dest='cat', default='', help="RECO category")
    parser.add_option("--year", dest='year', default='2016', help="Year")
    parser.add_option('--catDict', dest='catDict', default=None, help='Category JSON. Required.')
    parser.add_option('--config', dest='config', default=None, help='YAML config file describing the reco_mass/gen_mass file groups. Required.')
    # Scan of mass points: either an explicit list, or min/max/n
    parser.add_option('--window', dest='window', default=None, type='float', help='Half-width of the mass window, as a fraction of m ')
    parser.add_option('--logSpace', action='store_true', help='Use log-spaced mass points instead of linear spacing')
    parser.add_option('--massList', dest='massList', default='', help='Comma-separated explicit list of mass points to scan (overrides --minMass/--maxMass/--nMassPoints if given)')
    parser.add_option('--minMass', dest='minMass', default=100., type='float', help='Minimum mass point of the scan')
    parser.add_option('--maxMass', dest='maxMass', default=3000., type='float', help='Maximum mass point of the scan')
    parser.add_option('--nMassPoints', dest='nMassPoints', default=40, type='int', help='Number of mass points in the scan')
    parser.add_option('--doPlots', dest='doPlots', action='store_true', help='Also draw the efficiency and the numerator/denominator.')

    opt, args = parser.parse_args()

    if opt.config is None:
        sys.exit("[ERROR] --config <yaml file> is required.")
    if opt.catDict is None:
        sys.exit("[ERROR] --catDict <category json> is required: the efficiency is "
                  "per analysis category, so the category selection has to be applied.")
    if opt.cat is None:
        sys.exit("[ERROR] --cat <category name> is required (a key of --catDict).")

    return opt, args


# ----------------------------------------------------------------------
# Operators as they appear in the catDict, mapped to what DataFrame.query expects.
cat_ops = {'=': '==', '==': '==', '!=': '!=', '>': '>', '>=': '>=', '<': '<', '<=': '<='}


def build_cat_cut(cat_dict_path, cat):
    """Turn one category's cat_filter into a DataFrame.query string."""
    try:
        with open(cat_dict_path) as f:
            cat_dict = json.load(f)
    except Exception as e:
        sys.exit(f"[ERROR] could not read --catDict '{cat_dict_path}': {e}")

    if cat not in cat_dict:
        sys.exit(f"[ERROR] category '{cat}' not in {cat_dict_path}. "
                  f"Available: {sorted(cat_dict)}")

    terms = []
    for entry in cat_dict[cat].get('cat_filter', []):
        if len(entry) != 3:
            sys.exit(f"[ERROR] malformed cat_filter entry in {cat_dict_path}: {entry}")
        var, op, val = entry
        if op not in cat_ops:
            sys.exit(f"[ERROR] unsupported operator '{op}' in {cat_dict_path} for '{var}'.")
        if isinstance(val, bool):
            val = 'True' if val else 'False'
        terms.append(f"{var} {cat_ops[op]} {val}")

    if not terms:
        sys.exit(f"[ERROR] category '{cat}' in {cat_dict_path} has an empty cat_filter.")
    # The columns the cut needs are the first element of each cat_filter entry
    return " and ".join(terms), [entry[0] for entry in cat_dict[cat]['cat_filter']]


# ----------------------------------------------------------------------
def load_config(config_path):
    with open(config_path, 'r') as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        sys.exit(f"[ERROR] could not parse '{config_path}' as a YAML mapping.")
    return cfg


def get_var_map(section):
    var_map = {'gen_mass': 'gen_dipho_mass', 'weighty': 'genWeight'}
    for entry in section:
        if isinstance(entry, dict) and 'var' in entry:
            var_map.update(entry['var'] or {})
            break  # only the first 'var' entry is used
    return var_map


def get_columns(section):
    """Columns to read from the parquets, declared in the config beside 'var'.
    Returns None if the config does not say, meaning read every column."""
    for entry in section:
        if isinstance(entry, dict) and 'columns' in entry:
            return list(entry['columns'] or [])
    return None

# ----------------------------------------------------------------------
def resolve_file_list(files_entry):
    """Expand files entry with glob.glob and return a flat, 
    deduplicated, ordered list of files."""
    if isinstance(files_entry, str):
        files_entry = [files_entry]

    resolved = []
    for pattern in files_entry:
        matches = sorted(glob.glob(pattern))
        if not matches:
            print(f"     [WARN] no files matched pattern: '{pattern}'")
            continue
        for m in matches:
            if m not in resolved:
                resolved.append(m)
    return resolved

def resolve_xsec_key(paths):
    """Find which xsec_bkg key is contained in the file path(s) for this group."""
    xsec = []
    labels = []
    for p in paths:
        for key in xsec_bkg:
            if key in p:
                xsec.append(xsec_bkg[key])
                labels.append(key)
    if len(xsec) != len(paths):
        sys.exit(f"[ERROR] Number of xsec is different from number of paths")
    return xsec, labels


def read_sum_genw(path):
    """Read one file's sum of gen weights before selection from its metadata."""
    metadata = pq.ParquetFile(path).schema_arrow.metadata or {}
    if b'sum_genw_presel' not in metadata:
        sys.exit(f"[ERROR] no 'sum_genw_presel' in the metadata of\n"
                 f"        {path}\n"
                 f"        Run tools/backfill_sum_genw.py on the directory holding it.")
    return float(metadata[b'sum_genw_presel'])


def normalize_weight(df, weight_col, sum_genw):
    """Divide each sample's weights by its own sum of gen weights before selection."""
    df[weight_col] = df[weight_col] / df['label'].map(sum_genw)
    return df


def load_group(gtype, group, id1_col, id2_col, weight_col, normalize=False, verbose=False, extra_cut=None, columns=None):
    """Load one file-group, concatenate its files, apply the ggbox 
    parton-level selection (Generator_id1 == 21 and Generator_id2 == 21),
    then apply the group's own 'cut' if any. """
    if 'files' not in group:
        sys.exit(f"[ERROR] file group is missing the 'files' key: {group}")

    paths = resolve_file_list(group['files'])
    if not paths:
        sys.exit(f"[ERROR] no parquet files matched for group: {group}")

    xsec, labels = resolve_xsec_key(paths)

    # The category selection is ANDed onto whatever the group already asks for.
    cuts = [c for c in (group.get('cut'), extra_cut) if c]
    cut = " and ".join(f"({c})" for c in cuts) if cuts else None

    print(f"     found {len(paths)} file(s)" + (f", cut='{cut}'" if cut else ""))

    dfs = []
    sum_genw = {}
    for p, x, label in zip(paths, xsec, labels):
        if verbose: print(f"       - {p}")
        df = pd.read_parquet(p, columns=columns)
        if normalize:
            sum_genw[label] = sum_genw.get(label, 0.) + read_sum_genw(p)
        df['label'] = label
        df[weight_col] = df[weight_col] * x
        for col in (id1_col, id2_col):
            if col not in df.columns:
                sys.exit(f"[ERROR] column '{col}' not found in {p}. "
                          f"Available columns: {list(df.columns)}")
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True, sort=False)
    if normalize:
        df = normalize_weight(df, weight_col, sum_genw)
    df['first1'] = (df['flags1'] & 2**12 == 0) if gtype == 'gen' else (df['gen_lead_statusFlags'] & 2**12 == 0)
    df['first2'] = (df['flags2'] & 2**12 == 0) if gtype == 'gen' else (df['gen_sublead_statusFlags'] & 2**12 == 0)

    sel = (df[id1_col] == 21) & (df[id2_col] == 21)
    df = df.loc[sel].copy()

    if cut:
        try:
            df = df.query(cut)
        except Exception as e:
            sys.exit(f"[ERROR] failed to apply cut '{cut}' on group {group['files']}: {e}")

    return df


def load_ggbox_sample_from_config(gtype, section, id1_col, id2_col, weight_col, verbose=False, extra_cut=None, columns=None):
    """Load and apply cut on each file group, then concatenate them all together."""
    groups = [entry for entry in section if isinstance(entry, dict) and 'files' in entry]
    if not groups:
        sys.exit("[ERROR] expected at least one file group (with a 'files' key) in the config.")

    normalize = next(entry['normalize'] for entry in section if 'normalize' in entry)
    print(f"   * normalize: {normalize}")

    dfs = []
    for i, group in enumerate(groups):
        print(f"   * group {i + 1}/{len(groups)}")
        dfs.append(load_group(gtype, group, id1_col, id2_col, weight_col, normalize=normalize,
                              verbose=verbose, extra_cut=extra_cut, columns=columns))
 
    return pd.concat(dfs, ignore_index=True, sort=False)


# ----------------------------------------------------------------------
def weighted_sum_and_err(weights):
    """Sum of weights and its statistical uncertainty (sqrt(sum w^2))"""
    w = np.asarray(weights, dtype=float)
    s = w.sum()
    e = np.sqrt(np.sum(w ** 2))
    return s, e


# ----------------------------------------------------------------------
def compute_efficiency(df_reco, df_gen, reco_var, gen_var, masses, window=None):
    """
    Compute eff_cat(m) = sum_w_reco,cat[window] / sum_w_gen[window] for every m in `masses` 
    """
    for var in gen_var:
        if gen_var[var] not in df_gen.columns:
            sys.exit(f"[ERROR] column '{gen_var[var]}' (gen '{var}') not found in "
                        f"gen sample. Available columns: {list(df_gen.columns)}")
    for var in reco_var:
        if reco_var[var] not in df_reco.columns:
            sys.exit(f"[ERROR] column '{reco_var[var]}' (reco '{var}') not found in "
                        f"reco sample. Available columns: {list(df_reco.columns)}")

    gen_mass_col, gen_weight_col = gen_var['gen_mass'], gen_var['weight']
    reco_mass_col, reco_weight_col = reco_var['gen_mass'], reco_var['weight']

    rows = []
    if window:
        bins = [((1.0 - window) * m, (1.0 + window) * m) for m in masses]
    else:
        bins = zip(masses[:-1], masses[1:])

    for lo, hi in bins:
        gen_mask = (df_gen[gen_mass_col] > lo) & (df_gen[gen_mass_col] < hi)
        n_gen = len(df_gen.loc[gen_mask, gen_weight_col])
        den, den_err = weighted_sum_and_err(df_gen.loc[gen_mask, gen_weight_col])

        reco_mask = (df_reco[reco_mass_col] > lo) & (df_reco[reco_mass_col] < hi)
        n_reco = len(df_reco.loc[reco_mask, reco_weight_col])
        num, num_err = weighted_sum_and_err(df_reco.loc[reco_mask, reco_weight_col])

        if den > 0:
            eff = num / den
            rel_num = (num_err / num) if num > 0 else 0.
            rel_den = den_err / den
            eff_err = eff * np.sqrt(rel_num ** 2 + rel_den ** 2)
        else:
            eff, eff_err = np.nan, np.nan

        rows.append({
            'mLow': lo, 'mNom': (lo+hi)//2, 'mHigh': hi,
            'n_reco': n_reco, 'n_gen': n_gen,
            'num': num, 'num_err': num_err,
            'den': den, 'den_err': den_err,
            'eff': eff, 'eff_err': eff_err,
        })

    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
def add_spline_fit(eff_df):
    """
    Smooth the efficiency: fit a cubic spline through the (mNom, eff) points,
    weighted by 1/eff_err so noisy points pull the curve less.

    Points outside the fitted mass range are clamped to the boundary value
    (ext=3) rather than extrapolated, since the spline can swing to negative
    efficiencies outside the range the data constrained.
    """
    mask = eff_df['eff'].notna() & (eff_df['eff_err'] > 0)

    x = eff_df.loc[mask, 'mNom'].to_numpy(dtype=float)
    y = eff_df.loc[mask, 'eff'].to_numpy(dtype=float)
    w = 1.0 / eff_df.loc[mask, 'eff_err'].to_numpy(dtype=float)

    spline = UnivariateSpline(x, y, w=w, ext=3)
    # An efficiency cannot be negative: the spline can undershoot if eff is close to 0. Clip to avoid issues.
    eff_df['eff_fit'] = np.clip(spline(eff_df['mNom'].to_numpy(dtype=float)), 0., None)
    return eff_df


# ----------------------------------------------------------------------
def main():
    opt, args = get_options()

    print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GGBox efficiency ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")

    print(f" --> Reading config from {opt.config}")
    cfg = load_config(opt.config)

    reco_var = get_var_map(cfg["reco"])
    gen_var = get_var_map(cfg["gen"])

    print(f" --> Loading reco sample(s) from config['reco']")
    cat_cut, cat_cols = build_cat_cut(opt.catDict, opt.cat)
    print(f"   * category '{opt.cat}' -> {cat_cut}")
    reco_cols = get_columns(cfg["reco"])
    if reco_cols is not None: reco_cols = sorted(set(reco_cols) | set(cat_cols))
    df_reco = load_ggbox_sample_from_config("reco", cfg["reco"], "generator_id1", "generator_id2", reco_var['weight'], verbose=opt.verbose, extra_cut=cat_cut, columns=reco_cols)
    print(f"     {len(df_reco)} events pass Generator_id1==21 and Generator_id2==21"
          f" (+ per-group cuts)")

    print(f" --> Loading gen sample(s) from config['gen']")
    df_gen = load_ggbox_sample_from_config("gen", cfg["gen"], "Generator_id1", "Generator_id2", gen_var['weight'], verbose=opt.verbose, columns=get_columns(cfg["gen"]))
    print(f"     {len(df_gen)} events pass Generator_id1==21 and Generator_id2==21"
          f" (+ per-group cuts)")

    # Build the mass scan
    if opt.massList:
        masses = np.array([float(x) for x in opt.massList.split(',')])
    elif opt.logSpace:
        masses = np.logspace(np.log10(opt.minMass), np.log10(opt.maxMass), opt.nMassPoints)
    else:
        masses = np.linspace(opt.minMass, opt.maxMass, opt.nMassPoints)

    print(f" --> Scanning {len(masses)} mass points between "
          f"{masses.min():.1f} and {masses.max():.1f} GeV ")
    if opt.window: print(f"(window = +/-{opt.window*100:.0f}%)")

    eff_df = compute_efficiency(df_reco, df_gen, reco_var, gen_var, masses, opt.window)

    print(f" --> Smoothing the efficiency with a spline")
    eff_df = add_spline_fit(eff_df)

    outDir = "%s/results/outdir_%s/computeGGBoxEff/output"%(iwd__,opt.ext)
    if not os.path.isdir(outDir): os.makedirs(outDir)
    outCsv = "%s/ggbox_eff_%s_%s_%s_%s.csv"%(outDir,opt.ext,opt.proc,opt.year,opt.cat)

    eff_df.to_csv(outCsv, index=False)
    print(f" --> Efficiency table saved to {outCsv}")

    if opt.doPlots:
        plotDir = "%s/results/outdir_%s/computeGGBoxEff/plots"%(iwd__,opt.ext)
        if not os.path.isdir(plotDir): os.makedirs(plotDir)
        # tag by year and category, so per-category runs do not overwrite each other
        print(f" --> Plotting into {plotDir}")
        plotGGBoxEfficiency(outCsv, plotDir, "_%s_%s_%s"%(opt.proc,opt.year,opt.cat))
        plotGGBoxNumDen(outCsv, plotDir, "_%s_%s_%s"%(opt.proc,opt.year,opt.cat))

    print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")


if __name__ == '__main__':
    main()