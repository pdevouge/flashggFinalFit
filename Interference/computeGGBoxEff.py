import glob
import sys
from optparse import OptionParser

import numpy as np
import pandas as pd
import yaml

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

    parser.add_option('--config', dest='config', default=None,
                       help='YAML config file describing the reco_mass/gen_mass '
                            'file groups (and their per-group mass cuts). Required.')

    parser.add_option('--window', dest='window', default=0.10, type='float',
                       help='Half-width of the mass window, as a fraction of m '
                            '(default 0.10 -> [0.90*m, 1.10*m])')

    # Scan of mass points: either an explicit list, or min/max/n
    parser.add_option('--logSpace', action='store_true',
                     help='Use log-spaced mass points instead of linear spacing')
    parser.add_option('--massList', dest='massList', default='',
                       help='Comma-separated explicit list of mass points to scan '
                            '(overrides --minMass/--maxMass/--nMassPoints if given)')
    parser.add_option('--minMass', dest='minMass', default=100., type='float',
                       help='Minimum mass point of the scan')
    parser.add_option('--maxMass', dest='maxMass', default=3000., type='float',
                       help='Maximum mass point of the scan')
    parser.add_option('--nMassPoints', dest='nMassPoints', default=40, type='int',
                       help='Number of mass points in the scan')

    parser.add_option('--outCsv', dest='outCsv', default='ggbox_efficiency.csv',
                       help='Output CSV file')

    opt, args = parser.parse_args()

    if opt.config is None:
        sys.exit("[ERROR] --config <yaml file> is required (see module docstring "
                  "for the expected format).")

    return opt, args


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

# ----------------------------------------------------------------------
def resolve_file_list(files_entry):
    """Expand files entry entry with glob.glob and return a flat, 
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


def load_group(gtype, group, id1_col, id2_col, weight_col, verbose=False):
    """Load one file-group, concatenate its files, apply the ggbox 
    parton-level selection (Generator_id1 == 21 and Generator_id2 == 21),
    then apply the group's own 'cut' if any. """
    if 'files' not in group:
        sys.exit(f"[ERROR] file group is missing the 'files' key: {group}")

    paths = resolve_file_list(group['files'])
    if not paths:
        sys.exit(f"[ERROR] no parquet files matched for group: {group}")

    xsec, labels = resolve_xsec_key(paths)

    print(f"     found {len(paths)} file(s)"
          + (f", cut='{group['cut']}'" if group.get('cut') else ""))          

    dfs = []
    for p, x, label in zip(paths, xsec, labels):
        if verbose: print(f"       - {p}")
        df = pd.read_parquet(p)
        df['label'] = label
        df[weight_col] = df[weight_col] * x
        for col in (id1_col, id2_col):
            if col not in df.columns:
                sys.exit(f"[ERROR] column '{col}' not found in {p}. "
                          f"Available columns: {list(df.columns)}")
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True, sort=False)
    df['first1'] = (df['flags1'] & 2**12 == 0) if gtype == 'gen' else (df['gen_lead_statusFlags'] & 2**12 == 0)
    df['first2'] = (df['flags2'] & 2**12 == 0) if gtype == 'gen' else (df['gen_sublead_statusFlags'] & 2**12 == 0)

    sel = (df[id1_col] == 21) & (df[id2_col] == 21)
    df = df.loc[sel].copy()

    cut = group.get('cut')
    if cut:
        try:
            df = df.query(cut)
        except Exception as e:
            sys.exit(f"[ERROR] failed to apply cut '{cut}' on group {group['files']}: {e}")

    return df


def load_ggbox_sample_from_config(gtype, section, id1_col, id2_col, weight_col, verbose=False):
    """Load and apply cut on each file group, then concatenate them all together."""
    groups = [entry for entry in section if isinstance(entry, dict) and 'files' in entry]
    if not groups:
        sys.exit("[ERROR] expected at least one file group (with a 'files' key) in the config.")
 
    dfs = []
    for i, group in enumerate(groups):
        print(f"   * group {i + 1}/{len(groups)}")
        dfs.append(load_group(gtype, group, id1_col, id2_col, weight_col, verbose=verbose))
 
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
def main():
    opt, args = get_options()

    print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GGBox efficiency ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")

    print(f" --> Reading config from {opt.config}")
    cfg = load_config(opt.config)

    reco_var = get_var_map(cfg["reco"])
    gen_var = get_var_map(cfg["gen"])

    print(f" --> Loading reco sample(s) from config['reco']")
    df_reco = load_ggbox_sample_from_config("reco", cfg["reco"], "generator_id1", "generator_id2", reco_var['weight'], verbose=opt.verbose)
    print(f"     {len(df_reco)} events pass Generator_id1==21 and Generator_id2==21"
          f" (+ per-group cuts)")

    print(f" --> Loading gen sample(s) from config['gen']")
    df_gen = load_ggbox_sample_from_config("gen", cfg["gen"], "Generator_id1", "Generator_id2", gen_var['weight'], verbose=opt.verbose)
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
          f"{masses.min():.1f} and {masses.max():.1f} GeV "
          f"(window = +/-{opt.window*100:.0f}%)")

    eff_df = compute_efficiency(df_reco, df_gen, reco_var, gen_var, masses, opt.window)

    eff_df.to_csv(opt.outCsv, index=False)
    print(f" --> Efficiency table saved to {opt.outCsv}")

    print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")


if __name__ == '__main__':
    main()