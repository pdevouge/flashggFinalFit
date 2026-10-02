from optparse import OptionParser
import ROOT
from CombineHarvester.CombineTools.plotting import *
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)


def get_options():
    parser = OptionParser()
    # Repeatable: --input can be passed multiple times, one per mass range/analysis
    parser.add_option('--input', dest='inputJson', action='append', default=[],
                       help="Limits.json file. Can be given multiple times, "
                            "e.g. --input low.json --input high.json")
    parser.add_option('--label', dest='labels', action='append', default=[],
                       help="Legend label for the corresponding --input, in the same order, "
                            "e.g. --label '350-650 GeV' --label '500-1000 GeV'")
    parser.add_option('--color', dest='colors', action='append', default=[],
                       help="ROOT color code for the median/observed line of the corresponding --input")
    parser.add_option('--range', dest='ranges', action='append', default=[],
                       help="Restrict the corresponding --input to this mass window before drawing, given as 'xmin,xmax'"
                            "Pass an empty string '' to leave a given input untrimmed.")
    parser.add_option('--title', dest='title', default='Limit', help="Plot title")
    parser.add_option('--unblinded', dest='unblinded', action='store_true', help="Unblind limit plot")
    parser.add_option('--logy', dest='logy', action='store_true', help="Draw the y-axis with a log scale")
    return parser.parse_args()


def trim_graph_to_range(g, xmin, xmax):
    """Remove any points from graph g whose x value falls outside [xmin, xmax].
    Works in-place; safe for TGraph and TGraphAsymmErrors (RemovePoint keeps
    the associated asymmetric errors consistent)."""
    for i in reversed(range(g.GetN())):
        x = g.GetX()[i]
        if x < xmin or x > xmax:
            g.RemovePoint(i)


(opt, args) = get_options()

if not opt.inputJson:
    opt.inputJson = ['limits_default.json']

# Default label/color fallbacks if not fully specified
default_colors = [ROOT.kRed, ROOT.kBlue + 1, ROOT.kMagenta + 1, ROOT.kOrange + 1, ROOT.kCyan + 2, ROOT.kBlack]
labels = opt.labels + [''] * (len(opt.inputJson) - len(opt.labels))
colors = [int(c) for c in opt.colors] + default_colors[len(opt.colors):]

# Parse --range strings ("xmin,xmax") into (xmin, xmax) tuples, or None if not given/empty
raw_ranges = opt.ranges + [''] * (len(opt.inputJson) - len(opt.ranges))
ranges = []
for r in raw_ranges:
    if r:
        xmin_str, xmax_str = r.split(',')
        ranges.append((float(xmin_str), float(xmax_str)))
    else:
        ranges.append(None)

# Get limit TGraphs as a dictionary
draw = ['obs', 'exp0', 'exp1', 'exp2'] if opt.unblinded else ['exp0', 'exp1', 'exp2']
all_graphs = [StandardLimitsFromJSONFile(inp, draw=draw) for inp in opt.inputJson]

# Trim each dataset to its requested mass window
for graphs, rng in zip(all_graphs, ranges):
    if rng is None:
        continue
    xmin, xmax = rng
    for g in graphs.values():
        trim_graph_to_range(g, xmin, xmax)

# Style and pads
ModTDRStyle()
canv = ROOT.TCanvas('limit', 'limit')
pads = OnePad()

# Build a combined x-axis spanning the union of all mass ranges
x_min, x_max = None, None
for graphs in all_graphs:
    for g in graphs.values():
        n = g.GetN()
        if n == 0:
            continue
        xs = [g.GetX()[i] for i in range(n)]
        gmin, gmax = min(xs), max(xs)
        x_min = gmin if x_min is None else min(x_min, gmin)
        x_max = gmax if x_max is None else max(x_max, gmax)

axis = CreateAxisHist(list(all_graphs[0].values())[0])
axis.GetXaxis().SetLimits(x_min, x_max)
axis.GetXaxis().SetTitle('m_{X} (GeV)')
axis.GetYaxis().SetTitle('95% CL limit on #it{#sigma#times#bf{B}} [pb]')
pads[0].cd()
if opt.logy:
    pads[0].SetLogy(True)
axis.Draw('axis')

# Create a legend in the top right
legend = PositionedLegend(0.3, 0.2 + 0.035 * len(all_graphs), 3, 0.015)

# Set the standard green and yellow colors and draw
# Transparency lets overlapping green/yellow bands from different mass
band_alpha = 0.85

for i, (graphs, color) in enumerate(zip(all_graphs, colors)):
    trans_green = CreateTransparentColor(ROOT.kGreen, band_alpha)
    trans_yellow = CreateTransparentColor(ROOT.kYellow, band_alpha)

    overwrite_style = {
        'exp0': {'LineColor': color, 'LineWidth': 2, 'LineStyle': 1 if i == 0 else 2},
        'exp1': {'FillColor': trans_green},
        'exp2': {'FillColor': trans_yellow},
    }
    if 'obs' in graphs:
        overwrite_style['obs'] = {'LineColor': color, 'MarkerColor': color, 'LineWidth': 2}

    StyleLimitBand(graphs, overwrite_style_dict=overwrite_style)

for graphs in all_graphs:
    band_keys = [k for k in ['exp2', 'exp1'] if k in graphs]
    DrawLimitBand(pads[0], graphs, draw=band_keys, draw_legend=[])

for i, (graphs, label) in enumerate(zip(all_graphs, labels)):
    line_keys = [k for k in ['exp0', 'obs'] if k in graphs]

    legend_overwrite = {
        'exp0': {'Label': f'{label} expected' if label else 'Expected'},
    }
    if 'obs' in graphs:
        legend_overwrite['obs'] = {'Label': f'{label} observed' if label else 'Observed'}

    # Only add the +-1/+-2 sigma band rows to the legend once, for the first input (avoid repeat for every mass range)
    legend_keys = line_keys + (['exp1', 'exp2'] if i == 0 else [])
    if i == 0:
        legend_overwrite['exp1'] = {'Label': '#pm1#sigma Expected'}
        legend_overwrite['exp2'] = {'Label': '#pm2#sigma Expected'}

    DrawLimitBand(pads[0], graphs, draw=line_keys, draw_legend=legend_keys,
                  legend=legend, legend_overwrite=legend_overwrite)

legend.Draw()

# Re-draw the frame and tick marks
pads[0].RedrawAxis()
pads[0].GetFrame().Draw()

# Adjust the y-axis range across the full combined mass range
if opt.logy:
    FixBothRanges(pads[0], GetPadYMin(pads[0]), 0.1, GetPadYMax(pads[0]), 0.25)
else:
    FixBothRanges(pads[0], 0, 0, GetPadYMax(pads[0]), 0.25)

# Standard CMS logo
if opt.title:
    main_title, sub_title = opt.title.split(',')
else:
    main_title, sub_title = '', ''
DrawCMSLogo(pads[0], 'CMS', 'Internal', 0, 0.13, 0.035, 1.2, '', 0.8)
DrawCMSLogo(pads[0], main_title, sub_title, 11, 0.2, 0.035, 1.2, '', 0.8)

# Luminosity text
lumi = ROOT.TLatex()
lumi.SetNDC()
lumi.SetTextFont(42)
lumi.SetTextSize(0.04)
lumi.SetTextAlign(31)
lumi.DrawLatex(0.95, 0.96, "138 fb^{-1} (13.6 TeV)")

canv.Print('.pdf')
canv.Print('.png')