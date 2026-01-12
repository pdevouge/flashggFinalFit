#!/usr/bin/env python3
import ROOT
import argparse

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# ----------------------------
# Parse options
# ----------------------------
parser = argparse.ArgumentParser(description="Plot signal PDF with multiple parameter sets")
parser.add_argument("--inputWSFile", required=True, help="Input RooWorkspace file")
parser.add_argument("--workspaceName", default="wsig_13TeV", help="Name of the workspace")
parser.add_argument("--pdfName", default="model_s", help="Name of the signal PDF")
parser.add_argument("--xvar", default="CMS_hgg_mass", help="Observable name in workspace")
parser.add_argument(
    "--parameterSets",
    default="",
    help='Semicolon-separated parameter sets: "p1=v1,p2=v2; p1=v1b,p2=v2b"'
)
parser.add_argument("--snapshot", default=None, help="Snapshot name to load")
parser.add_argument("--mass", type=float, default=None, help="Set MH value if needed")
parser.add_argument("--nBins", type=int, default=1600, help="Number of histogram bins")
parser.add_argument("--out", default="signal.pdf", help="Output plot file")
args = parser.parse_args()

# ----------------------------
# Load workspace
# ----------------------------
f = ROOT.TFile.Open(args.inputWSFile)
w = f.Get(args.workspaceName)
if not w:
    raise RuntimeError("Workspace not found")

if args.snapshot:
    print(f"Loading snapshot: {args.snapshot}")
    w.loadSnapshot(args.snapshot)

# ----------------------------
# Get observable + PDF
# ----------------------------
x = w.var(args.xvar)
if not x:
    raise RuntimeError(f"Observable {args.xvar} not found")

pdf = w.pdf(args.pdfName)
if not pdf:
    raise RuntimeError(f"PDF {args.pdfName} not found")

# Set MH (optional, common to all curves)
if args.mass is not None and w.var("MH"):
    w.var("MH").setVal(args.mass)

# ----------------------------
# Prepare frame
# ----------------------------
xmin, xmax = x.getMin(), x.getMax()
print(f"Plotting in range [{xmin}, {xmax}]")

x.setRange(650,750)
frame = x.frame()

# ----------------------------
# Loop over parameter sets
# ----------------------------
param_sets = []
if args.parameterSets:
    param_sets = [s.strip() for s in args.parameterSets.split(";") if s.strip()]
else:
    param_sets = [""]  # default: current workspace values

colors = [
    ROOT.kBlack,
    ROOT.kRed + 1,
    ROOT.kBlue + 1,
    ROOT.kGreen + 2,
    ROOT.kMagenta + 1,
    ROOT.kOrange + 7,
]

legend = ROOT.TLegend(0.1, 0.75, 0.88, 0.88)
legend.SetTextSize(0.02)
legend.SetBorderSize(0)
legend.SetFillStyle(0)

for i, pset in enumerate(param_sets):
    # Restore snapshot before modifying parameters
    if args.snapshot:
        w.loadSnapshot(args.snapshot)

    label = "default"
    if pset:
        label = pset
        print(f"Setting parameters for curve {i}:")
        for assignment in pset.split(","):
            name, value = assignment.split("=")
            var = w.var(name.strip())
            if not var:
                raise RuntimeError(f"Parameter {name} not found in workspace")
            var.setVal(float(value))
            print(f"  {name} = {value}")

    color = colors[i % len(colors)]

    pdf.plotOn(
        frame,
        ROOT.RooFit.Binning(args.nBins),
        ROOT.RooFit.LineColor(color),
        ROOT.RooFit.Name(f"curve_{i}")
    )

    legend.AddEntry(frame.findObject(f"curve_{i}"), label, "l")

# ----------------------------
# Draw
# ----------------------------
c = ROOT.TCanvas("c", "c", 800, 600)
frame.Draw()
legend.Draw()
c.SaveAs(args.out)

print(f"Saved: {args.out}")
