from optparse import OptionParser
import ROOT
from CombineHarvester.CombineTools.plotting import *
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)


def get_options():
  parser = OptionParser()
  # Take inputs from config file
  parser.add_option('--input', dest='inputJson', default='limits_default.json', help="Limits.json")
  parser.add_option('--unblinded', dest='unblinded', action='store_true', help="Unblind limit plot")
  return parser.parse_args()
(opt,args) = get_options()


# Style and pads
ModTDRStyle()
canv = ROOT.TCanvas('limit', 'limit')
pads = OnePad()

# Get limit TGraphs as a dictionary
draw = ['obs', 'exp0', 'exp1', 'exp2'] if opt.unblinded else ['exp0', 'exp1', 'exp2']
graphs = StandardLimitsFromJSONFile(opt.inputJson, draw=draw)

# Create an empty TH1 from the first TGraph to serve as the pad axis and frame
axis = CreateAxisHist(list(graphs.values())[0])
axis.GetXaxis().SetTitle('m_{X} (GeV)')
axis.GetYaxis().SetTitle('95% CL limit on #it{#sigma#times#bf{B}} [pb]')
pads[0].cd()
axis.Draw('axis')

# Create a legend in the top left
legend = PositionedLegend(0.3, 0.2, 3, 0.015)

# Set the standard green and yellow colors and draw
StyleLimitBand(graphs)
DrawLimitBand(pads[0], graphs, legend=legend)
legend.Draw()

# Re-draw the frame and tick marks
pads[0].RedrawAxis()
pads[0].GetFrame().Draw()

# Adjust the y-axis range such that the maximum graph value sits 25% below
# the top of the frame. Fix the minimum to zero.
FixBothRanges(pads[0], 0, 0, GetPadYMax(pads[0]), 0.25)

# Standard CMS logo
DrawCMSLogo(pads[0], 'CMS', 'Internal', 0, 0.13, 0.035, 1.2, '', 0.8)
DrawCMSLogo(pads[0], 'Spin-0', '2022PostEE', 11, 0.2, 0.035, 1.2, '', 0.8)

# Add luminosity text at the top right
lumi = ROOT.TLatex()
lumi.SetNDC()
lumi.SetTextFont(42)
lumi.SetTextSize(0.04)
lumi.SetTextAlign(31)
lumi.DrawLatex(0.95, 0.96, "138 fb^{-1} (13.6 TeV)")

canv.Print('.pdf')
canv.Print('.png')