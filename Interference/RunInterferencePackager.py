# Script for submitting the packager for interference models -- mirrors
# Signal/RunPackager.py, which does the same job for signal workspaces.
import os, sys
import glob
from optparse import OptionParser
from collections import OrderedDict as od

# Import tools from ./tools
from commonTools import *
from commonObjects import *
from tools.submissionTools import *

def get_options():
  parser = OptionParser()
  parser.add_option('--cats', dest='cats', default='', help="Comma separated list of categories")
  parser.add_option('--proc', dest='proc', default='', help="Signal process")
  parser.add_option('--exts', dest='exts', default='', help="Comma separated list of per-era extensions to merge")
  parser.add_option('--outputExt', dest='outputExt', default='packaged', help="Output extension")
  parser.add_option('--mergeYears', dest='mergeYears', default=False, action="store_true", help="Use if merging categories across years")
  parser.add_option('--year', dest='year', default='2022preEE', help="If not merging then add year tag to file name")
  parser.add_option('--batch', dest='batch', default='local', help='Batch')
  parser.add_option('--queue', dest='queue', default='espresso', help='Queue: should not take long (microcentury will do)')
  parser.add_option('--jobOpts', dest='jobOpts', default='', help="Additional options to add to job submission. For Condor separate individual options with a colon (specify all within quotes e.g. \"option_xyz = abc+option_123 = 456\")")
  parser.add_option('--printOnly', dest='printOnly', default=False, action="store_true", help="Dry run: print submission files only")
  return parser.parse_args()
(opt,args) = get_options()

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING INTERFERENCE PACKAGER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
def leave():
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING INTERFERENCE PACKAGER (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  exit(0)

# Store all opts in orderedDict for submissionTools
options = od()
options['cats'] = opt.cats
options['proc'] = opt.proc
options['ext'] = opt.outputExt
options['mode'] = 'packageInterference'
options['modeOpts'] = '--exts %s'%opt.exts
if opt.mergeYears: options['modeOpts'] += ' --mergeYears'
else: options['modeOpts'] += ' --year %s'%opt.year
options['batch'] = opt.batch
options['queue'] = opt.queue
options['jobOpts'] = opt.jobOpts
options['printOnly'] = opt.printOnly
options['nCats'] = len(options['cats'].split(","))

print(" --> Packaging interference workspaces from: %s"%opt.exts)
print(" --> For analysis categories: %s"%options['cats'])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make directory to store job scripts and output
if not os.path.isdir("%s/results/outdir_%s"%(iwd__,options['ext'])): os.system("mkdir -p %s/results/outdir_%s"%(iwd__,options['ext']))

# Write submission files: style depends on batch system
writeSubFiles(options)
print("  --> Finished writing submission scripts")

# Submit scripts to batch system
if not options['printOnly']:
  submitFiles(options)
else:
  print("  --> Running with printOnly option. Will not submit scripts")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
leave()
