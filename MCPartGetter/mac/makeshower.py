#
# Example PyROOT script to run analysis module, ana_base.
# The usage is same for inherited analysis class instance.
#

# Load libraries
import os, ROOT, sys
from ROOT import *
from ROOT import gSystem
import time

gSystem.Load("libLArLite_Base")
gSystem.Load("libLArLite_Analysis")
gSystem.Load("libLArLite_LArUtil")
gSystem.Load("libMCPartInfo_MCPartGetter")
gSystem.Load("libBasicTool_GeoAlgo")

# Now import ana_processor & your class. For this example, ana_base.
from ROOT import larlite as fmwk

# Create ana_processor instance
my_proc=fmwk.ana_processor()

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)
#my_proc.set_io_mode(storage_manager.WRITE)
#my_proc.set_io_mode(storage_manager.BOTH)

# Specify what data to read
#my_proc.set_data_to_read(fmwk.DATA.kMCTruth)
#my_proc.set_data_to_read(fmwk.DATA.kMCParticle)

# Set input root file: this is decoder output root file.
# This time, we use a sample file prepared.
#my_proc.add_input_file("./../../../NevisDecoder/Decoder/mac/xmit_subrun_2014_01_13_1_dma_no_1.root")
#my_proc.add_input_file("./../../../NevisDecoder/Decoder/mac/xmit_subrun_2014_01_13_1_trigger.root")
my_proc.add_input_file(sys.argv[1])

# Specify ROOT TDirectory in the file if such structure is present (which is the case for DataScanner output)
#my_proc.set_input_rootdir("scanner")

# Set output root file: this is a separate root file in which your
# analysis module can store anything such as histograms, your own TTree, etc.
my_proc.set_output_file("makeshowers.root")

# Create analysis class instance. For this example, ana_base.
# To show how one can run multiple analysis modules at once,
# we make multiple ana_base instance.

makeshowers = fmwk.MakeShowers()
makeshowers.SetVerbose(False)
makeshowers.SetEcut(0.01)

my_proc.add_process(makeshowers)

# Let's run it.
t0 = int(round(time.time()*1000))
numEvts = 15
my_proc.run(0,numEvts)
t1 = int(round(time.time()*1000))
dt = (t1-t0)/1000. #seconds
print "time diff is {0} sec.".format(dt)
print "time per event is: {0} seconds".format(dt/numEvts)
'''
while my_proc.process_event():
    usrinput = raw_input("Hit Enter: next evt  ||  q: exit viewer\n")
    if ( usrinput == "q" ):
        break
'''
# done!
