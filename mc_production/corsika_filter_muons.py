 #!/usr/bin/env python

# import required icecube-related stuff
from icecube import icetray, dataclasses, simclasses, dataio
from icecube.icetray import I3Units
from I3Tray import I3Tray

# command line options required to configure the simulation
from argparse import ArgumentParser
from os.path import expandvars

import os, sys

from icecube import phys_services

from icecube import corsika

parser = ArgumentParser()
parser.add_argument("--eprimarymin", type=int, dest="eMIN", help="Minimum energy per particle or per nucleon for the cosmic-rays")
parser.add_argument("--corsikaVersion", type=str, dest="VERSION", help="version of corsika to run")
parser.add_argument("--CORSIKAseed", type=int, dest="CORSIKAseed", help="CORSIKA seed")
parser.add_argument("--RunCorsika", type=bool, dest="RunCorsika", help="Run CORSIKA or only generate INPUTS file")
parser.add_argument("--seed", type=int, dest="SEED", help="Initial seed for the random number generator")
parser.add_argument("--HistogramFilename", type=str, dest="HistoF", help="Histogram filename")
parser.add_argument("--CutoffType", type=str, dest="CUTOFF", help="Sets SPRIC=T (EnergyPerNucleon) or F (EnergyPerParticle)")
parser.add_argument("--UsePipe", type=bool, dest="UsePipe", help="Use pipe for corsika output")
parser.add_argument("--outputfile", type=str, dest="OUTPUT", help="Output filename")
parser.add_argument("--eprimarymax", type=int, dest="eMAX", help="Minimum energy per particle or per nucleon for the cosmic-rays")
parser.add_argument("--oversampling", type=int, dest="OVERSAMPLE", help="How many times the variation of impact parameters of rays are reproduced")
parser.add_argument("--procnum", type=int, dest="JOB", help="process number")
parser.add_argument("--compress", type=bool, dest="COMP", help="compress corsika output")
parser.add_argument("--pnorm", type=list, dest="pNORM", help="5-component relative contribution H, He, N, Al, Fe")
parser.add_argument("--RepoURL", type=str, dest="URL", help="URL of repository containing corsika tarballs")
parser.add_argument("--nproc", type=int, dest="nPROC", help="Number of processes for (RNG)")
parser.add_argument("--gcdfile", type=str, dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_argument("--nshowers", type=int, dest="SHOWER", help="Number of generated CR showers")
parser.add_argument("--pgam", type=list, dest="pGAM", help="5-component spectral indices H, He, N, Al, Fe")
parser.add_argument("--summaryfile", type=str, dest="SUMMARY", help="Filename of summary file")
parser.add_argument("--EnableHistogram", type=bool, dest="HISTO", help="Write a SanityChecker histogram file")

# parse cmd line args, bail out if anything is not understood
args = parser.parse_args()

tray = I3Tray()

# randomService = phys_services.I3SPRNGRandomService(
#     seed = args.SEED,
#     nstreams = 100000000,
#     streamnum = args.RUNNUMBER)

tray.context['I3RandomService'] = randomService

tray.AddModule('I3Reader', 'reader',  filenamelist = [args.GCDFILE])

tray.AddSegment(corsika.Corsika5ComponentGenerator,"makeCRs",
                nshowers=args.SHOWER,
                procnum=args.JOB,
                nproc=args.nPROC,
                seed=args.SEED,
                gcdfile=args.GCDFILE,
                outputfile=args.OUTPUT,
                RunCorsika=args.RunCorsika,
                sumaryfile=args.SUMMARY,
                pnorm=args.pNORM,
                pgam=args.pGAM,
                CORSIKAseed=args.CORSIKAseed,
                eprimarymax=args.eMAX,
                eprimarymin=args.eMIN,
                OverSampling=args.OVERSAMPLE,
                corsikaVersion=args.VERSION,
                CutoffType=args.CUTOFF,
                RepoURL=args.URL,
                UsePipe=args.UsePipe,
                compress=args.COMP,
                HistogramFilename=args.HistoF,
                EnableHistogram=args.HISTO
               )

tray.AddModule('I3Writer', 'writer',
    Streams=[icetray.I3Frame.Stream('S'), icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ],
    filename=args.OUTPUT)

tray.AddModule('TrashCan', 'YesWeCan')
print("Executing...")
tray.Execute()
print("Finish!")
tray.Finish()