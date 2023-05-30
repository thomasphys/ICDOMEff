#! /usr/bin/env python

"""
Process the I3 files.
"""

from __future__ import print_function, division  # 2to3

# Kludge to allow importing from parent directory for shared utility modules
import os
import sys
sys.path.append('/home/tmcelroy/icecube/domeff')
import inspect
import numpy as np
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from icecube.filterscripts.offlineL2 import Globals
#from icecube.filterscripts.offlineL2.Rehydration import Rehydration, Dehydration
from icecube.filterscripts.offlineL2.PhotonTables import InstallTables
from icecube.filterscripts.offlineL2.ClipStartStop import ClipStartStop
#from icecube.phys_services.which_split import which_split

# Setup logging
from icecube.icetray.i3logging import *
from icecube import dataio, icetray, finiteReco, gulliver, simclasses, dataclasses, photonics_service, phys_services, spline_reco #, MuonGun
from icecube.common_variables import direct_hits, hit_multiplicity, hit_statistics
from icecube.filterscripts.offlineL2.level2_HitCleaning_WIMP import WimpHitCleaning
from I3Tray import I3Tray, I3Units, load
#from filters_InIceSplit import in_ice, min_bias, SMT8, MPEFit, InIceSMTTriggered
from filters_InIceSplit_2015 import in_ice, min_bias, SMT8, MPEFit, InIceSMTTriggered, FiniteRecoFilter, muon_zenith
from general import get_truth_muon, get_truth_endpoint, count_hits, reco_endpoint, move_cut_variables, totaltimefilter,timestartfilter, tot_charge, movellhparams
from geoanalysis import calc_dist_to_border, calc_dist_to_border_mctruth
from domanalysis import dom_data
from writeEvent import EventWriter
from writeEventRecoPulses import EventWriterRecoPulses 
import argparse
# Reconstructions
from icecube.filterscripts.offlineL2 import Globals
from icecube.filterscripts.offlineL2.level2_Reconstruction_WIMP import FiniteReco
from icecube.filterscripts.offlineL2.level2_Reconstruction_Muon import SPE, MPE
from icecube.filterscripts.offlineL2.PhotonTables import InstallTables
from icecube import cramer_rao
from icecube.filterscripts.shadowfilter import ShadowFilter
import ROOT

load('libipdf')
load('libgulliver')
load('libgulliver-modules')
load('liblilliput')
load('libstatic-twc')
#load('libjeb-filter-2012')
load('libfilterscripts')

def printtag(frame,message) :
	print(message)
	return True

eventcount1 = 0.0
def countevents1(frame) :
	global eventcount1
	eventcount1 += 1.0

eventcount2 = 0.0
def countevents2(frame) :
	global eventcount2
        eventcount2 += 1.0

rand = ROOT.TRandom3(0)
def prescalemodule(frame,prescalefrac) :
  global rand
  if rand.Uniform() < prescalefrac :
    return True
  return False

def GetNewTrack(reco,truth,scale):

	#Get closest approach point for Reco
        l_reco = -(reco.pos.x*reco.dir.x+reco.pos.y*reco.dir.y+reco.pos.z*reco.dir.z)
        point1_reco = [reco.pos.x+(l_reco-100.)*reco.dir.x,reco.pos.y+(l_reco-100.)*reco.dir.y,reco.pos.z+(l_reco-100.)*reco.dir.z]
        point2_reco = [reco.pos.x+(l_reco+100.)*reco.dir.x,reco.pos.y+(l_reco+100.)*reco.dir.y,reco.pos.z+(l_reco+100.)*reco.dir.z]

        #Get closest approach point for truth
        l_truth = -(truth.pos.x*truth.dir.x+truth.pos.y*truth.dir.y+truth.pos.z*truth.dir.z)
        point1_truth = [truth.pos.x+(l_truth-100.)*truth.dir.x,truth.pos.y+(l_truth-100.)*truth.dir.y,truth.pos.z+(l_truth-100.)*truth.dir.z]
        point2_truth = [truth.pos.x+(l_truth+100.)*truth.dir.x,truth.pos.y+(l_truth+100.)*truth.dir.y,truth.pos.z+(l_truth+100.)*truth.dir.z]

        #compute vectors between points
        vector1 = [point1_truth[0]-point1_reco[0],point1_truth[1]-point1_reco[1],point1_truth[2]-point1_reco[2]]
        vector2 = [point2_truth[0]-point2_reco[0],point2_truth[1]-point2_reco[1],point2_truth[2]-point2_reco[2]]

        newpoint1 = [point1_reco[0]+scale*vector1[0],point1_reco[1]+scale*vector1[1],point1_reco[2]+scale*vector1[2]]
        newpoint2 = [point2_reco[0]+scale*vector2[0],point2_reco[1]+scale*vector2[1],point2_reco[2]+scale*vector2[2]]

        newdirection = [newpoint2[0]-newpoint1[0],newpoint2[1]-newpoint1[1],newpoint2[2]-newpoint1[2]]
        magnitude = np.sqrt(newdirection[0]**2+newdirection[1]**2+newdirection[2]**2)
        newdirection = [newdirection[0]/magnitude,newdirection[1]/magnitude,newdirection[2]/magnitude]

        centerpoint = [(newpoint1[0]+newpoint2[0])*0.5,(newpoint1[1]+newpoint2[1])*0.5,(newpoint1[2]+newpoint2[2])*0.5]
        newvertex = [centerpoint[0]-l_reco*newdirection[0],centerpoint[1]-l_reco*newdirection[1],centerpoint[2]-l_reco*newdirection[2]]

        return dataclasses.I3Direction(newdirection[0],newdirection[1],newdirection[2]), dataclasses.I3Position(newvertex[0],newvertex[1],newvertex[2])

def GetMuonTrack(MMCTrackList) :

        mindist = 100000.
        mindistI3Particle = None

        for track in MMCTrackList :
                xc = track.xc
                yc = track.yc
                zc = track.zc

                if np.sqrt(xc**2+yc**2+zc**2) < mindist :
                        mindistI3Particle = track.GetI3Particle()
                        mindist = np.sqrt(xc**2+yc**2+zc**2)

        return mindistI3Particle

def shifttrack(frame,trackname,pulseseries) :

	DC_Strings_all = [79,80,81,82,83,84,85,86]

	biasfactor = 0.0
	biasstep = 1.0/19.
	boost = 0.0
	booststep = 1.0/19.
	truthtrack = None
	if frame.Has("MMCTrackList"):
		truthtrack = GetMuonTrack(frame["MMCTrackList"])
	else :
		frame[trackname+'_shift'] = frame[trackname]

        newtrack = dataclasses.I3Particle()
	maxChargeDOM = None
	maxcharge = 0.0
	totalcharge_IC = 0.0
	totalcharge_DC = 0.0
	
	pulsemap = frame[pulseseries].apply(frame)
	for dom in pulsemap.keys() :
		charge = 0.0
		for pulse in pulsemap[dom] :
			charge += pulse.charge
		if charge > maxcharge :
			maxcharge = charge
			maxChargeDOM = dom
		if dom.string in DC_Strings_all:
			totalcharge_DC += charge
		else :
			totalcharge_IC += charge
		
	totalcharge = totalcharge_IC + totalcharge_DC
	shift = (booststep*9.0)*(totalcharge_DC/totalcharge) + (biasstep*3.0)*(totalcharge_IC/totalcharge)

	newtrack = frame[trackname]
	
        newtrack.dir,newtrack.pos = GetNewTrack(newtrack,truthtrack,shift)
        #newtrack.shape = dataclasses.I3Particle.InfiniteTrack
        frame[trackname+'_shift'] = newtrack
	return True

def standardcuts(frame) :
    
  if frame['SplineMPE'].pos.z < -450.0 : return False

  if frame['DistToBorder'].value < 50. : return False
  
  if frame['FiniteRecoLLHRatio'].value < 5.0 : return False

  if frame['rlogl'].value < 5 : return False

  if frame['ICNHits'].value  > 30 : return False

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--datadir', help='Directory of data files.',type=str,
        default = '/data/user/sanchezh/IC86_2015/Final_Level2_IC86_MPEFit_')
parser.add_argument('-g', '--gcd', help='Geometry file.', type = str,
        default = "${I3_DATA}/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE.i3.gz")
parser.add_argument('-o', '--output', help='Name of output file.', type=str,
        default = "out.root")
parser.add_argument('-z', '--zenithrange', help='Range of muon Zeniths', type = float,
        nargs = '+',  default = [-180.0,180.0])
parser.add_argument('-q', '--energyrange', help='Range of muon Energies', type = float,
        nargs = "+", default = [0.0, 9999999.00])
parser.add_argument('-s', '--sim', help='Is file simulation', type = bool, default = False)
parser.add_argument('-n', '--nevents', help='Number of events to process.', type = int, default = -1)
parser.add_argument('-t', '--datafiletype', help='Suffix of Datafiles', type = str, default = 'i3.bz2')
parser.add_argument('-r', '--runnum', help='number to identify target file', type = int, default = 0)
parser.add_argument('-p', '--pulsename', help='Name of new pulse list', type = str, default = 'SRTInIcePulsesDOMEff')
parser.add_argument('-m', '--maxdist', help='maximum distance to DOM to consider', type = float, default = 200.0)
parser.add_argument('-x', '--dstfile', help='dts file, should be .root',type=str,default = "")
parser.add_argument('-y','--prescale',help='prescale fraction',type=float,default = 1.0)

args = parser.parse_args()

dom_data_options = {}
#    options['pulses_name'] = 'SplitInIcePulses' 
dom_data_options['pulses_name'] = args.pulsename
dom_data_options['max_dist'] = args.maxdist

tray = I3Tray()

datafilename = "{0:0{1}d}".format(args.runnum,5)

datafilename2 = "{0:0{1}d}_hitdist".format(args.runnum,5)


# Read the files.
tray.AddModule('I3Reader', 'I3Reader',
               Filenamelist=[args.gcd, args.datadir+datafilename+args.datafiletype])

if args.prescale < 1.0 :
  tray.AddModule(prescalemodule,'prescale',prescalefrac=args.prescale)

tray.AddModule(countevents1,"count1")

#if "PFDSTnoSPE" in args.datadir :
#  tray.AddModule(printtag, 'printtag_pff',message = "pff data")

#  tray.AddModule(ClipStartStop, 'clipstartstop')

#  tray.AddModule(printtag, 'printtag_clip',message = "pass clipstartstop")

#  tray.AddSegment(Rehydration, 'rehydrator',
#                #dstfile=args.dstfile,
#                mc=False,
#                doNotQify=False,
#                pass2=True,
#                )
#  tray.AddModule(printtag, 'printtag_rehydrate',message = "pass rehydrate")

#tray.AddModule(printtag, 'printtag_newevent',message = "new event")
# Filter the ones with sub_event_stream == InIceSplit
tray.AddModule(in_ice, 'in_ice')
#tray.AddModule(printtag, 'printtag_in_ice',message = "passed in_ice")
# Make sure that the length of SplitInIcePulses is >= 8

if args.sim :

	tray.AddSegment(ShadowFilter, "MoonAndSun",
			mcseed=args.runnum)

tray.AddModule('TriggerCheck_13', 'TriggerCheck_13',
               I3TriggerHierarchy='I3TriggerHierarchy',
               InIceSMTFlag='InIceSMTTriggered',
               IceTopSMTFlag='IceTopSMTTriggered',
               InIceStringFlag='InIceStringTriggered',
               PhysMinBiasFlag='PhysMinBiasTriggered',
               PhysMinBiasConfigID=106,
               DeepCoreSMTFlag='DeepCoreSMTTriggered',
               DeepCoreSMTConfigID=1010)

	# Check that InIceSMTTriggered is true.
#tray.AddModule(printtag, 'printtag_trigcheck',message = "passed trigcheck")
tray.AddModule(InIceSMTTriggered, 'InIceSMTTriggered')
#tray.AddModule(printtag, 'printtag_smttrig',message = "passed smttrig")
tray.AddModule(SMT8, 'SMT8')
#tray.AddModule(printtag, 'printtag_smt8',message = "passed smt8")
#tray.AddModule(min_bias, 'min_bias')
#tray.AddModule(printtag, 'printtag_minbias',message = "passed minbias")

# Generate RTTWOfflinePulses_FR_WIMP, used to generate the finite reco reconstruction in data

#---- Generate filtered pulse series with same configuration as used by Nick and Sebastian -----
# Generate the SRTInIcePulses, which are used for running basic reconstruction algorithms on data

from icecube.STTools.seededRT.configuration_services import I3DOMLinkSeededRTConfigurationService
seededRTConfig = I3DOMLinkSeededRTConfigurationService(
    ic_ic_RTRadius=150.0 * I3Units.m,
    ic_ic_RTTime=1000.0 * I3Units.ns,
    treat_string_36_as_deepcore=False,
    useDustlayerCorrection=False,
    allowSelfCoincidence=True
)

# Notice that we name the pulse series SRTInIcePulsesDOMeff to avoid
# repeating frame objects in the frames that originally had reconstuctions
tray.AddModule('I3SeededRTCleaning_RecoPulseMask_Module', 'North_seededrt',
               InputHitSeriesMapName='SplitInIcePulses',
               OutputHitSeriesMapName=args.pulsename,
               STConfigService=seededRTConfig,
               SeedProcedure='HLCCoreHits',
               NHitsThreshold=2,
               MaxNIterations=3,
               Streams=[icetray.I3Frame.Physics],
               If=lambda f: True
               )

#tray.AddModule(EventWriterRecoPulses,"EventWriterRecoPulses",
#               FileName=args.output+datafilename2+'.h5',
#               pulseseries = args.pulsename)


# ---- Linefit and SPEfit ---------------------------------------------------
tray.AddSegment(SPE,'SPE',
                If = lambda f: True,
                Pulses=args.pulsename,
                suffix='DOMeff',
                LineFit= 'LineFit',
                SPEFitSingle = 'SPEFitSingle',
                SPEFit = 'SPEFit2',
                SPEFitCramerRao = 'SPEFit2CramerRao',
                N_iter = 2
                )

#tray.AddModule(printtag, 'printtag_SPE',message = "passed SPE")
# ---- MPEFit reconstruction ------------------------------------------------
tray.AddSegment(MPE, 'MPE',
                Pulses = args.pulsename,
                Seed = 'SPEFit2',
                #If = which_split(split_name='InIceSplit') & (lambda f:muon_wg(f)),
                If = lambda f: True,
                suffix='DOMeff',
                MPEFit = 'MPEFit',
                MPEFitCramerRao = 'MPEFitCramerRao'
                )

#tray.AddModule(printtag, 'printtag_MPE',message = "passed MPE")

tray.AddModule(MPEFit, 'MPEFit')

#tray.AddModule(printtag, 'printtag_MPEFit',message = "passed MPEfit")

# -----Spline Reco -------------------------------------------------------
#spline paths Madison
timingSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_prob_z20a10_V2.fits'
amplitudeSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_abs_z20a10_V2.fits'
stochTimingSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfHighEStoch_mie_prob_z20a10.fits'
stochAmplitudeSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfHighEStoch_mie_abs_z20a10.fits'
pulses = "SRTOfflinePulses"
EnEstis = ["SplineMPETruncatedEnergy_SPICEMie_AllDOMS_Muon",
           "SplineMPETruncatedEnergy_SPICEMie_DOMS_Muon",
           "SplineMPETruncatedEnergy_SPICEMie_AllBINS_Muon",
           "SplineMPETruncatedEnergy_SPICEMie_BINS_Muon",
            "SplineMPETruncatedEnergy_SPICEMie_ORIG_Muon"
           ]

# splineMPE with default configuration!
tray.AddSegment(spline_reco.SplineMPE, "SplineMPE",
                configuration="default",
                PulsesName= args.pulsename,
                TrackSeedList=["MPEFitDOMeff"],
                BareMuTimingSpline=timingSplinePath,
                BareMuAmplitudeSpline=amplitudeSplinePath,
                fitname="SplineMPE",
                )

tray.AddModule(shifttrack,"shifttrack",
		trackname = "SplineMPE",
		pulseseries = args.pulsename)

#tray.AddModule(printtag, 'printtag_splineMPE',message = "passed SplineMPE")
#Thomas, will apply this later, want to study impact.
tray.AddModule(muon_zenith, 'MuonZenithFilter',
               reco_fit='SplineMPE_shift')
    

# -----Finite Reco------------------------------------------------------------
tray.AddSegment(InstallTables, 'InstallPhotonTables')
#tray.AddSegment(FiniteReco,'FiniteReco',
#                suffix = 'DOMeff',
#                InputTrackName = 'SplineMPE',
#		Pulses = args.pulsename,
                #Pulses = 'RTTWOfflinePulses_FR_WIMP_DOMeff',
#                Pulses = 'SRTInIcePulses'
#		)


ic86 = [ 21, 29, 39, 38, 30, 40, 50, 59, 49, 58, 67, 66, 74, 73, 65, 72, 78, 48, 57, 47,
	46, 56, 63, 64, 55, 71, 70, 76, 77, 75, 69, 60, 68, 61, 62, 52, 44, 53, 54, 45,
	18, 27, 36, 28, 19, 20, 13, 12, 6, 5, 11, 4, 10, 3, 2, 83, 37, 26, 17, 8, 9, 16,
	25, 85, 84, 82, 81, 86, 35, 34, 24, 15, 23, 33, 43, 32, 42, 41, 51,
	31, 22, 14, 7, 1, 79, 80] # Taken from http://wiki.icecube.wisc.edu/index.php/Deployment_order

icetray.load('finiteReco', False)
tray.AddService('I3GulliverFinitePhPnhFactory', 'GulliverPhPnh',
		InputReadout = args.pulsename,
		PhotorecName = Globals.PhotonicsServiceFiniteReco,
		ProbName = 'PhPnhPhotorec', # Use photorec tables to calculate probabilities
		RCylinder = 200, # Radius around the track in which probabilities are considered
		SelectStrings = ic86,)

#first guess to start-stop points
tray.AddModule('I3StartStopPoint', 'GulliverVertexReco',
		Name = 'SplineMPE_shift', 
		InputRecoPulses = args.pulsename,
		ExpectedShape = 70, # Contained track, this way the start AND stop point are reconstructed
		CylinderRadius = 200, # Cylinder radius for the cut calculation,  take care to use Cylinder Radius ==200
		If = lambda f: True,)

#tray.AddModule(printtag, 'printtag_startstoppoint',message = "passed I3StartStopPoint")
#starting/stopping probability
tray.AddModule('I3StartStopLProb', 'GulliverFiniteRecoLlh',
	Name = 'SplineMPE_shift_Finite', # Name of the input track with _Finite added from I3StartStopPoint
	ServiceName = 'GulliverPhPnh', # Name of the service is the instance name from I3GulliverFinitePhPnhFactory
	If = lambda f: True,)

#tray.AddModule(printtag, 'printtag_LProb',message = "passed LProb")
	
# rename the finiteReco track
tray.AddModule('Rename', 'Gulliver_FiniteRecoRename',
		Keys = ['SplineMPE_shift_Finite', 'FiniteRecoFitDOMeff'],)
	
# rename the finitRecoLlh
tray.AddModule('Rename', 'Gulliver_FiniteRecoLlhRename',
		Keys = ['GulliverFiniteRecoLlh_StartStopParams', 'FiniteRecoLlhDOMeff'],)  # Name is instance name of I3StartStopLProb with _StartStopParams added
	                 		
# rename the finiteReco cuts
tray.AddModule('Rename', 'Gulliver' + '_FiniteRecoCutsRename' + 'DOMeff',
		Keys = ['SplineMPE_shift_FiniteCuts', 'FiniteRecoCutsDOMeff'],) # Name of the input track with _FiniteCuts added from I3StartStopPoint

tray.AddModule(FiniteRecoFilter, 'FiniteRecoFilter')
#tray.AddModule(printtag, 'printtag_FiniteRecoFilter',message = "passed FiniteRecoFilter")

tray.AddModule(movellhparams, "MoveLLHParams",
		llhparams = 'FiniteRecoLlhDOMeff',	      
)

#tray.AddModule(printtag, 'printtag_move',message = "passed MOveLLHParams")


# -----Endpoint---------------------------------------------------------------
# Add the reconstructed event endpoint to the frame.
tray.AddModule(reco_endpoint, 'reco_endpoint',
               endpoint_fit='FiniteRecoFitDOMeff'
               )

#tray.AddModule(printtag, 'printtag_reco_endpoint',message = "passed reco_endpoint")

tray.AddModule(tot_charge,'tot_charge',
                reco_fit='SplineMPE_shift',
                pulses = args.pulsename,
              )

#tray.AddModule(printtag, 'printtag_tot_charge',message = "passed tot_charge")

# DOManalysis
# This uses the MPEFit's to calculate TotalCharge, RecoDistance, etc.
tray.AddModule(dom_data, 'dom_data',
               #reco_fit='MPEFitDOMeff',
               reco_fit='SplineMPE_shift',
               options=dom_data_options
               )
#tray.AddModule(printtag, 'printtag_dom_data',message = "passed dom_data")
# General

# Calculate cut variables
tray.AddSegment(direct_hits.I3DirectHitsCalculatorSegment, 'I3DirectHits',
                PulseSeriesMapName=args.pulsename,
#                ParticleName='MPEFitDOMeff',
                ParticleName='SplineMPE_shift',
#                OutputI3DirectHitsValuesBaseName='MPEFitDOMeffDirectHits')
                OutputI3DirectHitsValuesBaseName='SplineMPEDirectHits'
                )

#tray.AddModule(printtag, 'printtag_directhits',message = "passed directhits")

tray.AddSegment(hit_multiplicity.I3HitMultiplicityCalculatorSegment, 'I3HitMultiplicity',
                PulseSeriesMapName=args.pulsename,
                OutputI3HitMultiplicityValuesName='HitMultiplicityValues'
                )

#tray.AddModule(printtag, 'printtag_hitmultiplicity',message = "pssed hitmultiplicity")
 
tray.AddSegment(hit_statistics.I3HitStatisticsCalculatorSegment, 'I3HitStatistics',
                PulseSeriesMapName=args.pulsename,
                OutputI3HitStatisticsValuesName='HitStatisticsValues'
                )
#tray.AddModule(printtag, 'printtag_hitstats',message = "passed hit_stats")
# Move the cut variables into the top level of the frame.
tray.AddModule(move_cut_variables, 'move_cut_variables',
              # direct_hits_name='MPEFitDOMeffDirectHits',
              # fit_params_name='MPEFitDOMeffFitParams')
               direct_hits_name='SplineMPEDirectHits',
               fit_params_name='SplineMPEFitParams'
               )
#tray.AddModule(printtag, 'printtag_movecuts',message = "passed move_cuts")
# Calculate ICAnalysisHits, DCAnalysisHits, ICNHits, and DCNHits
tray.AddModule(count_hits, 'count_hits',
               pulses_name=args.pulsename)
 
#tray.AddModule(printtag, 'printtag_hitcount',message = "passed hitcount")

# Geoanalysis
# Calculate the distance of each event to the detector border.
tray.AddModule(calc_dist_to_border, 'calc_dist_to_border')

tray.AddModule(standardcuts,'standardcutsmodule')

#tray.AddModule(printtag, 'printtag_dist to border',message = "pssed dist_to_border")

#if args.sim :
        # Count the number of in ice muons and get the truth muon
        #tray.AddModule(get_truth_muon, 'get_truth_muon')
        #tray.AddModule(get_truth_endpoint, 'get_truth_endpoint')
        #tray.AddModule(calc_dist_to_border_mctruth,'calc_dist_to_border_mctruth')
	#tray.AddModule(printtag, 'printtag_simstuff',message = "passed sim stuff")

# Write the data out to an HDF5 analysis file

tray.AddModule(EventWriter, 'EventWriter',
               FileName=args.output+datafilename+'.h5',
		splinempename='SplineMPE_shift')

#tray.AddModule(printtag, 'printtag_writer',message = "passed writer")
# Write out the data to an I3 file
#tray.AddModule('I3Writer', 'I3Writer',
#               FileName=args.output+datafilename+'.i3.gz',
               #SkipKeys=['InIceRecoPulseSeriesPattern.*'],
#               DropOrphanStreams=[icetray.I3Frame.DAQ],
#               Streams=[icetray.I3Frame.TrayInfo,icetray.I3Frame.DAQ,icetray.I3Frame.Physics,icetray.I3Frame.Simulation]
#               )
    
tray.AddModule('TrashCan', 'yeswecan')
if args.nevents > 0 :
  tray.Execute(opts['nevents'])
else :
  tray.Execute()
tray.Finish()
print("events %f / %f passed" % (eventcount1,eventcount2)) 
