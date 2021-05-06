"""
Functions for filtering out frames.

These are the base "cuts": they are always the same. The variable cuts
(the ones you want to fiddle with), are done later in the cuts.py script
(see the 'event_cuts' and 'dom_cuts' dictionaries in cut_options_example.py).
"""

from __future__ import print_function, division

import math

from icecube import dataclasses, phys_services


def in_ice(frame):
    """
    Check that sub_event_stream == 'InIceSplit'.
    """
    event_header = frame['I3EventHeader']
    return event_header.sub_event_stream == 'InIceSplit'


#def min_bias(frame, reco_fit, options):
def min_bias(frame) :
    """
    Check that condition_passed and prescale_passed for FilterMinBias_11 are
    both True. .. moved it to 'filterMinBias_13
    """
    #if frame.Has('QFilterMask'):
    if frame.Has('FilterMask'):
        #filter_mask = frame['QFilterMask_NullSplit0']
        filter_mask = frame['FilterMask']
        filter_min_bias = filter_mask['FilterMinBias_13']
        return filter_min_bias.condition_passed #and filter_min_bias.prescale_passed
#    if frame.Has('PhysMinBiasTriggered') :
#	return frame['PhysMinBiasTriggered'].value
    else:
        return False

def SMT8(frame):
    """
    Check that the length of SplitInIcePulses >= 8.
    """
    pulse_series = frame['SplitInIcePulses'].apply(frame)
  #  frame['SplitInIcePulseSeries'] = pulse_series
    return len(pulse_series) >= 8

# To avoid those frames where MPEFit had already been run, we changed the name of the reconstruction frame objects
# to MPEFitDOMeff

def MPEFit(frame):
    """
    Check that the fit_status of MPEFit is OK and that the zenith is
    between 40 and 70 degrees.
    """
    #  hasattr(go look at this object, see if this atrribute exists)
#    if hasattr(frame['MPEFit'],azimuth) is True:
    if frame.Has('MPEFitDOMeff'):
        mpe = frame['MPEFitDOMeff']
#        angle = math.degrees(mpe.dir.zenith) 
        return (mpe.fit_status == dataclasses.I3Particle.OK) #and 40 < angle < 70)
    else: 
        return False

   #if frame['MPEFit'] is False:
  # 	print "NoMPE in frame"

def muon_zenith(frame,reco_fit):
    """
    Check that the fit_status of the analysis reconstruction is OK and that the zenith is
    between 40 and 70 degrees.
    """
    muon = frame[reco_fit]
    angle = math.degrees(muon.dir.zenith)
    #return (muon.fit_status == dataclasses.I3Particle.OK and 40 < angle < 70)
    return (muon.fit_status == dataclasses.I3Particle.OK)

def FiniteRecoFilter(frame):
    """
    Check that the fit_status of FiniteReco is OK.
    """
    return frame['FiniteRecoFitDOMeff'].fit_status == dataclasses.I3Particle.OK


def InIceSMTTriggered(frame):
    """
    Check that InIceSMTTriggered is True.
    """
    return frame['InIceSMTTriggered'].value
