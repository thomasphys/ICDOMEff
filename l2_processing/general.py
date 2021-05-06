"""
General functions that don't fit into one category.
"""
from __future__ import print_function, division
import time
from icecube import dataclasses, finiteReco
from icecube.dataclasses import I3Constants
from icecube.phys_services import I3Calculator as calc

def timestartMPE(frame):
    t = time.time()
    frame['TimeStartMPE'] = dataclasses.I3Double(t)
def totaltimeMPE(frame):
    t = time.time()
    tstart = frame['TimeStartMPE'].value
    frame['TotalTimeMPE'] = dataclasses.I3Double(t-tstart)

def timestartSpline(frame):
    t = time.time()
    frame['TimeStartSpline'] = dataclasses.I3Double(t)
def totaltimeSpline(frame):

    t = time.time()
    tstart = frame['TimeStartSpline'].value
    frame['TotalTimeSpline'] = dataclasses.I3Double(t-tstart)

def timestartall(frame):
    t = time.time()
    frame['TimeStartall'] = dataclasses.I3Double(t)
def totaltimeall(frame):

    t = time.time()
    tstart = frame['TimeStartall'].value
    frame['TotalTimeall'] = dataclasses.I3Double(t-tstart)

def timestartfilter(frame):
    t = time.time()
    frame['TimeStartFilter'] = dataclasses.I3Double(t)
def totaltimefilter(frame):

    t = time.time()
    tstart = frame['TimeStartFilter'].value
    frame['TotalTimeFilter'] = dataclasses.I3Double(t-tstart)

def movellhparams(frame,llhparams):
# Moving FiniteReco LLH parameters to write them out'''

    frame['FiniteRecoLLHRatio'] = dataclasses.I3Double(float(frame[llhparams].LLHStoppingTrack) - float(frame[llhparams].LLHInfTrack))

def tot_charge(frame,reco_fit,pulses):
    n_ice_group = I3Constants.n_ice_group
    n_ice_phase = I3Constants.n_ice_phase
    pulse_series = frame[pulses].apply(frame)
    dom_geo = frame['I3Geometry'].omgeo.items()
    total_charge = 0
    for dom, geo in dom_geo:
        dom_position = geo.position
        mpe = frame[reco_fit]
        if dom in pulse_series.keys():
            for pulse in pulse_series[dom]:
                time_res = calc.time_residual(mpe, dom_position, pulse.time, n_ice_group, n_ice_phase)
                if time_res < 1000:
                    total_charge += pulse.charge
    frame['EventCharge'] = dataclasses.I3Double(total_charge)

def get_truth_muon(frame):
    """
    Count the number of in ice muons and get the most energetic one.

    Adds To Frame
    -------------
    NumInIceMuons : I3Double
        The number of in ice muons from the I3MCTree.

    TruthMuon : I3Particle
        The in ice muon with the largest energy. We take this muon as being
        the "truth".
    """

    tree = frame['I3MCTree']

    muons = []
    # Iterate over the in ice particles and find the muons.
#    for particle in tree.in_ice:
    for particle in tree:
        if particle.type_string in ('MuPlus', 'MuMinus') and particle.location_type_string == 'InIce':
            muons.append(particle)

    frame['NumInIceMuons'] = dataclasses.I3Double(len(muons))

    # Now find the in ice muon with the highest energy.
    max_energy_muon = muons[0]
    for muon in muons:
        if muon.energy > max_energy_muon.energy:
            max_energy_muon = muon

    frame['TruthMuon'] = max_energy_muon


def get_truth_endpoint(frame):
    """
    Get the endpoint of the truth muon track.
    Changed to calculating the endpoint directly instead of the now depreciated
    shift_along_track method.

    Adds to Frame
    -------------
    TruthEndpoint : I3Position
        The endpoint of the truth muon track.
    """

    muon = frame['TruthMuon']
#    truth_endpoint = muon.shift_along_track(muon.length)
    truth_endpoint = muon.pos+(muon.length*muon.dir)
    frame['TruthEndpoint'] = truth_endpoint


def count_hits(frame, pulses_name):
    """
    Count the number of pulses in the given pulse series that occur in specific
    regions of the detector.

    Parameters
    ---------
    pulses_name : str
        The key of the pulse series in the I3Frame.

    Adds To Frame
    -------------
    ICAnalysisHits : I3Double
        The number of hits within the IC analysis region.

    DCAnalysisHits : I3Double
        The number of hits within the DC analysis region.

    ICNHits : I3Double
        The number of hits outside the IC analysis region.

    DCNHits : I3Double
        The number of hits outside the DC analysis region.
    """

    # Strings in the analysis regions
    IC_strings = [26, 27, 37, 46, 45, 35, 17, 18, 19, 28, 38, 47, 56, 55, 54, 44, 34, 25]
    DC_strings = [81, 82, 83, 84, 85, 86]

    IC_analysis_hits = 0
    DC_analysis_hits = 0
    IC_nhits = 0
    DC_nhits = 0

    pulse_series = frame[pulses_name].apply(frame)

    # Iterato over the doms in the pulse series and count the hits in each
    # region
    for dom in pulse_series.keys():
        if dom.string in IC_strings and dom.om >= 40:
            IC_analysis_hits += 1
        if dom.string in DC_strings and dom.om >= 11:
            DC_analysis_hits += 1

        # We need to exclude hits on strings 36, 79, and 80
        if dom.string not in [36, 79, 80] + DC_strings + IC_strings:
            IC_nhits += 1
        if dom.string not in [36, 79, 80] + DC_strings:
            DC_nhits += 1

    frame['ICAnalysisHits'] = dataclasses.I3Double(IC_analysis_hits)
    frame['DCAnalysisHits'] = dataclasses.I3Double(DC_analysis_hits)
    frame['ICNHits'] = dataclasses.I3Double(IC_nhits)
    frame['DCNHits'] = dataclasses.I3Double(DC_nhits)


def reco_endpoint(frame, endpoint_fit):
    """
    Calculate the reconstructed endpoint of the event.
    Updated to calculate the endpoint from the position and direction instead of
    using the depreciated "shift_along_track" method.

    Parameters
    ----------
    endpoint_fit : str
        The fit used to calculate the endpoint.

    Adds To Frame
    -------------
    RecoEndpoint : I3Position
        The reconstructed endpoint of the event.
    """

    reco_fit = frame[endpoint_fit]
#    endpoint = reco_fit.shift_along_track(reco_fit.length)
    endpoint = reco_fit.pos+(reco_fit.length*reco_fit.dir)
    frame['RecoEndpoint'] = endpoint


def move_cut_variables(frame, direct_hits_name, fit_params_name):
    """
    Move NDirDoms calculated in direct_hits, the rlogl fit parameter for the
    provided fit, and the z coordinate of the reconstructed endpoint into the
    top level of the frame. This is done so cuts.py can access these values
    directly.

    We want to use the 'C' time window for direct hits, ie. [-15 ns, +75 ns].

    This function requires the finiteReco module be imported from icecube to
    access the fit parameter object in the frame.

    Parameters
    ----------
    direct_hits_name : str
        The OutputI3DirectHitsValuesBaseName.

    fit_params_name : str
        The fit parameters frame object name for the desired fit.

    Adds To Frame
    -------------
    NDirDoms : I3Double
    rlogl : I3Double
    RecoEndpointZ : I3Double
    """

    direct_hits = frame[direct_hits_name + 'C']
    frame['NDirDoms'] = dataclasses.I3Double(direct_hits.n_dir_doms)

    fit_params = frame[fit_params_name]
    frame['rlogl'] = dataclasses.I3Double(fit_params.rlogl)

    reco_z = frame['RecoEndpoint'].z
    frame['RecoEndpointZ'] = dataclasses.I3Double(reco_z)
