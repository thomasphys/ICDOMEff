"""
Functions used for analyzing the per-dom (as opposed to per-event) data.
"""

from __future__ import print_function, division  # 2to3

import math

from icecube import dataclasses, finiteReco
from icecube.phys_services import I3Calculator as calc
from icecube.dataclasses import I3Constants


#def om_partition(frame, output_name, options):
#    """
#    Partition the pulses.
#
#    Parameters
#    ----------
#    output_name : str
#    options : dict[str]
#
#    Adds To Frame
#    -------------
#    output_name.format(partition) : I3RecoPulseSeriesMap
#
#    """
#
#    # Initialize the RecoPulseSeriesMaps
#    for partition in range(options['partitions']):
#        key = output_name.format(partition)
#        frame[key] = dataclasses.I3RecoPulseSeriesMap()
#
#    # Get the pulse series
#    pulse_series = frame[options['pulses_name']].apply(frame)
#
#    for dom, pulse_vector in pulse_series.items():
#
#        # Find out which partition the pulse_vector should go in
#        partition_num = (dom.string + dom.om) % options['partitions']
#
#        # Put it in every partition except the one it is in.
#        for partition in range(options['partitions']):
#            if partition != partition_num:
#                key = output_name.format(partition)
#                frame[key][dom] = pulse_vector


def dom_data(frame, reco_fit, options):
    """
    Analyze and save the per-dom data using the provided fit.

    Parameters
    ----------
    reco_fit : str
        The key of the MPEFit reconstruction

    options : dict[str]

    Adds To Frame
    -------------
    TotalCharge : I3VectorDouble
    String : I3VectorDouble
    OM : I3VectorDouble
    DistAboveEndpoint : I3VectorDouble
    ImpactAngle : I3VectorDouble
    RecoDistance : I3VectorDouble

    Returns
    -------
    bool
        Whether any per-dom data was added to the frame. If no data was added,
        return False, because this frame contains no pertinent information.
        Otherwise return True.

    """

    n_ice_group = I3Constants.n_ice_group
    n_ice_phase = I3Constants.n_ice_phase

    IC_strings = [26, 27, 37, 46, 45, 35, 17, 18, 19, 28, 38, 47, 56, 55, 54, 44, 34, 25]
    DC_strings = [81, 82, 83, 84, 85, 86]

    reco_endpoint = frame['RecoEndpoint']

    # Get the pulse series
    pulse_series = frame[options['pulses_name']].apply(frame)

    mcpulses = None
    if frame.Has("I3MCPulseSeriesMap") :
        mcpulses = frame["I3MCPulseSeriesMap"]
        frame['DOM_MCPulses'] = dataclasses.I3VectorDouble()
    # Initialize the vectors containing the per DOM data
    frame['DOM_TotalCharge'] = dataclasses.I3VectorDouble()
    frame['DOM_TotalCharge_300ns'] = dataclasses.I3VectorDouble()
    frame['DOM_String'] = dataclasses.I3VectorDouble()
    frame['DOM_OM'] = dataclasses.I3VectorDouble()
    frame['DOM_DistAboveEndpoint'] = dataclasses.I3VectorDouble()
    frame['DOM_ImpactAngle'] = dataclasses.I3VectorDouble()
    frame['DOM_RecoDistance'] = dataclasses.I3VectorDouble()
    frame['DOM_TruthDistance'] = dataclasses.I3VectorDouble()
    frame['DOM_MinTimeResidual'] = dataclasses.I3VectorDouble()

    dom_geo = frame['I3Geometry'].omgeo.items()

    # Find all doms above the reconstructed z coord of endpoint and
    # within the specified distance interval of the track

    maxdist = 0.0
    closez = 0.0
    endpoint = 0.0
    alldom = 0.0
    wrongstring = 0.0

    for dom, geo in dom_geo:  # (OMKey, I3OMGeo)

        alldom += 1.0

        # We want to get DOMs that are in the IC/DC strings and below the dust
        # layer (40 and below for IC, 11 and below for DC).
        if (dom.string in IC_strings and dom.om >= 40) or (dom.string in DC_strings and dom.om >= 11):

            dom_position = geo.position
            #partition_num = (dom.string + dom.om) % options['partitions']
            #mpe = frame[reco_fit.format(partition_num)]  # MPEFit0...4
            #mpe = frame['MPEFit']
            # Get the reconstructed track using the frame key passed to this function
            reco_track = frame[reco_fit]
            # Find cherenkov distance from track to DOM
            reco_dist = calc.cherenkov_distance(reco_track, dom_position, n_ice_group, n_ice_phase)

            truth_dist = 50000.0
            true_track = None
            if frame.Has('MMCTrackList') :
               mctracks = frame['MMCTrackList']
               for track in mctracks:
                 temp_dist = calc.cherenkov_distance(track.GetI3Particle(), dom_position, n_ice_group, n_ice_phase)		
                 truth_dist = min(truth_dist,temp_dist)
                 if truth_dist == temp_dist :
                     true_track = track.GetI3Particle()
            

            # If the DOM is outside the maximum distance from the track or the track's position of
            # closest approach to the DOM is above the level of the DOM then skip this DOM
            if reco_dist > options['max_dist']:
		maxdist += 1.0
                continue
            # Keep if track is below DOM
            clos_app_pos = calc.closest_approach_position(reco_track, dom_position)
            if clos_app_pos.z > dom_position.z:
		closez += 1.0
                continue

            # Calculate the point at which direct Cherenkov light hitting the DOM would have been emitted
            # at and require that this be above the end point of the track
            cherenkov_pos = calc.cherenkov_position(reco_track, dom_position, n_ice_group, n_ice_phase)
            dist_above_endpoint = calc.distance_along_track(reco_track, reco_endpoint) - calc.distance_along_track(reco_track, cherenkov_pos)
            if dist_above_endpoint <= 0:
		endpoint += 1.0
                continue

            # We have survived to this point in the loop so now we need to store the data for this
            # DOM in the arrays we created
            frame['DOM_RecoDistance'].append(reco_dist)
            frame['DOM_TruthDistance'].append(truth_dist)
            frame['DOM_DistAboveEndpoint'].append(dist_above_endpoint)
            frame['DOM_String'].append(dom.string)
            frame['DOM_OM'].append(dom.om)

            # Calculate the angle at which direct Cherenkov light would impact the DOM
            perp_position = dataclasses.I3Position(dom_position.x, dom_position.y, clos_app_pos.z)
            delta = perp_position - clos_app_pos
            impact_param = delta.magnitude
            impact_angle = math.asin(impact_param / calc.closest_approach_distance(reco_track, dom_position))
            frame['DOM_ImpactAngle'].append(impact_angle)

            # Calculate the total charge and the time residual for the DOM
            total_charge = 0.0
            totalmc_charge = 0.0
            total_charge_300ns = 0.0
            min_time_residual = 1000.
           
            # If there are pulses, sum the charge of the ones with a time residual less than 1000 ns.
            if dom in pulse_series.keys():
                for pulse in pulse_series[dom]:
                    time_res = calc.time_residual(reco_track, dom_position, pulse.time, n_ice_group, n_ice_phase)
                    min_time_residual = min(min_time_residual,time_res)
                    if time_res < 1000.:
                        total_charge += pulse.charge
                    if time_res > -100. and time_res < 300. :
                        total_charge_300ns += pulse.charge
            if mcpulses :
	        if dom in mcpulses.keys() :
                     for pulse in mcpulses[dom]:
   			time_res = calc.time_residual(true_track, dom_position, pulse.time, n_ice_group, n_ice_phase)
                        if time_res < 1000.:
                             if pulse.source == 10 :
                                 totalmc_charge += pulse.charge
		    
	    #print("appending vectors")
            frame['DOM_TotalCharge'].append(total_charge)
            if frame.Has('DOM_MCPulses'):
                frame['DOM_MCPulses'].append(totalmc_charge)
            frame['DOM_TotalCharge_300ns'].append(total_charge_300ns) 
            frame['DOM_MinTimeResidual'].append(min_time_residual)
        else :
	    wrongstring += 1.0
    # After all that, if none of the DOMs made it through, get rid of this
    # frame.
    #print("maxdist = %f closez = %f and endpoint = %f wrongstring = %f total = %f" % (maxdist,closez,endpoint,wrongstring,alldom) )
    return len(frame['DOM_TotalCharge']) != 0
