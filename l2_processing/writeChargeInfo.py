#
# Module to output events in HDF5 format
#
from tables import open_file

from icecube import icetray

from event import *

import numpy as np

from icecube import dataclasses, finiteReco
from icecube.phys_services import I3Calculator as calc
from icecube.dataclasses import I3Constants
from icecube.weighting.fluxes import  GaisserH4a, GaisserH3a, GaisserH4a_IT, GaisserHillas, Hoerandel, Hoerandel5, Hoerandel_IT

class EventWriter(icetray.I3Module):
	def __init__(self, context):
		icetray.I3Module.__init__(self, context)
		self.AddParameter('FileName','Name of the file to write out','out.h5')
		self.eventId=0

	def Configure(self):
		self.h5file = open_file(self.GetParameter('FileName'), mode="w", title="DOM Calibration HDF5 File")
		# Create the table to store all the event objects
		self.events = self.h5file.create_table('/', 'events', ChargeInfo, "Calibration Events")
		# Create the table to store all the DOM hits
		# Create an index on the event ID
		self.events.cols.eventId.create_index()

	def Physics(self, frame):
		# Create an event entry
		event=self.events.row
		# Get the current event ID number
		event['eventId'] = frame['I3EventHeader'].event_id
		n_ice_group = I3Constants.n_ice_group
    		n_ice_phase = I3Constants.n_ice_phase
		dom_geo = frame['I3Geometry'].omgeo
		flux = GaisserH4a()

		# Loop over the DOM hits and add them to the DOM table
		pulseseriesmap = frame['SRTInIcePulsesDOMEff'].apply(frame)
		totalcharge = 0.0
		for dom in pulseseriesmap.keys():
			for pulse in pulseseriesmap[dom] :
				time_res = calc.time_residual(frame['SplineMPE'], dom_geo[dom].position, pulse.time, n_ice_group, n_ice_phase)
				if time_res > -200. and time_res < 1000. :
					totalcharge += pulse.charge
		event['zenith'] = frame['SplineMPE'].dir.zenith	
		event['totalcharge'] = totalcharge

		CorsikaWeight = frame['CorsikaWeightMap']
                pflux = flux(CorsikaWeight['PrimaryEnergy'],CorsikaWeight['PrimaryType'])
                energy_integral = CorsikaWeight['EnergyPrimaryMax']**(CorsikaWeight['PrimarySpectralIndex']+1.0)
                #print("energy_integral 1 = "+str(energy_integral))
                energy_integral = energy_integral - CorsikaWeight['EnergyPrimaryMin']**(CorsikaWeight['PrimarySpectralIndex']+1.0)
                #print("energy_integral 2 = "+str(energy_integral))
                energy_integral = energy_integral / (CorsikaWeight['PrimarySpectralIndex']+1.0)
                #print("energy_integral 3 = "+str(energy_integral))
                energy_weight = CorsikaWeight['PrimaryEnergy']**CorsikaWeight['PrimarySpectralIndex']
                #print('energy_weight1 = '+str(energy_weight))
                energy_weight = pflux*energy_integral/energy_weight*CorsikaWeight['AreaSum']
                event['weight'] = energy_weight/float(CorsikaWeight['NEvents'])

		#print('zenith = '+str(event['zenith'])+" totalcharge = "+str(event['totalcharge']))
		event.append()
		self.PushFrame(frame)                 # push the frame
		return

	def __del__(self):
		#self.events.flush()
		#self.doms.flush()
		#self.h5file.flush()
		self.h5file.close()
