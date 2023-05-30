#
# Module to output events in HDF5 format
#
from tables import open_file
import pickle
from icecube import icetray

from event import *

import numpy as np
from icecube.phys_services import I3Calculator as calc
from icecube.dataclasses import I3Constants

def GetMuonEnergy(MMCTrackList) :
	nmuons  = 0
	MaxMuonEnergy = 0.0
	maxemuon = None
	for track in MMCTrackList :
		xc = track.xc
		yc = track.yc
		zc = track.zc
		ec = track.Ec

		if np.sqrt(xc*xc+yc*yc+zc*zc) < 500. :
			nmuons += 1
			if ec > MaxMuonEnergy:
				MaxMuonEnergy = ec
				maxemuon = track.GetI3Particle()

	return nmuons, maxemuon

class EventWriterRecoPulses(icetray.I3Module):
	def __init__(self, context):
		icetray.I3Module.__init__(self, context)
		self.AddParameter('FileName','Name of the file to write out','out.h5')
		self.AddParameter('pulseseries','Nulse series to get DOM info from','I3Photons')

	def Configure(self):
		self.h5file = open_file(self.GetParameter('FileName'), mode="w", title="DOM Calibration HDF5 File")
		self.doms   = self.h5file.create_table('/', 'doms', recodata, "Hit DOMs")
		self.pulseseriesname = self.GetParameter('pulseseries')
		
	def Physics(self, frame):
		if not frame.Has(self.pulseseriesname) :
			return

		pulseseries = frame[self.pulseseriesname].apply(frame)
	
		nmuons, muon = GetMuonEnergy(frame['MMCTrackList'])

		if nmuons != 1 :
			return

		dom = self.doms.row
		domsused = frame['I3Geometry'].omgeo

		n_ice_group = I3Constants.n_ice_group
    		n_ice_phase = I3Constants.n_ice_phase

		for domkey in pulseseries.keys() :
			dom_position = domsused[domkey].position
			dist = calc.cherenkov_distance(muon, dom_position, n_ice_group, n_ice_phase) 
			totalcharge = 0.0
			for pulse in pulseseries[domkey] :
				totalcharge += pulse.charge
			dom['dist'] = dist
			dom['charge'] = totalcharge
			dom.append()

		self.PushFrame(frame)                 # push the frame
		return

	def __del__(self):
		self.h5file.close()
