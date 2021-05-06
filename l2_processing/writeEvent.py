#
# Module to output events in HDF5 format
#
from tables import open_file

from icecube import icetray

from event import *

def filtertriggerresult(state) :
	if state :
		return 1
	else :
		return 0

class EventWriter(icetray.I3Module):
	def __init__(self, context):
		icetray.I3Module.__init__(self, context)
		self.AddParameter('FileName','Name of the file to write out','out.h5')
		self.eventId=0

	def Configure(self):
		self.h5file = open_file(self.GetParameter('FileName'), mode="w", title="DOM Calibration HDF5 File")
		# Create the table to store all the event objects
		self.events = self.h5file.create_table('/', 'events', Event, "Calibration Events")
		# Create the table to store all the DOM hits
		self.doms   = self.h5file.create_table('/', 'doms', DOM, "Hit DOMs")
		# Create the table to store run information.
		self.run = self.h5file.create_table('/','runinfo',Run, "Run Info")
		# Create an index on the event ID
		self.doms.cols.eventId.create_index()
		self.events.cols.eventId.create_index()

	def Physics(self, frame):
		# Create an event entry
		event=self.events.row
		# Get the current event ID number
		event['eventId']=self.eventId
		event['eventId'] = frame['I3EventHeader'].event_id
		event['reco/time']   = frame['SplineMPE'].time
		event['reco/energy'] = frame['SplineMPE'].energy
		event['reco/speed']  = frame['SplineMPE'].speed
		event['reco/length'] = frame['SplineMPE'].length
		event['reco/pos/x'] = frame['SplineMPE'].pos.x
		event['reco/pos/y'] = frame['SplineMPE'].pos.y
		event['reco/pos/z'] = frame['SplineMPE'].pos.z
		event['mpe/time']   = frame['MPEFitDOMeff'].time
		event['mpe/energy'] = frame['MPEFitDOMeff'].energy
		event['mpe/speed']  = frame['MPEFitDOMeff'].speed
		event['mpe/length'] = frame['MPEFitDOMeff'].length
		event['mpe/pos/x'] = frame['MPEFitDOMeff'].pos.x
		event['mpe/pos/y'] = frame['MPEFitDOMeff'].pos.y
		event['mpe/pos/z'] = frame['MPEFitDOMeff'].pos.z
		event['spe/time']   = frame['SPEFit2DOMeff'].time
		event['spe/energy'] = frame['SPEFit2DOMeff'].energy
		event['spe/speed']  = frame['SPEFit2DOMeff'].speed
		event['spe/length'] = frame['SPEFit2DOMeff'].length
		event['spe/pos/x'] = frame['SPEFit2DOMeff'].pos.x
		event['spe/pos/y'] = frame['SPEFit2DOMeff'].pos.y
		event['spe/pos/z'] = frame['SPEFit2DOMeff'].pos.z
		event['line/time']   = frame['LineFitDOMeff'].time
		event['line/energy'] = frame['LineFitDOMeff'].energy
		event['line/speed']  = frame['LineFitDOMeff'].speed
		event['line/length'] = frame['LineFitDOMeff'].length
		event['line/pos/x'] = frame['LineFitDOMeff'].pos.x
		event['line/pos/y'] = frame['LineFitDOMeff'].pos.y
		event['line/pos/z'] = frame['LineFitDOMeff'].pos.z
		event['recoEndPoint/x'] = frame['RecoEndpoint'].x
		event['recoEndPoint/y'] = frame['RecoEndpoint'].y
		event['recoEndPoint/z'] = frame['RecoEndpoint'].z
		#event['truthEndPoint/x'] = frame['TruthEndpoint'].x
		#event['truthEndPoint/y'] = frame['TruthEndpoint'].y
		#event['truthEndPoint/z'] = frame['TruthEndpoint'].z
		event['reco/dir/zenith']  = frame['SplineMPE'].dir.zenith
		event['reco/dir/azimuth'] = frame['SplineMPE'].dir.azimuth
		event['mpe/dir/zenith']  = frame['MPEFitDOMeff'].dir.zenith
		event['mpe/dir/azimuth'] = frame['MPEFitDOMeff'].dir.azimuth
		event['spe/dir/zenith']  = frame['SPEFit2DOMeff'].dir.zenith
		event['spe/dir/azimuth'] = frame['SPEFit2DOMeff'].dir.azimuth
		event['line/dir/zenith']  = frame['LineFitDOMeff'].dir.zenith
		event['line/dir/azimuth'] = frame['LineFitDOMeff'].dir.azimuth
		event['startTime/tag'] = frame['I3EventHeader'].start_time
		event['endTime/tag'] = frame['I3EventHeader'].end_time
		event['dcHitsIn'] = frame['DCAnalysisHits'].value    
		event['dcHitsOut'] = frame['DCNHits'].value    
		event['icHitsIn'] = frame['ICAnalysisHits'].value    
		event['icHitsOut'] = frame['ICNHits'].value   
		event['borderDistance'] = frame['DistToBorder'].value
		#event['truthBorderDistance'] = frame['TruthDistToBorder'].value
		event['recoLogL'] = frame['rlogl'].value
		event['directHits'] = frame['NDirDoms'].value
		event['totalCharge'] = frame['EventCharge'].value
		event['stopLikeRatio'] = frame['FiniteRecoLLHRatio'].value

		filter_mask = []
		if frame.Has('FilterMask') :
			filter_mask = frame['FilterMask']
		elif frame.Has('QFilterMask') :
			filter_mask = frame['QFilterMask']
		elif frame.Has('QFilterMask_NullSplit0') :
			filter_mask = frame['QFilterMask_NullSplit0']

		event['filterMask/CascadeFilter_13']                   = filtertriggerresult(filter_mask['CascadeFilter_13'].condition_passed)
		event['filterMask/DeepCoreFilter_13']                  = filtertriggerresult(filter_mask['DeepCoreFilter_13'].condition_passed)
	#	vent['filterMask/DeepCoreFilter_TwoLayerExp_13']      = filtertriggerresult(filter_mask['DeepCoreFilter_TwoLayerExp_13'].condition_passed)
		#event['filterMask/EHEFilter_13']                       = filtertriggerresult(filter_mask['EHEFilter_13'].condition_passed)
		event['filterMask/FSSCandidate_13']                    = filtertriggerresult(filter_mask['FSSCandidate_13'].condition_passed)
		event['filterMask/FSSFilter_13']                       = filtertriggerresult(filter_mask['FSSFilter_13'].condition_passed)
		event['filterMask/FilterMinBias_13']                   = filtertriggerresult(filter_mask['FilterMinBias_13'].condition_passed)
		event['filterMask/FixedRateFilter_13']                 = filtertriggerresult(filter_mask['FixedRateFilter_13'].condition_passed)
		#event['filterMask/GCFilter_13']                        = filtertriggerresult(filter_mask['GCFilter_13'].condition_passed)
		event['filterMask/I3DAQDecodeException']               = filtertriggerresult(filter_mask['I3DAQDecodeException'].condition_passed)
		event['filterMask/IceTopSTA3_13']                      = filtertriggerresult(filter_mask['IceTopSTA3_13'].condition_passed)
		event['filterMask/IceTopSTA5_13']                      = filtertriggerresult(filter_mask['IceTopSTA5_13'].condition_passed)
		event['filterMask/IceTop_InFill_STA3_13']              = filtertriggerresult(filter_mask['IceTop_InFill_STA3_13'].condition_passed)
		event['filterMask/InIceSMT_IceTopCoincidence_13']      = filtertriggerresult(filter_mask['InIceSMT_IceTopCoincidence_13'].condition_passed)
		#event['filterMask/LID']                                = filtertriggerresult(filter_mask['LID'].condition_passed)
		event['filterMask/LowUp_13']                           = filtertriggerresult(filter_mask['LowUp_13'].condition_passed)
		event['filterMask/MoonFilter_13']                      = filtertriggerresult(filter_mask['MoonFilter_13'].condition_passed)
		event['filterMask/MuonFilter_13']                      = filtertriggerresult(filter_mask['MuonFilter_13'].condition_passed)
		#event['filterMask/OFUFilter_14']                       = filtertriggerresult(filter_mask['OFUFilter_14'].condition_passed)
		#event['filterMask/OnlineL2Filter_14']                  = filtertriggerresult(filter_mask['OnlineL2Filter_14'].condition_passed)
		#event['filterMask/SDST_FilterMinBias_13']              = filtertriggerresult(filter_mask['SDST_FilterMinBias_13'].condition_passed)
		event['filterMask/SDST_IceTopSTA3_13']                 = filtertriggerresult(filter_mask['SDST_IceTopSTA3_13'].condition_passed)
		event['filterMask/SDST_IceTop_InFill_STA3_13']         = filtertriggerresult(filter_mask['SDST_IceTop_InFill_STA3_13'].condition_passed)
		event['filterMask/SDST_InIceSMT_IceTopCoincidence_13'] = filtertriggerresult(filter_mask['SDST_InIceSMT_IceTopCoincidence_13'].condition_passed)
		event['filterMask/SlopFilter_13']                      = filtertriggerresult(filter_mask['SlopFilter_13'].condition_passed)
		event['filterMask/SunFilter_13']                       = filtertriggerresult(filter_mask['SunFilter_13'].condition_passed)
		event['filterMask/VEF_13']                             = filtertriggerresult(filter_mask['VEF_13'].condition_passed)

		event['triggerMask/InIceSMT']                          = filtertriggerresult(frame['InIceSMTTriggered'])
		event['triggerMask/IceTopSMT']                         = filtertriggerresult(frame['IceTopSMTTriggered'])
		event['triggerMask/InIceString']                       = filtertriggerresult(frame['InIceStringTriggered'])
		event['triggerMask/PhysMinBias']                       = filtertriggerresult(frame['PhysMinBiasTriggered'])
		event['triggerMask/DeepCoreSMT']                       = filtertriggerresult(frame['DeepCoreSMTTriggered'])

		#if frame.Has('SRTInIcePulses') :
	#		event['passl2cuts'] = 1
	#	else :
	#		event['passl2cuts'] = 0
		

		if frame.Has('CorsikaWeightMap'):
    
			event['corsika/primaryEnergy'] = frame['CorsikaWeightMap']['PrimaryEnergy']
			event['corsika/primaryType'] = frame['CorsikaWeightMap']['PrimaryType'] 
			event['corsika/primarySpectralIndex'] = frame['CorsikaWeightMap']['PrimarySpectralIndex']
			event['corsika/energyPrimaryMin'] = frame['CorsikaWeightMap']['EnergyPrimaryMin']
			event['corsika/energyPrimaryMax'] = frame['CorsikaWeightMap']['EnergyPrimaryMax']  
			event['corsika/areaSum'] = frame['CorsikaWeightMap']['AreaSum']
			event['corsika/nEvents'] = frame['CorsikaWeightMap']['NEvents']

		if self.eventId == 0 :
			runinfo = self.run.row
			if frame.Has('CorsikaWeightMap'):
				runinfo['nevents'] = frame['CorsikaWeightMap']['NEvents']
			runinfo['runId'] = frame['I3EventHeader'].run_id
			runinfo['subRunId'] = frame['I3EventHeader'].sub_run_id
			runinfo.append()

		event.append()

		# Loop over the DOM hits and add them to the DOM table
		for i in range(len(frame['DOM_TotalCharge'])):
			dom=self.doms.row
			dom['eventId']           = self.eventId
			dom['string']            = frame['DOM_String'][i]
			dom['om']                = frame['DOM_OM'][i]
			dom['totalCharge']       = frame['DOM_TotalCharge'][i]
			dom['totalCharge_300ns'] = frame['DOM_TotalCharge_300ns'][i]
			dom['distAboveEndpoint'] = frame['DOM_DistAboveEndpoint'][i]
			dom['impactAngle']       = frame['DOM_ImpactAngle'][i]
			dom['minTimeResidual']   = frame['DOM_MinTimeResidual'][i]
			dom['recoDist']          = frame['DOM_RecoDistance'][i]
			if frame.Has('MMCTrackList') :
				dom['truthDist']         = frame['DOM_TruthDistance'][i]
                        if frame.Has('DOM_MCPulses') :
                                dom['totalChargeMC'] = frame['DOM_MCPulses'][i]
			dom.append()
		# Increment the event ID
		self.eventId+=1
		self.PushFrame(frame)                 # push the frame
		return

	def __del__(self):
		#self.events.flush()
		#self.doms.flush()
		#self.h5file.flush()
		self.h5file.close()
