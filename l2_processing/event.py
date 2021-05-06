#
# HDF5 Table for DOMs which have been hit
#
from tables import *

class data(IsDescription):
	value = Float64Col()

class State(IsDescription):
	data           = StringCol(500)
	eff            = StringCol(500)
	flux           = StringCol(500)
	zenithmin      = Float64Col()
	zenithmax      = Float64Col()
	energymin      = Float64Col()
	impactanglemin = Float64Col()
	impactanglemax = Float64Col()
	trackendpoint  = Float64Col()
	cherdistmin    = Float64Col()
	cherdistmax    = Float64Col()
	binwidth       = Float64Col()

class DataPoint(IsDescription):
	meancharge    = Float64Col()
	sigmacharge   = Float64Col()
	meandistance  = Float64Col()
	sigmadistance = Float64Col()

class Data(IsDescription) :
	charge = Float64Col()
	dist = Float64Col()

class MC(IsDescription) :
	charge = Float64Col()
	dist = Float64Col()
	weight = Float64Col()

class Flux(IsDescription):
	totalcharge = Float64Col()
	energy      = Float64Col()
	zenith      = Float64Col()
	weight      = Float64Col()

class Position(IsDescription):
	x = Float64Col()
	y = Float64Col()
	z = Float64Col()

class Direction(IsDescription):
	zenith  = Float64Col()
	azimuth = Float64Col()

class Particle(IsDescription):
	id     = UInt32Col()  # Particle ID
	pos    = Position()   # Starting vertex of the particle
	dir    = Direction()  # Direction that the particle came from
	time   = Float64Col()
	energy = Float64Col() # Energy of the particle
	speed  = Float64Col()
	length = Float64Col()

class Run(IsDescription):
	nevents   = UInt32Col()
	startTime = Float64Col()
	endTime   = Float64Col()
	runId     = UInt32Col()
	subRunId  = UInt32Col()

class Corsika(IsDescription):
	primaryEnergy        = Float64Col()
	primaryType          = UInt32Col()
	primarySpectralIndex = Float64Col()
	energyPrimaryMin     = Float64Col()
	energyPrimaryMax     = Float64Col()
	areaSum              = Float64Col()
	nEvents              = Float64Col()

class Time(IsDescription) :
	tag = StringCol(37)

class FilterMask(IsDescription) :
	CascadeFilter_13                   = UInt32Col()
	DeepCoreFilter_13                  = UInt32Col()
	DeepCoreFilter_TwoLayerExp_13      = UInt32Col()
	EHEFilter_13                       = UInt32Col()
	FSSCandidate_13                    = UInt32Col()
	FSSFilter_13                       = UInt32Col()
	FilterMinBias_13                   = UInt32Col()
	FixedRateFilter_13                 = UInt32Col()
	GCFilter_13                        = UInt32Col()
	I3DAQDecodeException               = UInt32Col()
	IceTopSTA3_13                      = UInt32Col()
	IceTopSTA5_13                      = UInt32Col()
	IceTop_InFill_STA3_13              = UInt32Col()
	InIceSMT_IceTopCoincidence_13      = UInt32Col()
	LID                                = UInt32Col()
	LowUp_13                           = UInt32Col()
	MoonFilter_13                      = UInt32Col()
	MuonFilter_13                      = UInt32Col()
	OFUFilter_14                       = UInt32Col()
	OnlineL2Filter_14                  = UInt32Col()
	SDST_FilterMinBias_13              = UInt32Col()
	SDST_IceTopSTA3_13                 = UInt32Col()
	SDST_IceTop_InFill_STA3_13         = UInt32Col()
	SDST_InIceSMT_IceTopCoincidence_13 = UInt32Col()
	SlopFilter_13                      = UInt32Col()
	SunFilter_13                       = UInt32Col()
	VEF_13                             = UInt32Col()

class TriggerMask(IsDescription) :
	InIceSMT    = UInt32Col()
	IceTopSMT   = UInt32Col()
	InIceString = UInt32Col()
	PhysMinBias = UInt32Col()
	DeepCoreSMT = UInt32Col()

class Event(IsDescription):
	eventId                = UInt32Col()    # ID number of the event
	dcHitsIn               = Float64Col()    # Number of hits inside the deep core analysis region
	dcHitsOut              = Float64Col()    # Number of hits outside the deep core analysis region
	icHitsIn               = Float64Col()    # Number of hits inside the IceCube analysis region
	icHitsOut              = Float64Col()    # Number of hits outside the IceCube analysis region
	recoEndPoint           = Position()  # Coordinates of the reconstructed track endpoint
	truthEndPoint          = Position()
	recoLogL               = Float64Col()
	firstHit               = UInt32Col()    # Index of the first DOM hit for this event
	nHits                  = UInt32Col()    # Number of DOM hits in this event
	primary                = Particle()     # Polyplopia primary particle
	reco                   = Particle()     # Reconstructed particle track
	mpe					   = Particle()
	spe                    = Particle()
	line				   = Particle()
	corsika                = Corsika()
	startTime              = Time()
	endTime                = Time()
	totalCharge            = Float64Col()
	borderDistance         = Float64Col()
	truthBorderDistance    = Float64Col()
	directHits             = UInt32Col()
	stopLikeRatio          = Float64Col()
	filterMask             = FilterMask()
	triggerMask            = TriggerMask()

class DOM(IsDescription):
	eventId           = UInt32Col()   # ID number of the event this DOM belongs to
	string            = UInt16Col()       # String containing the DOM
	om                = UInt16Col()       # Opitcal module number of the DOM
	recoDist          = Float64Col()
	truthDist         = Float64Col()
	distAboveEndpoint = Float64Col()  # Distance of the DOM above the endpoint of the track
	impactAngle       = Float64Col()
	totalCharge       = Float64Col()  # Total charge seen by the DOM
        totalChargeMC     = Float64Col()
	totalCharge_300ns = Float64Col()  # Total charge seen by the DOM within a time residual of -100ns to 300 ns.
	minTimeResidual   = Float64Col()  # Minimum time residual for all pulses
