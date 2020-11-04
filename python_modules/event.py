#
# HDF5 Table for DOMs which have been hit
#
from tables import *

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
	corsika                = Corsika()
	startTime              = Time()
	endTime                = Time()
	totalCharge            = Float64Col()
	borderDistance         = Float64Col()
	truthBorderDistance    = Float64Col()
	directHits             = UInt32Col()
	stopLikeRatio          = Float64Col()

class DOM(IsDescription):
	eventId           = UInt32Col()   # ID number of the event this DOM belongs to
	string            = UInt16Col()       # String containing the DOM
	om                = UInt16Col()       # Opitcal module number of the DOM
	recoDist          = Float64Col()
	distAboveEndpoint = Float64Col()  # Distance of the DOM above the endpoint of the track
	impactAngle       = Float64Col()
	totalCharge       = Float64Col()  # Total charge seen by the DOM
	totalCharge_300ns = Float64Col()  # Total charge seen by the DOM within a time residual of -100ns to 300 ns.
	minTimeResidual   = Float64Col()  # Minimum time residual for all pulses
