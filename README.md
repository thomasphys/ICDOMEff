This repository contains all the coded needed to perform the UAlberta DOM efficiency analysis.

==================================================================================================
MC Production Files:

Corsika_filter_muons.py [David Gillcrist]
	Generates the initial corsikaevents using the inputs from config.json.

photon_syst.py [David Gillcrist]
	Performs the photon propagation with ClSim.

===================================================================================================
Processing Level2 Files and data files:

I3Script: process_splineMPE_2015.py [Combined and updated Nick Kulzac's scripts by Thomas McElroy]
Arguments:
	-d
		Directory of data files
        -g
		Geometry file, the default is ${I3_DATA}/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE.i3.gz
  	-o
		Name of output file.
       	-z
		Range of muon Zeniths in degrees.
  	-q
		Range of muon Energies.
  	-s
		True or False, Is file simulation?
  	-n	
		Number of events to process, default is -1 which will process all events.
  	-t 
		File type of data files, default is i3.bz2
 	-r 	
		Number to identify target file, for data this is subrun and for MC it is run number.
  	-p 
		Name of new pulse list, the default is SRTInIcePulsesDOMEff
	-m 
		Maximum distance to DOM to save, the default is 140.0


Example:

python process_splineMPE_2015.py  -g /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/../data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE.i3.gz -d /data/sim/IceCube/2016/generated/CORSIKA-in-ice/21269/IC86_2016_spe_templates_DOM_oversize1/level2/redo/eff090/0000000-0000999/Level2_eff090_IC86-2016_corsika_21269_ -r 0 -t .i3.zst -o -1 -s True


Dependancies:

general.py
	Contains functions used to pull data out of files and perform selection cuts on data. 

geoanalysis.py
	contains geometry functions like computing distance to detector boarder. 

domanalysis.py
	pulls the data out of 

event.py
	Data classes for making HDF5 files.

writeEvent.py
	Module to write data to HDF5 file.

===================================================================================================

DOM Efficiency Analysis Code.

CompilePlotData.py
	This file reads in the HDF5 files from the processing file and extracts the information needed for assessing the DOM efficiency and computing systematics.

Dependancies:

event.py
