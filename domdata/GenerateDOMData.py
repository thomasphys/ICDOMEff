import os, sys
from tables import open_file
import ROOT
import numpy as np
import argparse
from tables import open_file
from icecube.weighting.fluxes import  GaisserH4a, GaisserH3a, GaisserH4a_IT, GaisserHillas, Hoerandel, Hoerandel5, Hoerandel_IT
from icecube.weighting import weighting
from event import *
from I3Tray import OMKey
from array import array

weights_ic = []
weights_dc = []
reconstructedE_ic = []
reconstructedE_dc = []
zenith_ic = []
zenith_dc = []
EnergyTruth_ic = []
EnergyTruth_dc = []
ZenithTruth_ic = []
ZenithTruth_dc = []
totalcharge_ic = []
totalcharge_dc = []
recoEndpoint_ic = []
recoEndpoint_dc = []
borderDistance_ic = []
borderDistance_dc = []
stopLikeRatio_ic = []
stopLikeRatio_dc = []
recoLogL_ic = []
recoLogL_dc = []
directHits_ic = []
directHits_dc = []
HitsOut_ic = []
HitsOut_dc = [] 

ImpactAll_ic = ROOT.TH1F("ImpactAll_IC","",1000,-1.0,1.0)
ImpactAll_dc = ROOT.TH1F("ImpactAll_DC","",1000,-1.0,1.0)
Impact_seeMPE_ic = ROOT.TH1F("ImpactseeMPE_IC","",1000,-1.0,1.0)
Impact_seeMPE_dc = ROOT.TH1F("ImpactseeMPE_DC","",1000,-1.0,1.0)

Impact_vs_Zenith_ic = ROOT.TH2F("Impact_vs_Zenith_IC","",200,-1.0,1.0,200,-1.0,1.0)
Impact_vs_Zenith_dc = ROOT.TH2F("Impact_vs_Zenith_DC","",200,-1.0,1.0,200,-1.0,1.0)

TotalCharge_vs_Zenith_ic = ROOT.TH2F("TotalCharge_vs_Zenith_IC","",200,-1.0,1.0,1000,0.0,3000.0)
TotalCharge_vs_Zenith_dc = ROOT.TH2F("TotalCharge_vs_Zenith_DC","",200,-1.0,1.0,1000,0.0,3000.0)

binneddistance_dc = np.zeros(1,dtype=float)
binneddistanceerror_dc = np.zeros(1,dtype=float)
binnedcharge_dc = np.zeros(1,dtype=float)
binnedchargeerror_dc = np.zeros(1,dtype=float)
binnedcharge300_dc = np.zeros(1,dtype=float)
binnedcharge300error_dc = np.zeros(1,dtype=float)
binneddistance_ic = np.zeros(1,dtype=float)
binneddistanceerror_ic = np.zeros(1,dtype=float)
binnedcharge_ic = np.zeros(1,dtype=float)
binnedchargeerror_ic = np.zeros(1,dtype=float)
binnedcharge300_ic = np.zeros(1,dtype=float)
binnedcharge300error_ic = np.zeros(1,dtype=float)

def find_error_bootstrap(values,weights):
# this needs to be eddited, I do not thing this was programmed corerctly.

	total = len(values)
	sum_weights = sum([weights[i] for i in range(0,len(weights))])
	mu = sum([vales[i]*weights[i] for i in range(0,len(weights))])
	mu = mu/sum_weights
	std_mu = []
	for j in range(0,11) :
		means = []
		size = 0.1 + (0.065)*j
		sum_weights = 0.0
		for i in range(0,100):
			resampled = np.random.randint(low=0.0, high=total, size=size)
			sum_weights = sum([weights[i] for i in resampled])
			mu = sum([charge_list[i]*weights[i] for i in resampled])
			mu = mu/sum_weights
			means.append(mu)
			std_mu.append(np.std(means,ddof=1))

	return std_mu_limit

def calc_charge_info(values,weights):

    
	"""
	Calculate the mean distance of the bin, (average) charge, and error in average charge

	Parameters
	----------
	total_charge_dict: dict
	Contains the total_charge_dict data.

	Returns
	-------
	charge_info: dict
	Contains the mean distance, charge, and error of each thing.
	"""
		
	# We use the standard IC procedure to calculate the statistical error for the
	# weighted average
	# There are three main terms
	# 1) The  weighted sum of charges, and its variance
	wsc = sum([ weights[i]*values[i] for i in range(0,len(values))])
	sw = sum(weights) 
	# Defining the weighted average
	mu = wsc/sw
	var_wsc = sum([(weights[i]*values[i])**2 for i in range(0,len(values))]) 
	# 2) The sum of weights and its variance 
	var_sw = sum([weights[i]**2 for i in range(0,len(values))])
	# 3) The covariance associated to both sums of weights 
	cov = sum([values[i]*(weights[i]**2) for i in range(0,len(values))])
	# The error in the weighted average of charges is given by the variance of the quantity 
	# wsc/sw; a division of two sums of weights (for two sets of weights that are correlated). \
	std_mu = var_wsc/(wsc)**2
	std_mu += var_sw/(sw)**2
	std_mu -= 2.*cov/(sw*wsc)
	std_mu = std_mu**0.5
	std_mu *= mu

	return mu , std_mu

def ComputeMeanandError(value,weight) :

	nelements = len(value)
	for i in range(0,nelements) :
		mean += value[i]/nelements

	for i in range(0,nelements) :
		sigma += ((value[i]-mean)**2.0)/nelements

	return mean, sigma**0.5


def ComputeWeightedMeanandError(value,weight):

	'''
	This function computes the standard error of the mean for a weighted set following the bootstrap validated formulation here:
	https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#:~:text=Statistical%20properties,-The%20weighted%20sample&text=can%20be%20called%20the%20standard,weights%20except%20one%20are%20zero.

	'''

	nelements = len(value)
	if len(weight) != nelements :
		print("error lists are not of equal length")
		return 0.0,0.0

	mean = 0.0
	sumweights = 0.0
	sumweights_sqr = 0.0
	n_nonzero = 0.0
	sigma = 0.

	for i in range(0,nelements) :
		sumweights = sumweights+weight[i]
		sumweights_sqr = sumweights_sqr+(weight[i]/sumweights)**2.0
		if weight[i] > 0.0 : n_nonzero += 1.0

	#Compute mean
	mean2 = 0.0
	for i in range(0,nelements) :
		mean += weight[i]*value[i]/sumweights
		mean2 += value[i]/nelements

	if n_nonzero == 0.0 : 
		print("Sum of weights is zero")
		return 0.0,0.0	

	for i in range(0,nelements) :
		sigma += (weight[i]*(value[i]-mean)/sumweights)**2.0

	#Compute standard error squared
	sigma /= (n_nonzero-1.0)/n_nonzero

	return mean, sigma**0.5


def OutputRoot(filename) :

	global weights_ic
	global weights_dc
	global reconstructedE_ic
	global reconstructedE_dc
	global zenith_ic
	global zenith_dc
	global EnergyTruth_ic
	global EnergyTruth_dc
	global ZenithTruth_ic
	global ZenithTruth_dc
	global totalcharge_ic
	global totalcharge_dc
	global recoEndpoint_ic
	global recoEndpoint_dc
	global borderDistance_ic
	global borderDistance_dc
	global stopLikeRatio_ic
	global stopLikeRatio_dc
	global recoLogL_ic
	global recoLogL_dc
	global directHits_ic
	global directHits_dc
	global HitsOut_ic
	global HitsOut_dc
	global binneddistance_dc
	global binneddistanceerror_dc
	global binnedcharge_dc
	global binnedchargeerror_dc
	global binneddistance_ic
	global binneddistanceerror_ic
	global binnedcharge_ic
	global binnedchargeerror_ic
	global ImpactAll_ic
	global ImpactAll_dc
	global Impact_seeMPE_ic
	global Impact_seeMPE_dc
	global Impact_vs_Zenith_ic
	global Impact_vs_Zenith_dc
	global TotalCharge_vs_Zenith_ic
	global TotalCharge_vs_Zenith_dc

	x_data_ic = array('f',binneddistance_ic)
	x_error_ic = array('f',binneddistanceerror_ic)
	y_data_ic = array('f',binnedcharge_ic)
	y_error_ic = array('f',binnedchargeerror_ic)
	y300_data_ic = array('f',binnedcharge300_ic)
	y300_error_ic = array('f',binnedcharge300error_ic)
	x_data_dc = array('f',binneddistance_dc)
	x_error_dc = array('f',binneddistanceerror_dc)
	y_data_dc = array('f',binnedcharge_dc)
	y_error_dc = array('f',binnedchargeerror_dc)
	y300_data_dc = array('f',binnedcharge300_dc)
	y300_error_dc = array('f',binnedcharge300error_dc)

	fout = ROOT.TFile.Open(filename+".root","RECREATE")

	fout.cd()

	Charge_Distance_IC = ROOT.TGraphErrors(len(x_data_ic),x_data_ic,y_data_ic,x_error_ic,y_error_ic)
	Charge_Distance_DC = ROOT.TGraphErrors(len(x_data_dc),x_data_dc,y_data_dc,x_error_dc,y_error_dc)
	Charge300_Distance_IC = ROOT.TGraphErrors(len(x_data_ic),x_data_ic,y300_data_ic,x_error_ic,y300_error_ic)
	Charge300_Distance_DC = ROOT.TGraphErrors(len(x_data_dc),x_data_dc,y300_data_dc,x_error_dc,y300_error_dc)
	TotalCharge_IC = ROOT.TH1F("TotalCharge_IC","",1000,0,3000)
	Zenith_IC = ROOT.TH1F("Zenith_IC","",200,-1.0,1.0)
	RecEnergy_IC = ROOT.TH1F("RecEnergy_IC","",1000,min(reconstructedE_ic)*0.9,max(reconstructedE_ic)*1.1)
	EndPointZ_IC = ROOT.TH1F("EndPointX_IC","",1000,min(recoEndpoint_ic)*0.9,max(recoEndpoint_ic)*1.1)
	BoarderDist_IC = ROOT.TH1F("BoarderDist_IC","",1000,min(borderDistance_ic)*0.9,max(borderDistance_ic)*1.1)
	StopLikeRatio_IC = ROOT.TH1F("StopLikeRatio_IC","",1000,min(stopLikeRatio_ic )*0.9,max(stopLikeRatio_ic)*1.1)
	RecoLogL_IC = ROOT.TH1F("RecoLogL_IC","",1000,min(recoLogL_ic)*0.9,max(recoLogL_ic)*1.1)
	DirectHits_IC = ROOT.TH1F("DirectHits_IC","",71,-0.5,70.5)
	HitsOut_IC = ROOT.TH1F("HitsOut_IC","",1000,0,max(HitsOut_ic)*1.1)
	TotalCharge_DC = ROOT.TH1F("TotalCharge_DC","",1000,0,3000)
	Zenith_DC = ROOT.TH1F("Zenith_DC","",200,-1.0,1.0)
	RecEnergy_DC = ROOT.TH1F("RecEnergy_DC","",1000,min(reconstructedE_dc)*0.9,max(reconstructedE_dc)*1.1)
	EndPointZ_DC = ROOT.TH1F("EndPointX_DC","",1000,min(recoEndpoint_dc)*0.9,max(recoEndpoint_dc)*1.1)
	BoarderDist_DC = ROOT.TH1F("BoarderDist_DC","",1000,min(borderDistance_dc)*0.9,max(borderDistance_dc)*1.1)
	StopLikeRatio_DC = ROOT.TH1F("StopLikeRatio_DC","",1000,min(stopLikeRatio_dc)*0.9,max(stopLikeRatio_dc)*1.1)
	RecoLogL_DC = ROOT.TH1F("RecoLogL_DC","",1000,min(recoLogL_dc)*0.9,max(recoLogL_dc)*1.1)
	DirectHits_DC = ROOT.TH1F("DirectHits_DC","",(1+max(directHits_dc)-min(directHits_dc)),min(directHits_dc)-0.5,max(directHits_dc)+0.5)
	HitsOut_DC = ROOT.TH1F("HitsOut_DC","",1000,min(HitsOut_dc)*0.9,max(HitsOut_dc)*1.1)

	for i in range(0,len(weights_ic)) :
		TotalCharge_IC.Fill(totalcharge_ic[i],weights_ic[i])
		Zenith_IC.Fill(zenith_ic[i],weights_ic[i])
		RecEnergy_IC.Fill(reconstructedE_ic[i],weights_ic[i])
		EndPointZ_IC.Fill(recoEndpoint_ic[i],weights_ic[i])
		BoarderDist_IC.Fill(borderDistance_ic[i],weights_ic[i])
		StopLikeRatio_IC.Fill(stopLikeRatio_ic[i],weights_ic[i])
		RecoLogL_IC.Fill(recoLogL_ic[i],weights_ic[i])
		DirectHits_IC.Fill(directHits_ic[i],weights_ic[i])
		HitsOut_IC.Fill(HitsOut_ic[i],weights_ic[i])
	for i in range(0,len(weights_dc)) :
		TotalCharge_DC.Fill(totalcharge_dc[i],weights_dc[i])
		Zenith_DC.Fill(zenith_dc[i],weights_dc[i])
		RecEnergy_DC.Fill(reconstructedE_dc[i],weights_dc[i])
		EndPointZ_DC.Fill(recoEndpoint_dc[i],weights_dc[i])
		BoarderDist_DC.Fill(borderDistance_dc[i],weights_dc[i])
		StopLikeRatio_DC.Fill(stopLikeRatio_dc[i],weights_dc[i])
		RecoLogL_DC.Fill(recoLogL_dc[i],weights_dc[i])
		DirectHits_DC.Fill(directHits_dc[i],weights_dc[i])
		HitsOut_DC.Fill(HitsOut_dc[i],weights_dc[i])

	Charge_Distance_IC.Write("Charge_Distance_IC")
	Charge_Distance_DC.Write("Charge_Distance_DC")
	Charge300_Distance_IC.Write("Charge300_Distance_IC")
	Charge300_Distance_DC.Write("Charge300_Distance_DC")
	TotalCharge_IC.Write()
	Zenith_IC.Write()
	RecEnergy_IC.Write()
	EndPointZ_IC.Write()
	BoarderDist_IC.Write()
	StopLikeRatio_IC.Write()
	RecoLogL_IC.Write()
	DirectHits_IC.Write()
	HitsOut_IC.Write()
	TotalCharge_DC.Write()
	Zenith_DC.Write()
	RecEnergy_DC.Write()
	EndPointZ_DC.Write()
	BoarderDist_DC.Write()
	StopLikeRatio_DC.Write()
	RecoLogL_DC.Write()
	DirectHits_DC.Write()
	HitsOut_DC.Write()
	ImpactAll_ic.Write()
	ImpactAll_dc.Write()
	Impact_seeMPE_ic.Write()
	Impact_seeMPE_dc.Write()
	Impact_vs_Zenith_ic.Write()
	Impact_vs_Zenith_dc.Write()
	TotalCharge_vs_Zenith_ic.Write()
	TotalCharge_vs_Zenith_dc.Write()

	fout.Close()

def OutputHDF5(filename,args) :

	global binneddistance_dc
	global binneddistanceerror_dc
	global binnedcharge_dc
	global binnedchargeerror_dc
	global binneddistance_ic
	global binneddistanceerror_ic
	global binnedcharge_ic
	global binnedchargeerror_ic

	h5file = open_file(filename+".h5", mode="w", title="DOM Calibration HDF5 File")
	
	icecube = h5file.create_table('/', 'icecube', DataPoint, "IceCube Charge vs Distance data")
	deapcore = h5file.create_table('/', 'deepcore', DataPoint, "DeepCore Charge vs Distance data")
	state = h5file.create_table('/','state',State,"Information on the data cuts")

	nelements = len(binneddistance_ic)

	icrow = icecube.row
	dcrow = deapcore.row
	staterow = state.row

	#h5file.root._v_attrs.event_cuts = event_cuts
    #h5file.root._v_attrs.dom_cuts = dom_cuts

	staterow['data'] = args.data
	staterow['eff'] = args.eff
	staterow['flux'] = args.flux
	staterow['zenithmin'] = args.zenithrange[0]
	staterow['zenithmax'] = args.zenithrange[1]
	staterow['energymin'] = args.energyrange[0]
	staterow['impactanglemin'] = args.energyrange[1]
	staterow['impactanglemax'] = args.impactrange[0]
	staterow['trackendpoint'] = args.impactrange[1]
	staterow['cherdistmin'] = args.cherdist[0]
	staterow['cherdistmax'] = args.cherdist[1]
	staterow['binwidth'] = args.binwidth
	staterow.append()

	for i in range(0,nelements) :
		icrow['meandistance'] = binneddistance_ic[i]
		icrow['sigmadistance'] = binneddistanceerror_ic[i]
		icrow['meancharge'] = binnedcharge_ic[i]
		icrow['sigmacharge'] = binnedchargeerror_ic[i]
		icrow.append()
		dcrow['meandistance'] = binneddistance_dc[i]
		dcrow['sigmadistance'] = binneddistanceerror_dc[i]
		dcrow['meancharge'] = binnedcharge_dc[i]
		dcrow['sigmacharge'] = binnedchargeerror_dc[i]
		dcrow.append()

	h5file.close()



if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--data', help='Directory of data files.',type=str,
				default = '/data/user/sanchezh/IC86_2015/Final_Level2_IC86_MPEFit_*.h5')
	parser.add_argument('-e', '--eff', help='efficiency to be used or data for data.', type = str,
				default = "eff100")
	parser.add_argument('-o', '--output', help='Name of output file.', type=str,
				default = "out.root")
	parser.add_argument('-f', '--flux', help='Name of flux model.', type=str,
				default = "data")
	parser.add_argument('-z', '--zenithrange', help='Range of muon Zeniths', type = float,
				nargs = 2,  default = [40.,90.])
	parser.add_argument('-p', '--energyrange', help='Range of muon Energies', type = float,
				nargs = 2, default = [0.0, 9999999.00])
	parser.add_argument('-i','--impactrange',help='Range of DOM impact parameters to include', 
				type = float, nargs = "+", default = [0.0,180.0])
	parser.add_argument('-t','--trackendpoint',help='Distance from track end point to include',
				type = float, default = 100.)
	parser.add_argument('-c','--cherdist', help='Distance from track to include', type = float, 
				nargs = 2, default = [0.0,200.])
	parser.add_argument('-w','--binwidth', help='Width to bin distances', type = float, default = 20.0)
	parser.add_argument('-r','--residual', help='time residual region', type = float, 
				nargs = 2, default = [-15.0,75.0])
	parser.add_argument('-x','--hitsout',help='Max number of hits outside analysis region',type = int, default = 20)
	parser.add_argument('-b','--boarder',help='Distance from bottom and side of detector.',type = float, 
						nargs = 2, default = [-400.0, 100.0])
	parser.add_argument('-n', '--nhits', help='Min number of direct hit DOMs and Max number of Outside analysis hits', 
						type = int, nargs = 2, default = [5,20])
	parser.add_argument('-l', '--likelihood', help='Fit likelyhoods, FiniteReco Likelihood ratio and SplineMPE Rlogl', type = float,
						nargs = 2, default = [10.,10.])


	args = parser.parse_args()
	weightname = 'weight_'+args.flux

	# compute how many 20m bins to use.
	nbins = int(args.cherdist[1] / args.binwidth)

	bin_DomCharge_ic  = [[] for i in range(nbins)]
	bin_DomCharge300_ic  = [[] for i in range(nbins)]
	bin_weights_ic = [[] for i in range(nbins)]
	bin_distance_ic = [[] for i in range(nbins)]
	bin_DomCharge_dc  = [[] for i in range(nbins)]
	bin_DomCharge300_dc  = [[] for i in range(nbins)]
	bin_weights_dc = [[] for i in range(nbins)]
	bin_distance_dc = [[] for i in range(nbins)]

	DC_Strings = [81,82,83,84,85,86]
	IC_Strings = [17,18,19,25,26,27,28,34,38,44,47,56,54,55]

	files_dir = args.data
	file_list_aux = os.listdir(files_dir)
	file_list_h5 = [x for x in file_list_aux if '.h5' in x]
	file_list = []
	if args.flux == "data" :
		file_list = [x for x in file_list_h5 if (args.eff in x and os.path.getsize(files_dir+x) > 50000000 )]
	else :
		file_list = [x for x in file_list_h5 if (args.eff in x)]

	nfiles = len(file_list)

	flux = GaisserH4a()
	if args.flux == "GaisserH3a" : flux = GaisserH3a()
	elif args.flux == "GaisserH4a_IT" : flux = GaisserH4a_IT()
	elif args.flux == "GaisserHillas" : flux = GaisserHillas()
	elif args.flux == "Hoerandel" : flux = Hoerandel()
	elif args.flux == "Hoerandel5" : flux = Hoerandel5()
	elif args.flux == "Hoerandel_IT" : flux = Hoerandel_IT()

	eventcount = 0
	max_weight = 0.0
	totalevent = 0
	#generator = weighting.from_simprod(21269,False,'vm-simprod2.icecube.wisc.edu')
	#generator = weighting.icetop_mc_weights(21269,'/home/tmcelroy/icecube/domeff/datasetConfig.json')
	nfiles = len(file_list)
	
	for filename in file_list :
		h5file = open_file(files_dir+filename, mode="r")
		domtable = h5file.root.doms
		eventtable = h5file.root.events
		runtable = h5file.root.runinfo

		domindex = 0

		for event in eventtable.iterrows() :

			totalevent += 1
			#Energy Cut
			if event['reco/energy'] < args.energyrange[0] or event['reco/energy'] > args.energyrange[1] : 
				#print("Event killed by energy Cut")
				#print(event['reco/energy'])
				continue

			#Zenith Cut
			if event['reco/dir/zenith'] < args.zenithrange[0]*ROOT.TMath.Pi()/180. or event['reco/dir/zenith'] > args.zenithrange[1]*ROOT.TMath.Pi()/180. : 
				#print("Event Killed by Zenith Cut")
				#print(event['reco/dir/zenith'])
				continue

			#Stopping Point Cut
			if event['recoEndPoint/z'] < args.boarder[0] :
				#print("Event killed by Bottom Distance Cut")
				#print(event['recoEndPoint/z'])
				#print("mctruth = %f" %(event['truthEndPoint/z']))
				#print("event likelihood = %f" % (event['stopLikeRatio']))
				continue

			if event['borderDistance'] < args.boarder[1] :
				#print("Event killed by Detector Edge Cut")
				#print(event['borderDistance'])
				#print(event['truthBorderDistance'])
				#print("event likelihood = %f" % (event['stopLikeRatio']))
				continue 

			#Likelihood cuts
			if event['stopLikeRatio'] < args.likelihood[0] :
				#print("Event cut by likelihood ratio cut")
				#print("event likelihood = %f" % (event['stopLikeRatio']))
				continue

			 #Likelihood cuts
                        if event['recoLogL'] > args.likelihood[1] :
				#print("Event killed by Likelihood check")
				#print(event['recoLogL'])
                                continue

			#direct hists
			if event['directHits'] < args.nhits[0]:
				#print("Event killed by N Direct Hists Cut")
				#print(event['directHits'])
				continue

			#print("Event Passed")
			eventcount += 1
			weight = 1.0
			if args.flux != "data" :
				pflux = flux(event['corsika/primaryEnergy'],event['corsika/primaryType'])
				energy_integral = event['corsika/energyPrimaryMax']**(event['corsika/primarySpectralIndex']+1)
				energy_integral = energy_integral - event['corsika/energyPrimaryMin']**(event['corsika/primarySpectralIndex']+1)
				energy_integral = energy_integral / (event['corsika/primarySpectralIndex']+1)
				energy_weight = event['corsika/primaryEnergy']**event['corsika/primarySpectralIndex']
				energy_weight = pflux*energy_integral/energy_weight*event['corsika/areaSum']
				weight = energy_weight/(event['corsika/nEvents'])
				if weight > max_weight:
					max_weight = weight
				#gen = generator(event['corsika/primaryEnergy'],event['corsika/primaryType'])
				#weight = pflux/gen
				#print("weight = %f" % (weight))

			first_ic = True
			first_dc = True
			for dom in domtable.iterrows(domindex) :
				if dom['eventId'] < event['eventId'] :
					domindex += 1
					continue
				elif dom['eventId'] == event['eventId'] :
					domindex += 1
				else :
					#print("new event")
					break
				if dom['impactAngle'] < args.impactrange[0]*3.14/180. or  dom['impactAngle'] > args.impactrange[1]*3.14/180.:
					#print("DOM killed by Impact Angle Cut")
					#print(dom['impactAngle']) 
					continue
				if dom['distAboveEndpoint'] < args.trackendpoint :
					#print("DOM killed by Dist Above End Cut")
					#print(dom['distAboveEndpoint']) 
					continue
				if dom['recoDist'] < args.cherdist[0] or dom['recoDist'] > args.cherdist[1] :
					#print("DOM killed by Distance from Track Cut")
					#print(dom['recoDist']) 
					continue 
				i_dist = (int)(dom['recoDist']/args.binwidth)
				if i_dist < 0 or i_dist > nbins-1 :
					#print("DOM out of bin range")
					continue
				if dom['string'] in DC_Strings :
					#print("DC DOM Passed")
					#if event['dcHitsOut']> args.nhits[1] or event['dcHitsOut'] < 1 :
					#	continue;
					#I want same overall events used for both.
					if event['icHitsOut']> args.nhits[1] or event['icHitsOut'] < 1 :
						continue;
					#print("DC Distance Charge")
                                        #print(dom['recoDist'])
					#print(dom['totalCharge'])
					if first_dc :
						weights_dc.append(weight)
						reconstructedE_dc.append(event['reco/energy'])
						zenith_dc.append(ROOT.TMath.Cos(event['reco/dir/zenith']))
						totalcharge_dc.append(event['totalCharge'])
						recoEndpoint_dc.append(event['recoEndPoint/z'])
						borderDistance_dc.append(event['borderDistance'])
						stopLikeRatio_dc.append(event['stopLikeRatio'])
						recoLogL_dc.append(event['recoLogL'])
						directHits_dc.append(event['directHits'])
						HitsOut_dc.append(event['dcHitsOut'])
						TotalCharge_vs_Zenith_dc.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['totalCharge'])
						first_dc = False
					bin_DomCharge_dc[i_dist].append(dom['totalCharge'])
					bin_DomCharge300_dc[i_dist].append(dom['totalCharge_300ns'])
					bin_weights_dc[i_dist].append(weight)
					bin_distance_dc[i_dist].append(dom['recoDist'])
					ImpactAll_dc.Fill(ROOT.TMath.Cos(dom['impactAngle']),weight)
					Impact_vs_Zenith_dc.Fill(ROOT.TMath.Cos(dom['impactAngle']),ROOT.TMath.Cos(event['reco/dir/zenith']),weight)
					Impact_seeMPE_dc.Fill(ROOT.TMath.Cos(dom['impactAngle']),weight)
				if dom['string'] in IC_Strings :
					#print("IC DOM Passed")
					if event['icHitsOut'] > args.nhits[1] or event['icHitsOut']< 1:
						continue;
					#print("IC Distance Charge")
					#print(dom['recoDist'])
					#print(dom['totalCharge'])
					if first_ic :
						weights_ic.append(weight)
						reconstructedE_ic.append(event['reco/energy'])
						zenith_ic.append(ROOT.TMath.Cos(event['reco/dir/zenith']))
						totalcharge_ic.append(event['totalCharge'])
						recoEndpoint_ic.append(event['recoEndPoint/z'])
						borderDistance_ic.append(event['borderDistance'])
						stopLikeRatio_ic.append(event['stopLikeRatio'])
						recoLogL_ic.append(event['recoLogL'])
						directHits_ic.append(event['directHits'])
						HitsOut_ic.append(event['icHitsOut'])
						TotalCharge_vs_Zenith_ic.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['totalCharge'])
						first_ic = False
					bin_DomCharge_ic[i_dist].append(dom['totalCharge'])
					bin_DomCharge300_ic[i_dist].append(dom['totalCharge_300ns'])
					bin_weights_ic[i_dist].append(weight)
					bin_distance_ic[i_dist].append(dom['recoDist'])
					Impact_vs_Zenith_ic.Fill(ROOT.TMath.Cos(dom['impactAngle']),ROOT.TMath.Cos(event['reco/dir/zenith']),weight)
					ImpactAll_ic.Fill(ROOT.TMath.Cos(dom['impactAngle']),weight)
					if dom['totalCharge'] > 0.0 :
						Impact_seeMPE_ic.Fill(ROOT.TMath.Cos(dom['impactAngle']),weight)

		h5file.close()
	
	print("Total number of events = %d/%d" % (eventcount,totalevent))
	
	binneddistance_dc = np.zeros(nbins,dtype=float)
	binneddistanceerror_dc = np.zeros(nbins,dtype=float)
	binnedcharge_dc = np.zeros(nbins,dtype=float)
	binnedchargeerror_dc = np.zeros(nbins,dtype=float)
	binnedcharge300_dc = np.zeros(nbins,dtype=float)
	binnedcharge300error_dc = np.zeros(nbins,dtype=float)
	binneddistance_ic = np.zeros(nbins,dtype=float)
	binneddistanceerror_ic = np.zeros(nbins,dtype=float)
	binnedcharge_ic = np.zeros(nbins,dtype=float)
	binnedchargeerror_ic = np.zeros(nbins,dtype=float)
	binnedcharge300_ic = np.zeros(nbins,dtype=float)
	binnedcharge300error_ic = np.zeros(nbins,dtype=float)

	for i in range(0,len(bin_distance_dc)):
		if len(bin_weights_dc[i]) > 0 :
			binneddistance_dc[i] , binneddistanceerror_dc[i] = calc_charge_info(bin_distance_dc[i],bin_weights_dc[i])
			binnedcharge_dc[i], binnedchargeerror_dc[i] = calc_charge_info(bin_DomCharge_dc[i],bin_weights_dc[i])
			binnedcharge300_dc[i], binnedcharge300error_dc[i] = calc_charge_info(bin_DomCharge300_dc[i],bin_weights_dc[i])
		if len(bin_weights_ic[i]) > 0 :
			binneddistance_ic[i] , binneddistanceerror_ic[i] = calc_charge_info(bin_distance_ic[i],bin_weights_ic[i])
			binnedcharge_ic[i], binnedchargeerror_ic[i] = calc_charge_info(bin_DomCharge_ic[i],bin_weights_ic[i])
			binnedcharge300_ic[i], binnedcharge300error_ic[i] = calc_charge_info(bin_DomCharge300_ic[i],bin_weights_ic[i])
		

	#outfilenamelist = args.output.split(".",1)
	#if "root" in outfilenamelist[1] :
	OutputRoot(args.output)

	#elif "h5" in outfilenamelist[1] :
	OutputHDF5(args.output,args)
	#elif "pdf" in outfilenamelist[1] :
	#	print("not yet supported")


