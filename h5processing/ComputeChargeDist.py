import os, sys
from tables import open_file
import ROOT
import numpy as np
import argparse
from tables import open_file
from event import *
from array import array

binneddistance_dc = np.zeros(1,dtype=float)
binneddistanceerror_dc = np.zeros(1,dtype=float)
binnedcharge_dc = np.zeros(1,dtype=float)
binnedchargeerror_dc = np.zeros(1,dtype=float)
binneddistance_ic = np.zeros(1,dtype=float)
binneddistanceerror_ic = np.zeros(1,dtype=float)
binnedcharge_ic = np.zeros(1,dtype=float)
binnedchargeerror_ic = np.zeros(1,dtype=float)



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
	if len(values) == 0 :
		print("length 0")
		return 0.0, 1.0
	wsc = sum([ weights[i]*values[i] for i in range(0,len(values))])
	sw = sum(weights)
	if sw == 0.0 or wsc == 0.0 :
		return 0.0, 1.0 
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
	global binneddistance_dc
	global binneddistanceerror_dc
	global binnedcharge_dc
	global binnedchargeerror_dc
	global binneddistance_ic
	global binneddistanceerror_ic
	global binnedcharge_ic
	global binnedchargeerror_ic
	x_data_ic = array('f',binneddistance_ic)
	x_error_ic = array('f',binneddistanceerror_ic)
	y_data_ic = array('f',binnedcharge_ic)
	y_error_ic = array('f',binnedchargeerror_ic)
	x_data_dc = array('f',binneddistance_dc)
	x_error_dc = array('f',binneddistanceerror_dc)
	y_data_dc = array('f',binnedcharge_dc)
	y_error_dc = array('f',binnedchargeerror_dc)

	fout = ROOT.TFile.Open(filename+".root","RECREATE")

	fout.cd()

	Charge_Distance_IC = ROOT.TGraphErrors(len(x_data_ic),x_data_ic,y_data_ic,x_error_ic,y_error_ic)
	Charge_Distance_DC = ROOT.TGraphErrors(len(x_data_dc),x_data_dc,y_data_dc,x_error_dc,y_error_dc)

	Charge_Distance_IC.Write("Charge_Distance_IC")
	Charge_Distance_DC.Write("Charge_Distance_DC")

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

	print(icecube)

	nelements = len(binneddistance_ic)

	icrow = icecube.row
	dcrow = deapcore.row

	for i in range(0,nelements) :
		print(binneddistance_ic[i])
		print(binnedcharge_ic[i])
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
	parser.add_argument('-c','--cherdist', help='Distance from track to include', type = float, 
				nargs = 2, default = [0.0,200.])
	parser.add_argument('-w','--binwidth', help='Width to bin distances', type = float, default = 20.0)


	args = parser.parse_args()

	# compute how many 20m bins to use.
	nbins = int(args.cherdist[1] / args.binwidth)

	bin_DomCharge_ic  = [[] for i in range(nbins)]
	bin_weights_ic = [[] for i in range(nbins)]
	bin_distance_ic = [[] for i in range(nbins)]
	bin_DomCharge_dc  = [[] for i in range(nbins)]
	bin_weights_dc = [[] for i in range(nbins)]
	bin_distance_dc = [[] for i in range(nbins)]


	files_dir = args.data
	file_list_aux = os.listdir(files_dir)
	file_list = list()
#	for (dirpath, dirnames, filenames) in os.walk(files_dir):
#		for eff in args.eff :
#			file_list += [os.path.join(dirpath,x) for x in filenames if '.h5' in x and eff in x]
	for filename in file_list_aux :
		if '.h5' in filename and args.eff in filename :
			file_list.append(filename)
    #remove duclicates
	file_list = list(set(file_list))

	nfiles = len(file_list)

	print(nfiles)

	jobfilelist = os.listdir("/data/user/tmcelroy/domeff/jobscripts/")
	
	for filename in file_list :
		h5file = 0

		if args.output in filename :
			continue
		try :
			h5file = open_file(files_dir+filename, mode="r")
		except : continue

		print(filename)
		weight = 1.0
                for jobfilename in jobfilelist :
                        if weight < 1.0 :
                                break
                        if filename.split(".",5)[0] in jobfilename :
                                jobfile = open("/data/user/tmcelroy/domeff/jobscripts/"+jobfilename,"r")
                                lines = jobfile.readlines()
                                for line in lines :
                                        if "-y" in line :
                                                if(float(line.split(" ",100)[-1]) > 0.1) :
                                                        weight = 0.1
                                jobfile.close()


		#print(filename)
		ictable = h5file.root.icecube
		dctable = h5file.root.deepcore

		columbnames = names = ictable.coldescrs.keys()
		#if "truthDist" not in columbnames :
		#	continue;
		#print(columbnames)	
		for dom in ictable.iterrows() :
	 
			i_dist = (int)(dom['dist']/args.binwidth)
			if i_dist < 0 or i_dist > nbins-1 :
				continue
			bin_distance_ic[i_dist].append(dom['dist'])
			#print(dom['charge'])
			bin_DomCharge_ic[i_dist].append(dom['charge'])
			if "weight" in columbnames :
				bin_weights_ic[i_dist].append(dom['weight']*weight)
			else :
				bin_weights_ic[i_dist].append(1.0)

		for dom in dctable.iterrows() :
	 
			i_dist = (int)(dom['dist']/args.binwidth)
			if i_dist < 0 or i_dist > nbins-1 :
				continue
			bin_distance_dc[i_dist].append(dom['dist'])
			#print(dom['charge'])
			bin_DomCharge_dc[i_dist].append(dom['charge'])
			if "weight" in columbnames :
				bin_weights_dc[i_dist].append(dom['weight']*weight)
			else :
				bin_weights_dc[i_dist].append(1.0)

		h5file.close()
	
	binneddistance_dc = np.zeros(nbins,dtype=float)
	binneddistanceerror_dc = np.zeros(nbins,dtype=float)
	binnedcharge_dc = np.zeros(nbins,dtype=float)
	binnedchargeerror_dc = np.zeros(nbins,dtype=float)
	binneddistance_ic = np.zeros(nbins,dtype=float)
	binneddistanceerror_ic = np.zeros(nbins,dtype=float)
	binnedcharge_ic = np.zeros(nbins,dtype=float)
	binnedchargeerror_ic = np.zeros(nbins,dtype=float)

	for i in range(0,len(bin_distance_dc)):
		if len(bin_weights_dc[i]) > 0 :
			binneddistance_dc[i] , binneddistanceerror_dc[i] = calc_charge_info(bin_distance_dc[i],bin_weights_dc[i])
			binnedcharge_dc[i], binnedchargeerror_dc[i] = calc_charge_info(bin_DomCharge_dc[i],bin_weights_dc[i])
		else :
			print("wight length 0")	
		if len(bin_weights_ic[i]) > 0 :
			binneddistance_ic[i] , binneddistanceerror_ic[i] = calc_charge_info(bin_distance_ic[i],bin_weights_ic[i])
			binnedcharge_ic[i], binnedchargeerror_ic[i] = calc_charge_info(bin_DomCharge_ic[i],bin_weights_ic[i])

	OutputHDF5(args.output,args)
	OutputRoot(args.output)

