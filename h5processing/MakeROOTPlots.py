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

TruthVsRecoDistance_IC = ROOT.TH2F("TruthVsRecoDistance_IC","",200,0.0,200.0,200,0.0,200.0)
TruthVsRecoDistance_DC = ROOT.TH2F("TruthVsRecoDistance_DC","",200,0.0,200.0,200,0.0,200.0)

TruePE_Vs_MeasuredPE_IC = ROOT.TH2F("TruePE_Vs_MeasuredPE_IC","",200,0.0,200.0,100,0.0,2.0)
TruePE_Vs_MeasuredPE_DC = ROOT.TH2F("TruePE_Vs_MeasuredPE_DC","",200,0.0,200.0,100,0.0,2.0)

TotalCharge_IC = ROOT.TH1F("TotalCharge_IC","",1000,0,3000)
Zenith_IC = ROOT.TH1F("Zenith_IC","",200,-1.0,1.0)
RecEnergy_IC = ROOT.TH1F("RecEnergy_IC","",1000,0.*0.9,1000000*1.1)
EndPointZ_IC = ROOT.TH1F("EndPointX_IC","",1000,-500.*1.1,500.*1.1)
BoarderDist_IC = ROOT.TH1F("BoarderDist_IC","",1000,0.0,500.*1.1)
StopLikeRatio_IC = ROOT.TH1F("StopLikeRatio_IC","",1000,-100.,100)
RecoLogL_IC = ROOT.TH1F("RecoLogL_IC","",1000,-100.,100.)
DirectHits_IC = ROOT.TH1F("DirectHits_IC","",71,-0.5,70.5)
HitsOut_IC = ROOT.TH1F("HitsOut_IC","",1000,0,100*1.1)

TotalCharge_DC = ROOT.TH1F("TotalCharge_DC","",1000,0,3000)
Zenith_DC = ROOT.TH1F("Zenith_DC","",200,-1.0,1.0)
RecEnergy_DC = ROOT.TH1F("RecEnergy_DC","",1000,0.*0.9,1000000*1.1)
EndPointZ_DC = ROOT.TH1F("EndPointX_DC","",1000,-500.*1.1,500.*1.1)
BoarderDist_DC = ROOT.TH1F("BoarderDist_DC","",1000,0.0,500.*1.1)
StopLikeRatio_DC = ROOT.TH1F("StopLikeRatio_DC","",1000,-100.,100)
RecoLogL_DC = ROOT.TH1F("RecoLogL_DC","",1000,-100.,100.)
DirectHits_DC = ROOT.TH1F("DirectHits_DC","",71,-0.5,70.5)
HitsOut_DC = ROOT.TH1F("HitsOut_DC","",1000,0,100*1.1)

ImpactAll_ic = ROOT.TH1F("ImpactAll_IC","",1000,-1.0,1.0)
ImpactAll_dc = ROOT.TH1F("ImpactAll_DC","",1000,-1.0,1.0)
Impact_seeMPE_ic = ROOT.TH1F("ImpactseeMPE_IC","",1000,-1.0,1.0)
Impact_seeMPE_dc = ROOT.TH1F("ImpactseeMPE_DC","",1000,-1.0,1.0)

Impact_vs_Zenith_ic = ROOT.TH2F("Impact_vs_Zenith_IC","",200,-1.0,1.0,200,-1.0,1.0)
Impact_vs_Zenith_dc = ROOT.TH2F("Impact_vs_Zenith_DC","",200,-1.0,1.0,200,-1.0,1.0)

TotalCharge_vs_Zenith_ic = ROOT.TH2F("TotalCharge_vs_Zenith_IC","",200,-1.0,1.0,1000,0.0,3000.0)
TotalCharge_vs_Zenith_dc = ROOT.TH2F("TotalCharge_vs_Zenith_DC","",200,-1.0,1.0,1000,0.0,3000.0)

StoppingZ_vs_Zenith_ic = ROOT.TH2F("StoppingZ_vs_Zenith_ic","",200,-1.0,1.0,1000,-500.,500.)
StoppingZ_vs_Zenith_dc = ROOT.TH2F("StoppingZ_vs_Zenith_dc","",200,-1.0,1.0,1000,-500.,500.)
BorderDist_vs_Zenith_ic = ROOT.TH2F("BorderDist_vs_Zenith_ic","",200,-1.0,1.0,1000,-500.,500.)
BorderDist_vs_Zenith_dc = ROOT.TH2F("BorderDist_vs_Zenith_dc","",200,-1.0,1.0,1000,-500.,500.)
NChannel_vs_Zenith_ic = ROOT.TH2F("NChannel_vs_Zenith_ic","",200,-1.0,1.0,5000,0,5000.)
NChannel_vs_Zenith_dc = ROOT.TH2F("NChannel_vs_Zenith_dc","",200,-1.0,1.0,5000,0,5000.)
StopLike_vs_Zenith_ic = ROOT.TH2F("StopLike_vs_Zenith_ic","",200,-1.0,1.0,500,0,500.)
StopLike_vs_Zenith_dc = ROOT.TH2F("StopLike_vs_Zenith_dc","",200,-1.0,1.0,500,0,500.)
rLogL_vs_Zenith_ic = ROOT.TH2F("rLogL_vs_Zenith_ic","",200,-1.0,1.0,500,0,500.)
rLogL_vs_Zenith_dc = ROOT.TH2F("rLogL_vs_Zenith_dc","",200,-1.0,1.0,500,0,500.)


zenith_all = ROOT.TH1F("zenith_all","",100,0.0,1.0)
zenith_all_line = ROOT.TH1F("zenith_all_line","",100,0.0,1.0)
zenith_all_spe = ROOT.TH1F("zenith_all_spe","",100,0.0,1.0)
zenith_all_mpe = ROOT.TH1F("zenith_all_mpe","",100,0.0,1.0)
zenith_zenith = ROOT.TH1F("zenith_zenith","",100,0.0,1.0)
zenith_endpointz = ROOT.TH1F("zenith_endpointz","",100,0.0,1.0)
zenith_boarder = ROOT.TH1F("zenith_boarder","",100,0.0,1.0)
zenith_stopratio = ROOT.TH1F("zenith_stopratio","",100,0.0,1.0)
zenith_recoLogL = ROOT.TH1F("zenith_recoLogL","",100,0.0,1.0)
zenith_directHits = ROOT.TH1F("zenith_directHits","",100,0.0,1.0)

TimeResidual_IC = []
TimeResidual_DC = []

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--data', help='Directory of data files.',type=str,
				default = '/data/user/sanchezh/IC86_2015/Final_Level2_IC86_MPEFit_*.h5')
	parser.add_argument('-e', '--eff', help='efficiency to be used or data for data.', type = str, nargs = '+',
				default = ["eff100"])
	parser.add_argument('-o', '--output', help='Name of output file.', type=str,
				default = "out.root")
	parser.add_argument('-f', '--flux', help='Name of flux model.', type=str,
				default = "data")
	parser.add_argument('-z', '--zenithrange', help='Range of muon Zeniths', type = float,
				nargs = 2,  default = [40.,70.])
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
	parser.add_argument('-y', '--skim',help="skim fraction",type = float, default = 1.0)


	args = parser.parse_args()
	weightname = 'weight_'+args.flux

	# compute how many 20m bins to use.
	nbins = int(args.cherdist[1] / args.binwidth)

	DC_Strings = [81,82,83,84,85,86]
	IC_Strings = [17,18,19,25,26,27,28,34,38,44,47,56,54,55]

	files_dir = args.data
	file_list_aux = os.listdir(files_dir)
	file_list = list()
	for (dirpath, dirnames, filenames) in os.walk(files_dir):
		for eff in args.eff :
			file_list += [os.path.join(dirpath,x) for x in filenames if '.h5' in x and eff in x]
    #remove duclicates
	file_list = list(set(file_list))
#	if args.flux == "data" :
#		file_list = [x for x in file_list_h5 if (args.eff in x and os.path.getsize(files_dir+x) > 12000000 )]
#	else :
#		file_list = [x for x in file_list_h5 if (args.eff in x)]

	nfiles = len(file_list)

	print(nfiles)

	fout = ROOT.TFile.Open(args.output+".root","RECREATE")

	fout.cd()

	flux = GaisserH4a()
	if args.flux == "GaisserH3a" : flux = GaisserH3a()
	elif args.flux == "GaisserH4a_IT" : flux = GaisserH4a_IT()
	elif args.flux == "GaisserHillas" : flux = GaisserHillas()
	elif args.flux == "Hoerandel" : flux = Hoerandel()
	elif args.flux == "Hoerandel5" : flux = Hoerandel5()
	elif args.flux == "Hoerandel_IT" : flux = Hoerandel_IT()

	eventcount = 0
	#generator = weighting.from_simprod(21269,False,'vm-simprod2.icecube.wisc.edu')
	#generator = weighting.icetop_mc_weights(21269,'/home/tmcelroy/icecube/domeff/datasetConfig.json')
	nfiles = len(file_list)

	rand = ROOT.TRandom3()

	for filename in file_list :
		h5file = 0
		try :
			h5file = open_file(filename, mode="r")
		except : continue
		
		domtable = h5file.root.doms
		eventtable = h5file.root.events
		runtable = h5file.root.runinfo

		domindex = 0
		eventindex = -1

		for event in eventtable.iterrows() :


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
				weight = weight/nfiles

			zenith_all.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),weight)

			eventindex += 1
			if args.flux == "data" :
				if rand.Uniform()>args.skim : continue

			if args.flux == "data" :			
				if event['filterMask/SunFilter_13'] : continue
				if event['filterMask/MoonFilter_13'] : continue 

			#Energy Cut
			if event['totalCharge'] < args.energyrange[0] or event['totalCharge'] > args.energyrange[1] : 
				continue

			#Zenith Cut
			if event['reco/dir/zenith'] < args.zenithrange[0]*ROOT.TMath.Pi()/180. or event['reco/dir/zenith'] > args.zenithrange[1]*ROOT.TMath.Pi()/180. : 
				continue
		
			#Stopping Point Cut
			if event['recoEndPoint/z'] < args.boarder[0] :
				continue

			if event['borderDistance'] < args.boarder[1] :
				continue 

			#Likelihood cuts
			if event['stopLikeRatio'] < args.likelihood[0] :
				continue

			 #Likelihood cuts
            		if event['recoLogL'] > args.likelihood[1] :
                		continue

			#direct hists
			if event['directHits'] < args.nhits[0]:
				continue

			if event['icHitsOut']> args.nhits[1] or event['icHitsOut'] < 1 :
                                        continue;

			first_ic = True
			first_dc = True
			for dom in domtable.iterrows(domindex) :


				if dom['eventId'] < int(eventindex) : #event['eventId'] :
					domindex += 1
					continue
				elif dom['eventId'] == int(eventindex) : #event['eventId'] :
					domindex += 1
				else :
					break
				#print("in dom loop 2")
				if dom['impactAngle'] < args.impactrange[0]*3.14/180. or  dom['impactAngle'] > args.impactrange[1]*3.14/180.:
					continue
				if dom['distAboveEndpoint'] < args.trackendpoint :
					continue
				if dom['recoDist'] < args.cherdist[0] or dom['recoDist'] > args.cherdist[1] :
					continue
				if event['icHitsOut']> args.nhits[1] or event['icHitsOut'] < 1 :
					continue; 

				if dom['string'] in DC_Strings :
					if first_dc :
						TotalCharge_DC.Fill(event['totalCharge'],weight)
						Zenith_DC.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),weight)
						RecEnergy_DC.Fill(event['reco/energy'],weight)
						EndPointZ_DC.Fill(event['recoEndPoint/z'],weight)
						BoarderDist_DC.Fill(event['borderDistance'],weight)
						StopLikeRatio_DC.Fill(event['stopLikeRatio'],weight)
						RecoLogL_DC.Fill(event['recoLogL'],weight)
						DirectHits_DC.Fill(event['directHits'],weight)
						HitsOut_DC.Fill(event['dcHitsOut'],weight)
						TotalCharge_vs_Zenith_dc.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['totalCharge'],weight)
						StoppingZ_vs_Zenith_dc.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['recoEndPoint/z'],weight)
						BorderDist_vs_Zenith_dc.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['borderDistance'],weight)
						NChannel_vs_Zenith_dc.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['icHitsOut']+event['icHitsIn'],weight)
						StopLike_vs_Zenith_dc.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['stopLikeRatio'],weight)
						rLogL_vs_Zenith_dc.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['recoLogL'],weight)
						first_dc = False
                                        if 'truthDist' in domtable.colnames:
						TruthVsRecoDistance_DC.Fill(dom['recoDist'],dom['truthDist'],weight*event['totalCharge'])
                                        if 'totalChargeMC' in domtable.colnames :
                                                if event['totalCharge'] > 0.0 :
                                                    TruePE_Vs_MeasuredPE_DC.Fill(dom['recoDist'],dom['totalChargeMC']/event['totalCharge'],weight*event['totalCharge'])
					ImpactAll_dc.Fill(ROOT.TMath.Cos(dom['impactAngle']),weight)
					Impact_vs_Zenith_dc.Fill(ROOT.TMath.Cos(dom['impactAngle']),ROOT.TMath.Cos(event['reco/dir/zenith']),weight)
					Impact_seeMPE_dc.Fill(ROOT.TMath.Cos(dom['impactAngle']),weight)
				if dom['string'] in IC_Strings :
					if first_ic :
						TotalCharge_IC.Fill(event['totalCharge'],weight)
						Zenith_IC.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),weight)
						RecEnergy_IC.Fill(event['reco/energy'],weight)
						EndPointZ_IC.Fill(event['recoEndPoint/z'],weight)
						BoarderDist_IC.Fill(event['borderDistance'],weight)
						StopLikeRatio_IC.Fill(event['stopLikeRatio'],weight)
						RecoLogL_IC.Fill(event['recoLogL'],weight)
						DirectHits_IC.Fill(event['directHits'],weight)
						HitsOut_IC.Fill(event['dcHitsOut'],weight)
						TotalCharge_vs_Zenith_ic.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['totalCharge'],weight)
						TotalCharge_vs_Zenith_ic.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['totalCharge'],weight)
						StoppingZ_vs_Zenith_ic.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['recoEndPoint/z'],weight)
						BorderDist_vs_Zenith_ic.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['borderDistance'],weight)
						NChannel_vs_Zenith_ic.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['icHitsOut']+event['icHitsIn'],weight)
						StopLike_vs_Zenith_ic.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['stopLikeRatio'],weight)
						rLogL_vs_Zenith_ic.Fill(ROOT.TMath.Cos(event['reco/dir/zenith']),event['recoLogL'],weight)
						first_ic = False
					if 'truthDist' in domtable.colnames:
                                                TruthVsRecoDistance_IC.Fill(dom['recoDist'],dom['truthDist'],weight*event['totalCharge'])
                                        if 'totalChargeMC' in domtable.colnames :
                                                if event['totalCharge'] > 0.0 :
                                                    TruePE_Vs_MeasuredPE_IC.Fill(dom['recoDist'],dom['totalChargeMC']/event['totalCharge'],weight*event['totalCharge'])
					Impact_vs_Zenith_ic.Fill(ROOT.TMath.Cos(dom['impactAngle']),ROOT.TMath.Cos(event['reco/dir/zenith']),weight)
					ImpactAll_ic.Fill(ROOT.TMath.Cos(dom['impactAngle']),weight)
					if dom['totalCharge'] > 0.0 :
						Impact_seeMPE_ic.Fill(ROOT.TMath.Cos(dom['impactAngle']),weight)

		h5file.close()
	print("nevents = "+str(TotalCharge_IC.Integral()))
	fout.cd()

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
	StoppingZ_vs_Zenith_dc.Write()
	StoppingZ_vs_Zenith_ic.Write()
	BorderDist_vs_Zenith_dc.Write()
	BorderDist_vs_Zenith_ic.Write()
	NChannel_vs_Zenith_dc.Write()
	NChannel_vs_Zenith_ic.Write()
	StopLike_vs_Zenith_dc.Write()
	StopLike_vs_Zenith_ic.Write()
	rLogL_vs_Zenith_dc.Write()
	rLogL_vs_Zenith_ic.Write()
	zenith_all.Write()
	zenith_zenith.Write()
	zenith_endpointz.Write()
	zenith_boarder.Write()
	zenith_stopratio.Write()
	zenith_recoLogL.Write()
	zenith_directHits.Write()
	zenith_all_line.Write()
	zenith_all_spe.Write()
	zenith_all_mpe.Write()
	TruthVsRecoDistance_IC.Write()
	TruthVsRecoDistance_DC.Write()
	TruePE_Vs_MeasuredPE_IC.Write()
	TruePE_Vs_MeasuredPE_DC.Write()
	fout.Close()


