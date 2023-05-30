import numpy as np
from array import array as arr
from tables import open_file
from iminuit import Minuit
import argparse
import matplotlib.pyplot as plt

currenteff = []
currentdata = []
currentratio = []
MCDataRatio = []

mindist = 55.
maxdist = 165.

def GetYequivError(distbin,data) :

	deriv = data[distbin+1].get('meancharge') - data[distbin-1].get('meancharge')
	deriv = deriv / (data[distbin+1].get('meandistance') - data[distbin-1].get('meandistance'))

	yequiv = deriv*data[distbin+1].get('sigmadistance')
	return (data[distbin+1].get('sigmacharge')**2.0 + yequiv**2.0)**0.5

def XinterpolationsandYequivError(dist,data):

    distbin = 0
    while distbin<len(data) and data[distbin].get('meandistance') < dist :
        distbin += 1

    distbin = min(distbin,len(data)-1)

    x_weight = (dist-data[distbin-1].get('meandistance'))/(data[distbin].get('meandistance')-data[distbin-1].get('meandistance'))

    deriv = data[distbin].get('meancharge')-data[distbin-1].get('meancharge')
    deriv = deriv/(data[distbin].get('meandistance')-data[distbin-1].get('meandistance'))

    yequiv = deriv*((1.0-x_weight)*data[distbin-1].get('sigmadistance') + x_weight*data[distbin].get('sigmadistance'))

    intyerror = (1.0-x_weight)*data[distbin-1].get('sigmacharge') + x_weight*data[distbin].get('sigmacharge')

    yerror = (yequiv*yequiv+intyerror*intyerror)**0.5

    charge = (1.0-x_weight)*data[distbin-1].get('meancharge') + x_weight*data[distbin].get('meancharge')

    return charge , yerror

def SimCharge(eff,dist,sim): 
	#Get efficiency curves to interpolate between.
	binmin = int((eff-0.9)/0.1)
	binmax = binmin+1

	# interpolate x  values for lower curve
	chargelow, errorlow   = XinterpolationsandYequivError(dist,sim[binmin])
	chargehigh, errorhigh = XinterpolationsandYequivError(dist,sim[binmax])

	#interpolate between two curves.
	y_weight = (eff-(0.9+0.1*binmin))/0.1

	charge = (1.0-y_weight)*chargelow + y_weight*chargehigh
	chargesig = (1.0-y_weight)*errorlow + y_weight*errorhigh

	return charge, chargesig

def chargedist(dist,amp,base,tau) :
	return base + amp*exp(-dist/tau)

def chargedistChi2(amp,base,tau) :
    chisq = 0.0
    for i in range(len(currentdata)) :
        if currentdata[i]['meandistance'] < mindist or currentdata[i]['meandistance'] > maxdist :
            continue
        data = currentdata[i].get('meancharge')
        error = GetYequivError(i,currentdata)
        exp = chargedist(currentdata[i].get('meandistance'),amp,base,tau)
        chisq += ((data-exp)**2.0)/(error**2.0)
    return chisq

def calcChi2(eff):
	chisq = 0.0
	for i in range(1,len(currentdata)-1) :
		simval ,  simerror = SimCharge(eff,currentdata[i].get('meandistance'),currenteff)
		dataval = currentdata[i].get('meancharge')
		dataerror = GetYequivError(i,currentdata)
		chisq += ((simval-dataval)**2.0)/(dataerror*dataerror+simerror*simerror)
	return chisq

def constChi2(ratio=1.0):
    chisq = 0.0
    for mcratio in currentratio :
        if mcratio['meandistance'] < mindist or mcratio['meandistance'] > maxdist :
            continue
        chisq += ((ratio - mcratio['meancharge'])**2.0)/(mcratio['sigmacharge']**2.0)
    return chisq

def linear(slope, intercept, eff) :
	return intercept + eff*slope

def linearChi2(slope, intercept):
	chisq = 0.0

	for mcdataratio in MCDataRatio :
		chisq += ((mcdataratio['scaledcharge']-linear(slope, intercept, mcdataratio['eff']))**2.0)/(mcdataratio['error']**2.0)
	return chisq


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', help='Directory of data files.',nargs="+",type=str,
			default = ['',''])
parser.add_argument('-e', '--eff', help='Ordered list of efficiency simulations to use, 0.9,1.0,1.1,1.2',nargs="+", type = str, default =["",""])
parser.add_argument('-l','--dist',help='distance range',type = float, nargs = "+",default=[55,165])
args = parser.parse_args()

mindist = args.dist[0]
maxdist = args.dist[1]


Eff_ic_x = []
Eff_ic_xerr = []
Eff_ic_y = []
Eff_ic_yerr = []
Eff_dc_x = []
Eff_dc_xerr = []
Eff_dc_y = []
Eff_dc_yerr = []

datafile_orig = open_file(args.data[0], mode="r")
datafile_new = open_file(args.data[1], mode="r")
mcfile_orig  = open_file(args.eff[0],mode="r")
eff_file = open_file(args.eff[1],mode="r")

data_ic = []
data_dc = []
datanew_ic = []
datanew_dc = []
mcorig_ic = []
mcorig_dc = []

for x in datafile_orig.root.icecube.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    data_ic.append(template)
for x in datafile_orig.root.deepcore.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    data_dc.append(template)

for x in datafile_new.root.icecube.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    datanew_ic.append(template)
for x in datafile_new.root.deepcore.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    datanew_dc.append(template)

for x in mcfile_orig.root.icecube.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    mcorig_ic.append(template)
for x in mcfile_orig.root.deepcore.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    mcorig_dc.append(template)

eff_ic = []
eff_ic_norm = []
eff_dc = []
eff_dc_norm = []

for x in eff_file.root.icecube.iterrows():
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    eff_ic.append(template)
    charge, error = XinterpolationsandYequivError(x['meandistance'],data_ic)
    charge_new, error_new = XinterpolationsandYequivError(x['meandistance'],datanew_ic)
    charge_mc, error_mc = XinterpolationsandYequivError(x['meandistance'],mcorig_ic)
    template['meancharge']=(x['meancharge']/charge_new)/(charge_mc/charge)
    if args.data[0] == args.data[1] :
        template['sigmacharge']=(x['meancharge']/charge_mc)*((x['sigmacharge']/x['meancharge'])**2.0+(error_mc/charge_mc)**2.0)**0.5
    else :
        template['sigmacharge']=(template['meancharge'])*((x['sigmacharge']/x['meancharge'])**2.0+(error/charge)**2.0+(error_mc/charge_mc)**2.0+(error_new/charge_new)**2.0)**0.5
    Eff_ic_x.append(template['meandistance'])
    Eff_ic_y.append(template['meancharge'])
    Eff_ic_xerr.append(0.0)
    Eff_ic_yerr.append(template['sigmacharge'])
    eff_ic_norm.append(template)
N=0
for x in eff_file.root.deepcore.iterrows():
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    eff_dc.append(template)
    charge, error = XinterpolationsandYequivError(x['meandistance'],data_dc)
    charge_new, error_new = XinterpolationsandYequivError(x['meandistance'],datanew_dc)
    charge_mc, error_mc = XinterpolationsandYequivError(x['meandistance'],mcorig_dc)
    template['meancharge']=(x['meancharge']/charge_new)/(charge_mc/charge)
    if args.data[0] == args.data[1] :
        template['sigmacharge']=(x['meancharge']/charge_mc)*((x['sigmacharge']/x['meancharge'])**2.0+(error_mc/charge_mc)**2.0)**0.5
    else :
        template['sigmacharge']=(template['meancharge'])*((x['sigmacharge']/x['meancharge'])**2.0+(error/charge)**2.0+(error_mc/charge_mc)**2.0+(error_new/charge_new)**2.0)**0.5
    Eff_dc_x.append(template['meandistance'])
    Eff_dc_y.append(template['meancharge'])
    Eff_dc_xerr.append(0.0)
    Eff_dc_yerr.append(template['sigmacharge'])
    eff_dc_norm.append(template)

currenteff = eff_ic
currentdata = data_ic

#compute the DOM efficency for IceCube by fitting the average scaled charge and then fitting line and taking intercept.

currentratio = eff_ic_norm
kwargs = dict(ratio=0.8, error_ratio=0.4,limit_ratio=(0.1,1.9),errordef = 1.0)
m = Minuit(constChi2,
            ratio = 0.8,
            error_ratio = 0.4,
            limit_ratio = (0.9,1.20),
            errordef = 1.0
            )
m.migrad()
m.hesse()
print("IC ratio = "+str(m.values["ratio"])+" +/- "+str(m.errors["ratio"]))

#compute the DOM efficency for IceCube by fitting the average scaled charge and then fitting line and taking intercept.
currentratio = eff_dc_norm
m = Minuit(constChi2,
	    ratio=0.8,
	    error_ratio = 0.4,
	    limit_ratio=(0.9,1.2),
	    errordef = 1.0
	  )
m.migrad()
m.hesse()
print("DC ratio = "+str(m.values["ratio"])+" +/- "+str(m.errors["ratio"]))

