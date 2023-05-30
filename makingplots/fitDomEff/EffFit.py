exec(open("/home/tmcelroy/icecube/domeff/CBcm.py").read())
import numpy as np
from array import array as arr
from tables import open_file
from iminuit import Minuit
import argparse
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.font_manager as font_manager

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

hfont = {'fontname':'serif'}

legendfont = font_manager.FontProperties(family='serif',
                                   weight='normal',
                                   style='normal', size=12)

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
parser.add_argument('-d', '--data', help='Directory of data files.',type=str,
			default = '/data/user/sanchezh/IC86_2015/Final_Level2_IC86_MPEFit_*.h5')
parser.add_argument('-e', '--eff', help='Ordered list of efficiency simulations to use, 0.9,1.0,1.1,1.2', type = str,
			nargs = '+', default =["","","",""])
parser.add_argument('-r','--roundtrip',type = str,default="")
parser.add_argument('-o','--out',help='output specifier', type = str, default="")
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
roundtrip_ic_x = []
roundtrip_ic_xerr = []
roundtrip_ic_y = []
roundtrip_ic_yerr = []
roundtrip_dc_x = []
roundtrip_dc_xerr = []
roundtrip_dc_y = []
roundtrip_dc_yerr = []
ChargeRatio_vs_Eff_ic_x = [0.9,1.0,1.1,1.2]
ChargeRatio_vs_Eff_ic_xerr = [0.0,0.0,0.0,0.0]
#ChargeRatio_vs_Eff_ic_x = [0.9,1.0,1.1,1.2]
#ChargeRatio_vs_Eff_ic_xerr = [0.0,0.0,0.0,0.0]
ChargeRatio_vs_Eff_ic_y = []
ChargeRatio_vs_Eff_ic_yerr = []
ChargeRatio_vs_Eff_dc_x = [0.9,1.0,1.1,1.2]
ChargeRatio_vs_Eff_dc_xerr = [0.0,0.0,0.0,0.0]
#ChargeRatio_vs_Eff_ic_x = [0.9,1.0,1.1,1.2]
#ChargeRatio_vs_Eff_ic_xerr = [0.0,0.0,0.0,0.0]
ChargeRatio_vs_Eff_dc_y = []
ChargeRatio_vs_Eff_dc_yerr = []

datafile = open_file(args.data, mode="r")
eff_file = [open_file(args.eff[0], mode="r"),open_file(args.eff[1], mode="r"),open_file(args.eff[2], mode="r"),open_file(args.eff[3], mode="r")]
data_ic = []
data_dc = []

for x in datafile.root.icecube.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    data_ic.append(template)
for x in datafile.root.deepcore.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    data_dc.append(template)

plot_roundtrip = False
IC_RoundTrip_eff = 1.076
DC_RoundTrip_eff = 1.076*0.963
roundtrip_ic = []
roundtrip_ic_norm = []
roundtrip_dc = []
roundtrip_dc_norm = []
roundtripfile = None
if len(args.roundtrip) > 0 :
    plot_roundtrip = True
    roundtripfile = open_file(args.roundtrip,mode="r")
    for x in roundtripfile.root.icecube.iterrows():
        print(x)
        template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
        template['meandistance']=x['meandistance']
        template['sigmadistance']=x['sigmadistance']
        template['meancharge']=x['meancharge']
        template['sigmacharge']=x['sigmacharge']
        roundtrip_ic.append(template)
        charge, error = XinterpolationsandYequivError(x['meandistance'],data_ic)
        template['meancharge']=x['meancharge']/charge
        template['sigmacharge']=(x['meancharge']/charge)*((x['sigmacharge']/x['meancharge'])**2.0+(error/charge)**2.0)**0.5
        roundtrip_ic_x.append(template['meandistance'])
        roundtrip_ic_y.append(template['meancharge'])
        roundtrip_ic_xerr.append(0.0)
        roundtrip_ic_yerr.append(template['sigmacharge'])
        roundtrip_ic_norm.append(template)
    N=0
    for x in roundtripfile.root.deepcore.iterrows():
        template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
        template['meandistance']=x['meandistance']
        template['sigmadistance']=x['sigmadistance']
        template['meancharge']=x['meancharge']
        template['sigmacharge']=x['sigmacharge']
        roundtrip_dc.append(template)
        charge, error = XinterpolationsandYequivError(x['meandistance'],data_dc)
        template['meancharge']=x['meancharge']/charge
        template['sigmacharge']=(x['meancharge']/charge)*((x['sigmacharge']/x['meancharge'])**2.0+(error/charge)**2.0)**0.5
        roundtrip_dc_x.append(template['meandistance'])
        roundtrip_dc_y.append(template['meancharge'])
        roundtrip_dc_xerr.append(0.0)
        roundtrip_dc_yerr.append(template['sigmacharge'])
        roundtrip_dc_norm.append(template)

eff_ic = []
eff_ic_norm = []
eff_dc = []
eff_dc_norm = []
for file in eff_file :
    eff_ic.append([])
    eff_dc.append([])
    eff_ic_norm.append([])
    eff_dc_norm.append([])
    Eff_ic_x.append([])
    Eff_dc_x.append([])
    Eff_ic_xerr.append([])
    Eff_dc_xerr.append([])
    Eff_ic_y.append([])
    Eff_dc_y.append([])
    Eff_ic_yerr.append([])
    Eff_dc_yerr.append([])
    for x in file.root.icecube.iterrows():
        template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
        template['meandistance']=x['meandistance']
        template['sigmadistance']=x['sigmadistance']
        template['meancharge']=x['meancharge']
        template['sigmacharge']=x['sigmacharge']
        eff_ic[-1].append(template)
        charge, error = XinterpolationsandYequivError(x['meandistance'],data_ic)
        template['meancharge']=x['meancharge']/charge
        template['sigmacharge']=(x['meancharge']/charge)*((x['sigmacharge']/x['meancharge'])**2.0+(error/charge)**2.0)**0.5
        Eff_ic_x[-1].append(template['meandistance'])
        Eff_ic_y[-1].append(template['meancharge'])
        Eff_ic_xerr[-1].append(0.0)
        Eff_ic_yerr[-1].append(template['sigmacharge'])
        eff_ic_norm[-1].append(template)
    N=0
    for x in file.root.deepcore.iterrows():
        template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
        template['meandistance']=x['meandistance']
        template['sigmadistance']=x['sigmadistance']
        template['meancharge']=x['meancharge']
        template['sigmacharge']=x['sigmacharge']
        eff_dc[-1].append(template)
        charge, error = XinterpolationsandYequivError(x['meandistance'],data_dc)
        template['meancharge']=x['meancharge']/charge
        template['sigmacharge']=(x['meancharge']/charge)*((x['sigmacharge']/x['meancharge'])**2.0+(error/charge)**2.0)**0.5
        Eff_dc_x[-1].append(template['meandistance'])
        Eff_dc_y[-1].append(template['meancharge'])
        Eff_dc_xerr[-1].append(0.0)
        Eff_dc_yerr[-1].append(template['sigmacharge'])
        eff_dc_norm[-1].append(template)

currenteff = eff_ic
currentdata = data_ic

#compute the DOM efficency for IceCube by fitting the average scaled charge and then fitting line and taking intercept.
eff_value = 0.9
nominal_eff = 0.94
for eff in eff_ic_norm :
    currentratio = eff
    kwargs = dict(ratio=0.8, error_ratio=0.4,limit_ratio=(0.1,1.9),errordef = 1.0)
    #m = Minuit(constChi2,**kwargs)
    m = Minuit(constChi2,
               ratio = 0.8,
               error_ratio = 0.4,
               limit_ratio = (0.9,1.20),
               errordef = 1.0
               )
    m.migrad()
    m.hesse()
    #m.minos()
    #MCDataRatio.append({'eff':eff_value*nominal_eff,'scaledcharge':m.values["ratio"],'error':m.errors["ratio"]})
    ChargeRatio_vs_Eff_ic_y.append(m.values["ratio"])
    #use sigma as error
    sigma = 0.0
    npoints = 0
    for ratio in currentratio:
        if ratio['meandistance'] < mindist or ratio['meandistance'] > maxdist :
            continue
        sigma += (ChargeRatio_vs_Eff_ic_y[-1] - ratio['meancharge'])**2.0
        npoints += 1
    sigma = np.sqrt(sigma/npoints)
    MCDataRatio.append({'eff':eff_value*nominal_eff,'scaledcharge':m.values["ratio"],'error':sigma})
    #ChargeRatio_vs_Eff_ic_yerr.append(m.errors["ratio"])
    ChargeRatio_vs_Eff_ic_yerr.append(sigma)

    eff_value += 0.1

ChargeRatio_roundtrip_ic_y = 0.0
ChargeRatio_roundtrip_dc_y = 0.0
ChargeRatio_roundtrip_ic_yerr = 0.0
ChargeRatio_roundtrip_dc_yerr = 0.0
if plot_roundtrip:
    currentratio = roundtrip_ic_norm
    m = Minuit(constChi2,
               ratio = 1.0,
               error_ratio = 0.4,
               limit_ratio = (0.9,1.20),
               errordef = 1.0
               )
    m.migrad()
    m.hesse()
    ChargeRatio_roundtrip_ic_y = m.values["ratio"]
    sigma = 0.0
    npoints = 0
    for ratio in currentratio:
        if ratio['meandistance'] < mindist or ratio['meandistance'] > maxdist :
            continue
        sigma += (ChargeRatio_roundtrip_ic_y - ratio['meancharge'])**2.0
        npoints += 1
    sigma = np.sqrt(sigma)/npoints
    ChargeRatio_roundtrip_ic_yerr = sigma

    currentratio = roundtrip_dc_norm
    m = Minuit(constChi2,
               ratio = 1.0,
               error_ratio = 0.4,
               limit_ratio = (0.9,1.20),
               errordef = 1.0
               )
    m.migrad()
    m.hesse()
    ChargeRatio_roundtrip_dc_y = m.values["ratio"]
    #use sigma as error
    sigma = 0.0
    npoints = 0
    for ratio in currentratio:
        if ratio['meandistance'] < mindist or ratio['meandistance'] > maxdist :
            continue
        sigma += (ChargeRatio_roundtrip_dc_y - ratio['meancharge'])**2.0
        npoints += 1
    sigma = np.sqrt(sigma)/npoints
    #ChargeRatio_roundtrip_dc_yerr = m.errors["ratio"]
    ChargeRatio_roundtrip_dc_yerr = sigma

#kwargs = dict(slope=1.0, error_slope=0.1, inter=0.0, error_inter=0.1)
#m = Minuit(linearChi2,**kwargs)
m = Minuit(linearChi2,
			slope = 1.0,
			error_slope=0.1,
			intercept=0.0,
			error_intercept=0.1,
			errordef = 1.0
			)
print("Fit IceCube DOM Efficiency from scaled charge")
m.migrad()
slope_ic = m.values["slope"]
inter_ic = m.values["intercept"]
print("DOM efficiency: %f" % ((1.0-inter_ic)/slope_ic))
m.hesse()
#m.minos()
e_slope_ic = m.errors["slope"]
e_inter_ic = m.errors["intercept"]
nominal_ic_chi = linearChi2(slope_ic,inter_ic)
Chi2ErrorMatric_ic = list()
min_intercept = 10000.
max_intercept = 0.
ic_intercept_lower = 0.0
ic_intercept_upper = 0.0
ic_slope_lower = 0.0
ic_slope_upper = 0.0

for i in range(1000):
    Chi2ErrorMatric_ic.append(list())
    for j in range(1000):
        intercept = inter_ic+(-1.+(float(j)/500.))*e_inter_ic
        slope = slope_ic+(-1.+(float(i)/500.))*e_slope_ic
        chi2 = linearChi2(slope,intercept)
        #print("i = "+ str(i) + " j = " + str(j) + " chi original = " + str(nominal_ic_chi) + " new chi = "+str(chi2))
        if (chi2-nominal_ic_chi) < 1.0 :
            Chi2ErrorMatric_ic[-1].append(chi2)
            eff = ((1.0-intercept)/slope)
            if eff < min_intercept :
                ic_intercept_upper = intercept
                ic_slope_upper = slope
                min_intercept = eff
            if eff > max_intercept :
                ic_intercept_lower = intercept
                ic_slope_lower = slope
                max_intercept = eff

        else :
            Chi2ErrorMatric_ic[-1].append(0.0)
print(m.values)
print(m.errors)
print("DOM Efficiency Error: %f" % (((1.0-inter_ic)/slope_ic)*((e_slope_ic/slope_ic)**2.0+(e_inter_ic/(1.0-inter_ic))**2.0)**0.5))
print("DOM Efficiency Error2: %f" % ((max_intercept-min_intercept)*0.5))

min_intercept_IC = min_intercept
max_intercept_IC = max_intercept

#compute the DOM efficency for IceCube by fitting the average scaled charge and then fitting line and taking intercept.
eff_value = 0.9
nominal_eff = 0.94
N=0
MCDataRatio = []
for eff in eff_dc_norm :
    currentratio = eff
	#kwargs = dict(ratio=0.8, error_ratio=0.4,limit_ratio=(0.1,1.9))
	#m = Minuit(constChi2,**kwargs)
    m = Minuit(constChi2,
				ratio=0.8,
				error_ratio = 0.4,
				limit_ratio=(0.9,1.2),
				errordef = 1.0
				)
    m.migrad()
    m.hesse()
	#m.minos()
    #MCDataRatio.append({'eff':eff_value*nominal_eff,'scaledcharge':m.values["ratio"],'error':m.errors["ratio"]})
    ChargeRatio_vs_Eff_dc_y.append(m.values["ratio"])
    sigma = 0.0
    npoints = 0
    for ratio in currentratio:
        if ratio['meandistance'] < mindist or ratio['meandistance'] > maxdist :
            continue
        sigma += (ChargeRatio_vs_Eff_dc_y[-1] - ratio['meancharge'])**2.0
        npoints += 1
    sigma = np.sqrt(sigma)/npoints
    #print("fit error = " + str(m.errors["ratio"]) + " standar dev of mean = " + str(sigma))
    #ChargeRatio_vs_Eff_dc_yerr.append(m.errors["ratio"])
    MCDataRatio.append({'eff':eff_value*nominal_eff,'scaledcharge':m.values["ratio"],'error':sigma})
    ChargeRatio_vs_Eff_dc_yerr.append(sigma)
    eff_value += 0.1

#kwargs = dict(slope=1.0, error_slope=0.1, inter=0.0, error_inter=0.1)
#m = Minuit(linearChi2,**kwargs)
m = Minuit(linearChi2,
			slope = 1.0,
			error_slope = 0.1,
			intercept = 0.0,
			error_intercept = 0.1,
			errordef = 1.0
			)
print("Fit DeepCore DOM Efficiency from scaled charge")
m.migrad()
slope_dc = m.values["slope"]
inter_dc = m.values["intercept"]
print("DOM efficiency: %f" % ((1.0-inter_dc)/slope_dc))
m.hesse()
#m.minos()
e_slope_dc = m.errors["slope"]
e_inter_dc = m.errors["intercept"]

nominal_dc_chi = linearChi2(slope_dc,inter_dc)
Chi2ErrorMatric_dc = list()
min_intercept = 10000.
max_intercept = 0.
dc_intercept_lower = 0.0
dc_intercept_upper = 0.0
dc_slope_lower = 0.0
dc_slope_upper = 0.0

for i in range(1000):
    Chi2ErrorMatric_dc.append(list())
    for j in range(1000):
        intercept = inter_dc+(-1.+(float(j)/500.))*e_inter_dc
        slope = slope_dc+(-1.+(float(i)/500.))*e_slope_dc
        chi2 = linearChi2(slope,intercept)
        #print("i = "+ str(i) + " j = " + str(j) + " chi original = " + str(nominal_dc_chi) + " new chi = "+str(chi2))
        if (chi2-nominal_dc_chi) < 1.0 :
            Chi2ErrorMatric_dc[-1].append(chi2)
            eff = ((1.0-intercept)/slope)
            if eff < min_intercept :
                dc_intercept_upper = intercept
                dc_slope_upper = slope
                min_intercept = eff
            if eff > max_intercept :
                dc_intercept_lower = intercept
                dc_slope_lower = slope
                max_intercept = eff
        else :
            Chi2ErrorMatric_dc[-1].append(0.0)

print(m.values)
print(m.errors)
print("DOM Efficiency Error: %f" % (((1.0-inter_dc)/slope_dc)*((e_slope_dc/slope_dc)**2.0+(e_inter_dc/(1.0-inter_dc))**2.0)**0.5))
print("DOM Efficiency Error2: %f" % ((max_intercept-min_intercept)*0.5))

min_intercept_DC = min_intercept
max_intercept_DC = max_intercept

datafile.close()
for file in eff_file :
	file.close()

_xrange = []
i=0
while 20.*(i+1) < maxdist :
	if 20.*(i+1) > mindist :
		_xrange.append(20.*(i+1))
	i+=1

fig, ax = plt.subplots()
ax.errorbar(np.array(Eff_ic_x[0]),np.array(Eff_ic_y[0]),yerr=np.array(Eff_ic_yerr[0]),color='C4',marker='o',label=r'$\epsilon_{DOM}$ = 0.846')
fitline_090_ic_y = np.ones(len(_xrange))*ChargeRatio_vs_Eff_ic_y[0]
ax.plot(np.array(_xrange),fitline_090_ic_y,color='C4')
_090_error = ChargeRatio_vs_Eff_ic_yerr[0]
ax.fill_between(np.array(_xrange),fitline_090_ic_y-_090_error,fitline_090_ic_y+_090_error,alpha=0.2,color='C4')
ax.errorbar(np.array(Eff_ic_x[1]),np.array(Eff_ic_y[1]),yerr=np.array(Eff_ic_yerr[1]),color='C1',marker='o',label=r'$\epsilon_{DOM}$ = 0.940')
fitline_100_ic_y = np.ones(len(_xrange))*ChargeRatio_vs_Eff_ic_y[1]
ax.plot(np.array(_xrange),fitline_100_ic_y,color='C1')
_100_error = ChargeRatio_vs_Eff_ic_yerr[1]
ax.fill_between(np.array(_xrange),fitline_100_ic_y-_100_error,fitline_100_ic_y+_100_error,alpha=0.2,color='C1')
ax.errorbar(np.array(Eff_ic_x[2]),np.array(Eff_ic_y[2]),yerr=np.array(Eff_ic_yerr[2]),color='C2',marker='o',label=r'$\epsilon_{DOM}$ = 1.034')
fitline_110_ic_y = np.ones(len(_xrange))*ChargeRatio_vs_Eff_ic_y[2]
ax.plot(np.array(_xrange),fitline_110_ic_y,color='C2')
_110_error = ChargeRatio_vs_Eff_ic_yerr[2]
ax.fill_between(np.array(_xrange),fitline_110_ic_y-_110_error,fitline_110_ic_y+_110_error,alpha=0.2,color='C2')
ax.errorbar(np.array(Eff_ic_x[3]),np.array(Eff_ic_y[3]),yerr=np.array(Eff_ic_yerr[3]),color='C3',marker='o',label=r'$\epsilon_{DOM}$ = 1.128')
fitline_120_ic_y = np.ones(len(_xrange))*ChargeRatio_vs_Eff_ic_y[3]
ax.plot(np.array(_xrange),fitline_120_ic_y,color='C3')
_120_error = ChargeRatio_vs_Eff_ic_yerr[3]
ax.fill_between(np.array(_xrange),fitline_120_ic_y-_120_error,fitline_120_ic_y+_120_error,alpha=0.2,color='C3')
if plot_roundtrip:
    ax.errorbar(np.array(roundtrip_ic_x),np.array(roundtrip_ic_y),yerr=np.array(roundtrip_ic_yerr),color='C0',marker='o',label=r'$\epsilon_{DOM}$ = 1.013')
    fitline_roundtrip_ic_y = np.ones(len(_xrange))*ChargeRatio_roundtrip_ic_y
    ax.plot(np.array(_xrange),fitline_roundtrip_ic_y,color='C0')
    _roundtrip_error = ChargeRatio_roundtrip_ic_yerr
    ax.fill_between(np.array(_xrange),fitline_roundtrip_ic_y-_roundtrip_error,fitline_roundtrip_ic_y+_roundtrip_error,alpha=0.2,color='C0')
legend = ax.legend(loc = "lower right",prop=legendfont,ncol=2)
ax.set_xlabel("Distance from Track (m)",**hfont)
ax.set_ylabel("Q-Ratio (MC/Data)",**hfont)
fig.savefig('IceCube_ChargeRatio_'+args.out+'.png',bbox_inches='tight')

fig, ax = plt.subplots()
ax.errorbar(np.array(Eff_dc_x[0]),np.array(Eff_dc_y[0]),yerr=np.array(Eff_dc_yerr[0]),color='C4',marker='o',label="0.9*0.94 DOM Eff")
fitline_090_dc_y = np.ones(len(_xrange))*ChargeRatio_vs_Eff_dc_y[0]
ax.plot(np.array(_xrange),fitline_090_dc_y,color='C4')
_090_error = ChargeRatio_vs_Eff_dc_yerr[0]
ax.fill_between(np.array(_xrange),fitline_090_dc_y-_090_error,fitline_090_dc_y+_090_error,alpha=0.2,color='C4')
ax.errorbar(np.array(Eff_dc_x[1]),np.array(Eff_dc_y[1]),yerr=np.array(Eff_dc_yerr[1]),color='C1',marker='o',label="1.0*0.94 DOM Eff")
fitline_100_dc_y = np.ones(len(_xrange))*ChargeRatio_vs_Eff_dc_y[1]
ax.plot(np.array(_xrange),fitline_100_dc_y,color='C1')
_100_error = ChargeRatio_vs_Eff_dc_yerr[1]
ax.fill_between(np.array(_xrange),fitline_100_dc_y-_100_error,fitline_100_dc_y+_100_error,alpha=0.2,color='C1')
ax.errorbar(np.array(Eff_dc_x[2]),np.array(Eff_dc_y[2]),yerr=np.array(Eff_dc_yerr[2]),color='C2',marker='o',label="1.1*0.94 DOM Eff")
fitline_110_dc_y = np.ones(len(_xrange))*ChargeRatio_vs_Eff_dc_y[2]
ax.plot(np.array(_xrange),fitline_110_dc_y,color='C2')
_110_error = ChargeRatio_vs_Eff_dc_yerr[2]
ax.fill_between(np.array(_xrange),fitline_110_dc_y-_110_error,fitline_110_dc_y+_110_error,alpha=0.2,color='C2')
ax.errorbar(np.array(Eff_dc_x[3]),np.array(Eff_dc_y[3]),yerr=np.array(Eff_dc_yerr[3]),color='C3',marker='o',label="1.2*0.94 DOM Eff")
fitline_120_dc_y = np.ones(len(_xrange))*ChargeRatio_vs_Eff_dc_y[3]
ax.plot(np.array(_xrange),fitline_120_dc_y,color='C3')
_120_error = ChargeRatio_vs_Eff_dc_yerr[3]
ax.fill_between(np.array(_xrange),fitline_120_dc_y-_120_error,fitline_120_dc_y+_120_error,alpha=0.2,color='C3')
if plot_roundtrip:
    ax.errorbar(np.array(roundtrip_dc_x),np.array(roundtrip_dc_y),yerr=np.array(roundtrip_dc_yerr),color='C0',marker='o',label="1.027*0.94 DOM Eff")
    fitline_roundtrip_dc_y = np.ones(len(_xrange))*ChargeRatio_roundtrip_dc_y
    ax.plot(np.array(_xrange),fitline_roundtrip_dc_y,color='C0')
    _roundtrip_error = ChargeRatio_roundtrip_dc_yerr
    ax.fill_between(np.array(_xrange),fitline_roundtrip_dc_y-_roundtrip_error,fitline_roundtrip_dc_y+_roundtrip_error,alpha=0.2,color='C0')
legend = ax.legend(loc = "lower right",prop=legendfont,ncol=2)
ax.set_xlabel("Distance from Track (m)",**hfont)
ax.set_ylabel("Q-Ratio (MC/Data)",**hfont)
fig.savefig('DeepCore_ChargeRatio_'+args.out+'.png',bbox_inches='tight')

_x = np.array([0.9,1.0,1.1,1.2])*nominal_eff
_y = _x*slope_ic+inter_ic
print(e_slope_ic)
print(e_inter_ic)
_y_top_slope = _x*(slope_ic+e_slope_ic)+inter_ic
_y_bot_slope = _x*(slope_ic-e_slope_ic)+inter_ic
_y_top_slope = _x*(slope_ic)+inter_ic-e_inter_ic
_y_bot_slope = _x*(slope_ic)+inter_ic+e_inter_ic

top = []
bot = []
for i  in range(len(_y_top_slope)) :
	top.append(max(_y_top_slope[i],_y_top_slope[i]))
	bot.append(min(_y_bot_slope[i],_y_bot_slope[i]))
#_y_top = np.array(top) 
_y_top = _x*ic_slope_upper + ic_intercept_upper
#_y_bot = np.array(bot)
_y_bot = _x*ic_slope_lower + ic_intercept_lower

_x2 = np.array([0.9*nominal_eff,(1.-inter_ic)/slope_ic])
_y2 = np.array([1.0,1.0])

_x3 = np.array([(1.-inter_ic)/slope_ic,(1.-inter_ic)/slope_ic])
_y3 = np.array([ChargeRatio_vs_Eff_ic_y[0],1.0])

fig, ax = plt.subplots()
ax.plot(_x,_y,color = 'C0')
ax.fill_between(_x,_y_bot,_y_top,alpha=0.2,color = 'C0')
ax.plot(_x2,_y2,color = 'C4')
ax.plot(_x3,_y3,color = 'C4')
_xint = [min_intercept_IC,max_intercept_IC]
_y1int = [1.0,1.0]
_y2int = [ChargeRatio_vs_Eff_ic_y[0],ChargeRatio_vs_Eff_ic_y[0]]
ax.fill_between(_xint,_y1int,_y2int,alpha=0.2,color = 'C4')
ax.errorbar(np.array(ChargeRatio_vs_Eff_ic_x)*nominal_eff,np.array(ChargeRatio_vs_Eff_ic_y),yerr=np.array(ChargeRatio_vs_Eff_ic_yerr),color = "black",marker=',',capsize=3,linestyle="None")
if plot_roundtrip:
    ax.errorbar(np.array([IC_RoundTrip_eff])*nominal_eff,np.array([ChargeRatio_roundtrip_ic_y]),yerr=np.array([ChargeRatio_roundtrip_ic_yerr]),color = "C1",marker=',',capsize=3,linestyle="None")
ax.set_xlabel(r'$\epsilon_{DOM}$',**hfont)
ax.set_ylabel(r'$\mu_{q-ratio}$',**hfont)
fig.savefig('IceCube_eff_'+args.out+'.png',bbox_inches='tight')

_x = np.array([0.9,1.0,1.1,1.2])*nominal_eff
_y = _x*slope_dc+inter_dc
print(e_slope_dc)
print(e_inter_dc)
_y_top_slope = _x*(slope_dc+e_slope_dc)+inter_dc
_y_bot_slope = _x*(slope_dc-e_slope_dc)+inter_dc
_y_top_slope = _x*(slope_dc)+inter_dc-e_inter_dc
_y_bot_slope = _x*(slope_dc)+inter_dc+e_inter_dc

top = []
bot = []
for i  in range(len(_y_top_slope)) :
	top.append(max(_y_top_slope[i],_y_top_slope[i]))
	bot.append(min(_y_bot_slope[i],_y_bot_slope[i]))
#_y_top = np.array(top) 
_y_top = _x*dc_slope_upper + dc_intercept_upper
#_y_bot = np.array(bot)
_y_bot = _x*dc_slope_lower + dc_intercept_lower

_x2 = np.array([0.9*nominal_eff,(1.-inter_dc)/slope_dc])
_y2 = np.array([1.0,1.0])

_x3 = np.array([(1.-inter_dc)/slope_dc,(1.-inter_dc)/slope_dc])
_y3 = np.array([ChargeRatio_vs_Eff_dc_y[0],1.0])

if plot_roundtrip:
    ChargeRatio_vs_Eff_dc_y.append(ChargeRatio_roundtrip_dc_y)
    ChargeRatio_vs_Eff_dc_yerr.append(ChargeRatio_roundtrip_dc_yerr)
    ChargeRatio_vs_Eff_dc_x.append(DC_RoundTrip_eff)
    ChargeRatio_vs_Eff_dc_xerr.append(0.0)

fig, ax = plt.subplots()
ax.plot(_x,_y,color = 'C0')
ax.fill_between(_x,_y_bot,_y_top,alpha=0.2,color = 'C0')
ax.plot(_x2,_y2,color = 'C4')
ax.plot(_x3,_y3,color = 'C4')
_xint = [min_intercept_DC,max_intercept_DC]
_y1int = [1.0,1.0]
_y2int = [ChargeRatio_vs_Eff_dc_y[0],ChargeRatio_vs_Eff_dc_y[0]]
ax.fill_between(_xint,_y1int,_y2int,alpha=0.2,color = 'C4')
ax.errorbar(np.array(ChargeRatio_vs_Eff_dc_x)*nominal_eff,np.array(ChargeRatio_vs_Eff_dc_y),yerr=np.array(ChargeRatio_vs_Eff_dc_yerr),color = "black",marker=',',capsize=3,linestyle="None")
ax.set_xlabel(r'$\epsilon_{DOM}$',**hfont)
ax.set_ylabel(r'$\mu_{q-ratio}$',**hfont)
fig.savefig('DeepCore_eff_'+args.out+'.png',bbox_inches='tight')

