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

def fitconst(mindist,maxdist) :
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
    sigma = 0.0
    npoints = 0
    for ratio in currentratio:
        if ratio['meandistance'] < mindist or ratio['meandistance'] > maxdist :
            continue
        sigma += (m.values["ratio"] - ratio['meancharge'])**2.0
        npoints += 1
    sigma = np.sqrt(sigma)/npoints
    return m.values["ratio"], sigma

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', help='Directory of data files.',type=str,nargs = '+',
			default = ['/data/user/sanchezh/IC86_2015/Final_Level2_IC86_MPEFit_*.h5'])
parser.add_argument('-e', '--eff', help='Ordered list of efficiency simulations to use, 0.9,1.0,1.1,1.2', type = str,
			nargs = '+', default =["",""])
parser.add_argument('-r','--roundtrip',type = str,default="")
parser.add_argument('-o','--out',help='output specifier', type = str, default="")
parser.add_argument('-l','--dist',help='distance range',type = float, nargs = "+",default=[55,165])
args = parser.parse_args()

mindist = args.dist[0]
maxdist = args.dist[1]

Eff1_ic_x = []
Eff1_ic_xerr = []
Eff1_ic_y = []
Eff1_ic_yerr = []
Eff1_dc_x = []
Eff1_dc_xerr = []
Eff1_dc_y = []
Eff1_dc_yerr = []

Eff2_ic_x = []
Eff2_ic_xerr = []
Eff2_ic_y = []
Eff2_ic_yerr = []
Eff2_dc_x = []
Eff2_dc_xerr = []
Eff2_dc_y = []
Eff2_dc_yerr = []

datafile_1 = open_file(args.data[0], mode="r")
datafile_2 = open_file(args.data[1], mode="r")
eff_file_1 = open_file(args.eff[0], mode="r")
eff_file_2 = open_file(args.eff[1], mode="r")
data1_ic = []
data1_dc = []
data2_ic = []
data2_dc = []

for x in datafile_1.root.icecube.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    data1_ic.append(template)
for x in datafile_1.root.deepcore.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    data1_dc.append(template)

for x in datafile_2.root.icecube.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    data2_ic.append(template)
for x in datafile_2.root.deepcore.iterrows() :
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    data2_dc.append(template)

eff1_ic = []
eff1_ic_norm = []
eff1_dc = []
eff1_dc_norm = []
for x in eff_file_1.root.icecube.iterrows():
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    eff1_ic.append(template)
    charge, error = XinterpolationsandYequivError(x['meandistance'],data1_ic)
    template['meancharge']=x['meancharge']/charge
    template['sigmacharge']=(x['meancharge']/charge)*((x['sigmacharge']/x['meancharge'])**2.0+(error/charge)**2.0)**0.5
    Eff1_ic_x.append(template['meandistance'])
    Eff1_ic_y.append(template['meancharge'])
    Eff1_ic_xerr.append(0.0)
    Eff1_ic_yerr.append(template['sigmacharge'])
    eff1_ic_norm.append(template)
    N=0
for x in eff_file_1.root.deepcore.iterrows():
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    eff1_dc.append(template)
    charge, error = XinterpolationsandYequivError(x['meandistance'],data1_dc)
    template['meancharge']=x['meancharge']/charge
    template['sigmacharge']=(x['meancharge']/charge)*((x['sigmacharge']/x['meancharge'])**2.0+(error/charge)**2.0)**0.5
    Eff1_dc_x.append(template['meandistance'])
    Eff1_dc_y.append(template['meancharge'])
    Eff1_dc_xerr.append(0.0)
    Eff1_dc_yerr.append(template['sigmacharge'])
    eff1_dc_norm.append(template)

eff2_ic = []
eff2_ic_norm = []
eff2_dc = []
eff2_dc_norm = []
for x in eff_file_2.root.icecube.iterrows():
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    eff2_ic.append(template)
    charge, error = XinterpolationsandYequivError(x['meandistance'],data2_ic)
    template['meancharge']=x['meancharge']/charge
    template['sigmacharge']=(x['meancharge']/charge)*((x['sigmacharge']/x['meancharge'])**2.0+(error/charge)**2.0)**0.5
    Eff2_ic_x.append(template['meandistance'])
    Eff2_ic_y.append(template['meancharge'])
    Eff2_ic_xerr.append(0.0)
    Eff2_ic_yerr.append(template['sigmacharge'])
    eff2_ic_norm.append(template)
    N=0
for x in eff_file_2.root.deepcore.iterrows():
    template = {'meandistance':0.0,'sigmadistance':0.0,'meancharge': 0.0,'sigmacharge':0.0}
    template['meandistance']=x['meandistance']
    template['sigmadistance']=x['sigmadistance']
    template['meancharge']=x['meancharge']
    template['sigmacharge']=x['sigmacharge']
    eff2_dc.append(template)
    charge, error = XinterpolationsandYequivError(x['meandistance'],data2_dc)
    template['meancharge']=x['meancharge']/charge
    template['sigmacharge']=(x['meancharge']/charge)*((x['sigmacharge']/x['meancharge'])**2.0+(error/charge)**2.0)**0.5
    Eff2_dc_x.append(template['meandistance'])
    Eff2_dc_y.append(template['meancharge'])
    Eff2_dc_xerr.append(0.0)
    Eff2_dc_yerr.append(template['sigmacharge'])
    eff2_dc_norm.append(template)

#compute the DOM efficency for IceCube by fitting the average scaled charge and then fitting line and taking intercept.
currentratio = eff1_ic_norm
currentdata = data1_ic
ratio1_ic,sigma1_ic=fitconst(mindist,maxdist) 
currentratio = eff2_ic_norm
currentdata = data2_ic
ratio2_ic,sigma2_ic=fitconst(mindist,maxdist)

currentratio = eff1_dc_norm
currentdata = data1_dc
ratio1_dc,sigma1_dc=fitconst(mindist,maxdist)
currentratio = eff2_dc_norm
currentdata = data2_dc
ratio2_dc,sigma2_dc=fitconst(mindist,maxdist)

_xrange = []
i=0
while 20.*(i+1) < maxdist :
	if 20.*(i+1) > mindist :
		_xrange.append(20.*(i+1))
	i+=1

fig, ax = plt.subplots()
ax.errorbar(np.array(Eff1_ic_x),np.array(Eff1_ic_y),yerr=np.array(Eff1_ic_yerr),color='C4',marker='o',label='Round Trip Data')
fitline1_ic_y = np.ones(len(_xrange))*ratio1_ic
ax.plot(np.array(_xrange),np.ones(len(_xrange))*ratio1_ic,color='C4')
ax.fill_between(np.array(_xrange),fitline1_ic_y-sigma1_ic,fitline1_ic_y+sigma1_ic,alpha=0.2,color='C4')
ax.errorbar(np.array(Eff2_ic_x),np.array(Eff2_ic_y),yerr=np.array(Eff2_ic_yerr),color='C0',marker='o',label='Systematic Data Set')
fitline2_ic_y = np.ones(len(_xrange))*ratio2_ic
ax.plot(np.array(_xrange),np.ones(len(_xrange))*ratio2_ic,color='C0')
ax.fill_between(np.array(_xrange),fitline2_ic_y-sigma2_ic,fitline2_ic_y+sigma2_ic,alpha=0.2,color='C0')
legend = ax.legend(loc = "lower right",prop=legendfont,ncol=2)
plt.ylim([0.9,1.15])
ax.set_xlabel("Distance from Track (m)",**hfont)
ax.set_ylabel("Q-Ratio (MC/Data)",**hfont)
fig.savefig('IceCube_ChargeRatio_'+args.out+'.png',bbox_inches='tight')

fig, ax = plt.subplots()
ax.errorbar(np.array(Eff1_dc_x),np.array(Eff1_dc_y),yerr=np.array(Eff1_dc_yerr),color='C4',marker='o',label='Round Trip Data')
fitline1_dc_y = np.ones(len(_xrange))*ratio1_dc
ax.plot(np.array(_xrange),np.ones(len(_xrange))*ratio1_dc,color='C4')
ax.fill_between(np.array(_xrange),fitline1_dc_y-sigma1_dc,fitline1_dc_y+sigma1_dc,alpha=0.2,color='C4')
ax.errorbar(np.array(Eff2_dc_x),np.array(Eff2_dc_y),yerr=np.array(Eff2_dc_yerr),color='C0',marker='o',label='Systematic Data Set')
fitline2_dc_y = np.ones(len(_xrange))*ratio2_dc
ax.plot(np.array(_xrange),np.ones(len(_xrange))*ratio2_dc,color='C0')
ax.fill_between(np.array(_xrange),fitline2_dc_y-sigma2_dc,fitline2_dc_y+sigma2_dc,alpha=0.2,color='C0')
legend = ax.legend(loc = "lower right",prop=legendfont,ncol=2)
plt.ylim([0.9,1.15])
ax.set_xlabel("Distance from Track (m)",**hfont)
ax.set_ylabel("Q-Ratio (MC/Data)",**hfont)
fig.savefig('DeepCore_ChargeRatio_'+args.out+'.png',bbox_inches='tight')

