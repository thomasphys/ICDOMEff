from collections import namedtuple
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt



Result = namedtuple('Result', ['icvalue','dcvalue'])

def plot_student_results(scores_by_test):
    fig, ax1 = plt.subplots(figsize=(9, 7), constrained_layout=True)

    ax1.set_xlabel('Deviation from Standard Result')

    original_result = Result(1.01147,1.01147*0.963)
    origninal_slope = Result(0.5800588,0.828)


    test_names = list(scores_by_test.keys())
    percentilesic = [(1./origninal_slope.icvalue)*(1.-score.icvalue) for score in scores_by_test.values()]
    percentilesdc = [(1./origninal_slope.dcvalue)*(1.-score.dcvalue) for score in scores_by_test.values()]

    percentilesic[2] = scores_by_test['Distance Near'].icvalue - original_result.icvalue
    percentilesdc[2] = scores_by_test['Distance Near'].dcvalue - original_result.dcvalue
    percentilesic[3] = scores_by_test['Distance Far'].icvalue - original_result.icvalue
    percentilesdc[3] = scores_by_test['Distance Far'].dcvalue - original_result.dcvalue

    rectsdc = ax1.scatter(percentilesdc,test_names,marker = "x",color='black',label="DeepCore")
    rectsic = ax1.scatter(percentilesic,test_names,marker = ".",color='red',label="IceCube")
    ax1.legend()

    ax1.set_xlim([-0.05, .05])
    ax1.set_xticks([-0.05,-0.04, -0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05])
    ax1.xaxis.grid(True, linestyle='--', which='major',
                   color='grey', alpha=.25)
    ax1.axvline(50, color='grey', alpha=0.25)  # median position



    # Set the right-hand Y-axis ticks and labels

scores_by_test = {
    'Energy Range High': Result(1.00528677,0.9866),
    'Energy Range Low': Result(0.9911605,0.991649),
    'Distance Near': Result(1.017300,0.962674),
    'Distance Far': Result(1.009110,0.973606),
    'Zenith Upper': Result(1.02215,1.002268),
    'Zenith Lower': Result(0.978544,1.0123476),
    'CRM Gaisser H3a': Result(1.0001606,1.000162028),
    'CRM Gaisser H4a IT': Result(0.99951304,0.999018),
    'CRM GaisserHillas': Result(0.999974332,0.99956),
    'CRM Hoerandel': Result(0.99700307,0.99596188),
    'CRM Hoerandel 5': Result(0.99770066,0.9974428),
    'CRM Hoerandel IT': Result(0.9970283,0.995966899),
    'Reco MPE':Result(0.99733816,0.9992981667),
    'Reco SPE':Result(0.998247015,0.9999633122),
    #'Reco Linefit':Result(0.950969,0.942375),
    #'DOM Impact Range 0.0-0.25':Result(0.967195,0.914115),
    #'DOM Impact Range 0.25-0.5':Result(0.952859,0.929823),
    #'DOM Impact Range 0.5-0.75':Result(0.965917,0.916817),
    #'DOM Impact Range 0.75-1.0':Result(1.014641,0.877826),
    'Zenith Dist Reweighted':Result(0.9986,0.995014),
    'DOM Impact Reweighted':Result(0.998235,1.00027),
    'Endpoint Z Height Reweighted':Result(0.996597,1.0013769),
    'Endpoint Border Distance Reweighted':Result(0.99513311,0.9913336),
    'Distance From Endpoint Reweighted':Result(0.99823,1.002197),
    'Full Reweight':Result(0.99819,0.99000),
    'Track Reconstruction Bias':Result(0.9713796,0.97377128)
}

ICValue = 1.011724
DCValue = 0.973495
slope_IC = 0.5800588
slope_DC = 0.828

system_by_test = {}
names = ['Energy Range High','Energy Range Low']
for i in range(len(names)) :
    temp = Result((1./slope_IC)*(1.0-scores_by_test[names[i]].icvalue),(1./slope_DC)*(1.0-scores_by_test[names[i]].dcvalue))
    system_by_test[names[i]] = Result(temp.icvalue,temp.dcvalue+(1.3/12.)*temp.icvalue)

#energysigma_ic = (max([system_by_test[name].icvalue for name in names])-min([system_by_test[name].icvalue for name in names]))*0.5
#energysigma_dc = (max([system_by_test[name].dcvalue for name in names])-min([system_by_test[name].dcvalue for name in names]))*0.5

energysigma_ic = abs(system_by_test[names[0]].icvalue)
energysigma_dc = abs(system_by_test[names[0]].dcvalue)

print("Energy Systematic")
print(energysigma_ic)
print(energysigma_dc)


names = ['Distance Near','Distance Far']
for i in range(len(names)):
    temp = Result(ICValue-scores_by_test[names[i]].icvalue,DCValue-scores_by_test[names[i]].dcvalue)
    system_by_test[names[i]] = system_by_test[names[i]] = Result(temp.icvalue,temp.dcvalue+(1.3/12.)*temp.icvalue)

distancesigma_ic = (max([system_by_test[name].icvalue for name in names])-min([system_by_test[name].icvalue for name in names]))*0.5
distancesigma_dc = (max([system_by_test[name].dcvalue for name in names])-min([system_by_test[name].dcvalue for name in names]))*0.5
print('Distance Range')
print(distancesigma_ic)
print(distancesigma_dc)

names = ['Zenith Upper','Zenith Lower']
for i in range(len(names)) :
    temp = Result((1./slope_IC)*(1.0-scores_by_test[names[i]].icvalue),(1./slope_DC)*(1.0-scores_by_test[names[i]].dcvalue))
    system_by_test[names[i]] = Result(temp.icvalue,temp.dcvalue+(1.3/12.)*temp.icvalue)

zenithsigma_ic = (max([system_by_test[name].icvalue for name in names])-min([system_by_test[name].icvalue for name in names]))*0.5
zenithsigma_dc = (max([system_by_test[name].dcvalue for name in names])-min([system_by_test[name].dcvalue for name in names]))*0.5
print('Zenith Range')
print(zenithsigma_ic)
print(zenithsigma_dc)


names = ['CRM Gaisser H3a','CRM Gaisser H4a IT','CRM GaisserHillas','CRM Hoerandel','CRM Hoerandel 5','CRM Hoerandel IT']
for i in range(len(names)) :
    temp = Result((1./slope_IC)*(1.0-scores_by_test[names[i]].icvalue),(1./slope_DC)*(1.0-scores_by_test[names[i]].dcvalue))
    system_by_test[names[i]] = Result(temp.icvalue,temp.dcvalue+(1.3/12.)*temp.icvalue)

modelsigma_ic = (max([system_by_test[name].icvalue for name in names])-min([system_by_test[name].icvalue for name in names]))*0.5
modelsigma_dc = (max([system_by_test[name].dcvalue for name in names])-min([system_by_test[name].dcvalue for name in names]))*0.5
print('Cosmic raymodel')
print(modelsigma_ic)
print(modelsigma_dc)

names = ['Reco MPE','Reco SPE']
for i in range(len(names)) :
    temp = Result((1./slope_IC)*(1.0-scores_by_test[names[i]].icvalue),(1./slope_DC)*(1.0-scores_by_test[names[i]].dcvalue))
    system_by_test[names[i]] = Result(temp.icvalue,temp.dcvalue+(1.3/12.)*temp.icvalue)

recosigma_ic = (max([system_by_test[name].icvalue for name in names])-min([system_by_test[name].icvalue for name in names]))*0.5
recosigma_dc = (max([system_by_test[name].dcvalue for name in names])-min([system_by_test[name].dcvalue for name in names]))*0.5
print('Distance Range')
print(recosigma_ic)
print(recosigma_dc)

names = ['Zenith Dist Reweighted','DOM Impact Reweighted','Endpoint Z Height Reweighted','Endpoint Border Distance Reweighted','Distance From Endpoint Reweighted','Full Reweight']
for i in range(len(names)) :
    temp = Result((1./slope_IC)*(1.0-scores_by_test[names[i]].icvalue),(1./slope_DC)*(1.0-scores_by_test[names[i]].dcvalue))
    system_by_test[names[i]] = Result(temp.icvalue,temp.dcvalue+(1.3/12.)*temp.icvalue)

distributionsigma_ic =  (max([system_by_test[name].icvalue for name in names])-min([system_by_test[name].icvalue for name in names]))*0.5
distributionsigma_dc =  (max([system_by_test[name].dcvalue for name in names])-min([system_by_test[name].dcvalue for name in names]))*0.5
print('distributions')
print(distributionsigma_ic)
print(distributionsigma_dc)

names = ['Track Reconstruction Bias']
for i in range(len(names)) :
    temp = Result((1./slope_IC)*(1.0-scores_by_test[names[i]].icvalue),(1./slope_DC)*(1.0-scores_by_test[names[i]].dcvalue))
    system_by_test[names[i]] = Result(temp.icvalue,temp.dcvalue+(1.3/12.)*temp.icvalue)
trackbiassigma_ic = system_by_test['Track Reconstruction Bias'].icvalue
trackbiassigma_dc = system_by_test['Track Reconstruction Bias'].dcvalue

print('total')
print(np.sqrt(energysigma_ic**2.0+distancesigma_ic**2.0+zenithsigma_ic**2.0+modelsigma_ic**2.0+recosigma_ic**2.0+distributionsigma_ic**2.0+trackbiassigma_ic**2.0))
print(np.sqrt(energysigma_dc**2.0+distancesigma_dc**2.0+zenithsigma_dc**2.0+modelsigma_dc**2.0+recosigma_dc**2.0+distributionsigma_dc**2.0+trackbiassigma_dc**2.0))

names = ['Zenith Dist Reweighted','DOM Impact Reweighted','Endpoint Z Height Reweighted','Endpoint Border Distance Reweighted','Full Reweight','Reco MPE','Reco SPE','CRM Gaisser H3a','CRM Gaisser H4a IT','CRM GaisserHillas','CRM Hoerandel','CRM Hoerandel 5','CRM Hoerandel IT','Zenith Upper','Zenith Lower','Distance Near','Distance Far','Energy Range High','Energy Range Low']
print(np.sqrt(sum([abs(system_by_test[name].icvalue)**2.0 for name in names])))
print(np.sqrt(sum([abs(system_by_test[name].dcvalue)**2.0 for name in names])))

plot_student_results(scores_by_test)
plt.savefig("ICSystematicResults.png",bbox_inches='tight')
