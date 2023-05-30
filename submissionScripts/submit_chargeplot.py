#!/usr/bin/env python
import os, sys
import subprocess

opts = {}

opts["gcd"] = sys.argv[1].replace("\n","")
opts["data"] = sys.argv[2].replace("\n","")
opts["out"] = sys.argv[3].replace("\n","")
opts["subdir"] = sys.argv[4].replace("\n","")

subdirectorysplit = opts["subdir"].split("/",100)
pathtosub = opts["out"]+"hd5charge/"
for subdir in subdirectorysplit :
	if subdir == "" or subdir == " " : continue
	if not os.path.exists(pathtosub+subdir) :
		os.mkdir(pathtosub+subdir)
	pathtosub = pathtosub+subdir+"/"

scratch = '/scratch/tmcelroy/domeff'

files_dir = opts["data"]
folderlist = files_dir.split("/",1000)
folder = folderlist[-3] + '_' + folderlist[-2]
#file_list_aux = os.listdir(files_dir)
#file_list = [x for x in file_list_aux if '.i3.zst' in x]

totaljobs = 1000

#totaljobs = len(file_list)
#filecutsuff = file_list[0].replace('.i3.zst', '')
#filenamelist = filecutsuff.split("_",20)
#filenameprefix = filenamelist[0]
#for i in range(len(filenamelist)-2) :
#	filenameprefix = filenameprefix + "_" + filenamelist[i+1]
#filenameprefix = filenameprefix + "_"
filenameprefix = "Level2_IC86-2020_corsika_21269_"
eff_factor = ""
mc = "False"
if "090" in opts["subdir"].split("/",1000)[0] :
	eff_factor = "eff090"
	mc = "True" 
if "100" in opts["subdir"].split("/",1000)[0] :
        eff_factor = "eff100"
	mc = "True"
if "110" in opts["subdir"].split("/",1000)[0] :
        eff_factor = "eff110"
	mc = "True"
if "120" in opts["subdir"].split("/",1000)[0] :
        eff_factor = "eff120"
	mc = "True"
if "0.963" in opts["subdir"].split("/",1000)[0] :
        eff_factor = "eff0.963"
        mc = "True"

if eff_factor == "eff0.963" :
	filenameprefix = "Level2_IC86-2020_corsika_21269_"
elif len(eff_factor) > 0 :
	filenameprefix = "Level2_fixed_"+eff_factor+"_IC86-2020_corsika_21269_"
elif "Run" in opts["subdir"].split("/",1000)[0] :
	filenameprefix = "Level2pass2_IC86.2014_data_Run00126238_Subrun00000000_000" 
print("filenameprefix = "+filenameprefix)

#startnumber = 9999999999
#for file in file_list :
#	filecutsuff_temp = file.replace('.i3.zst', '')
#	filenamelist_temp = filecutsuff_temp.split("_",20)
#	if int(filenamelist_temp[-1])<startnumber :
#		startnumber = int(filenamelist_temp[-1])
startnumber = 0
if len(eff_factor) >0:
	startnumber = int(opts["subdir"].split("/",10)[-1].split("-",10)[0])

job_string = '''#!/bin/bash 

eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/metaprojects/combo/V00-00-04/env-shell.sh /home/tmcelroy/icecube/domeff/process_splineMPE_2015_chargeplot.py -g {} -d {} -r $1 -t .i3.zst -o {} -s {}

'''.format(opts["gcd"],files_dir+filenameprefix,pathtosub+filenameprefix,mc)
procesfilename = 'domeff_process_' + folder + '.sh'
with open('/home/tmcelroy/icecube/domeff/jobscripts/' + procesfilename, 'w') as ofile:
#with open(opts["out"] + '/jobscripts/' + procesfilename, 'w') as ofile:
	ofile.write(job_string)
	subprocess.Popen(['chmod','777','/home/tmcelroy/icecube/domeff/jobscripts/' + procesfilename])

submit_string = '''
executable = /home/tmcelroy/icecube/domeff/jobscripts/{}

transfer_input_files = domanalysis.py,event.py,general.py,geometry.py,process_splineMPE_2015_chargeplot.py,writeChargeInfo.py

#+AccountingGroup="sanctioned.$ENV(USER)"

Arguments = $$([$(Process)+{}])
output = /home/tmcelroy/icecube/domeff/out/DOMeff_process_{}_$(Process).out
error = /home/tmcelroy/icecube/domeff/error/DOMeff_process_{}_$(Process).err
log = /scratch/tmcelroy/domeff/log/DOMeff_process_{}_$(Process).log

Universe = vanilla
request_memory = 4GB
request_cpus = 1

notification = never

+TransferOutput=""

queue {}
'''.format(procesfilename,str(startnumber),str(startnumber),str(startnumber),str(startnumber),str(totaljobs))

submissionfilename = 'domeff_process_' + folder + '.submit'

with open('/home/tmcelroy/icecube/domeff/submissionscripts/' + submissionfilename, 'w') as ofile:
	ofile.write(submit_string)

submit = subprocess.Popen(['condor_submit','/home/tmcelroy/icecube/domeff/submissionscripts/' + submissionfilename])
