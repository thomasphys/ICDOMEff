#!/usr/bin/env python
import os, sys
import subprocess

opts = {}

opts["gcd"] = sys.argv[1]
opts["data"] = sys.argv[2]
opts["run"] = sys.argv[3]
opts["nevents"] = sys.argv[4]
opts["out"] = sys.argv[5]
opts["sim"] = sys.argv[7]
opts["subdir"] = sys.argv[6]

if len(sys.argv) > 8 :
	opts["pre"] = str(float(sys.argv[8])*0.1)
else : 
	opts["pre"] = 1.0

scratch = '/scratch/tmcelroy/domeff'

# check that directorys exist for output
subdirectorysplit = opts["subdir"].split("/",100)
pathtosub = opts["out"]+"datahd5/"
if sys.argv[7] :
	pathtosub = opts["out"]+"hd5/"
for subdir in subdirectorysplit :
	if subdir == "" or subdir == " " : continue
	if not os.path.exists(pathtosub+subdir) :
		os.mkdir(pathtosub+subdir)
	pathtosub = pathtosub+subdir+"/"

files_dir = opts["data"]
folderlist = files_dir.split("/",1000)
folder = folderlist[len(folderlist)-2] + '_' + folderlist[len(folderlist)-1]
file_list_aux = os.listdir(files_dir)
file_list_bz2 = [x for x in file_list_aux if ( '.i3.bz2' in x and '_IT' not in x)]
file_list_gz = [x for x in file_list_aux if ( '.i3.gz' in x and '_IT' not in x)]
file_list_zst = [x for x in file_list_aux if ( '.i3.zst' in x and '_IT' not in x)]

file_list = file_list_zst
ext = '.i3.zst'
if len(file_list_bz2) > len(file_list_gz) and len(file_list_bz2) > len(file_list_zst) :
	file_list = file_list_bz2
	ext = '.i3.bz2'
elif len(file_list_gz) > len(file_list_zst) :
	file_list = file_list_gz
	ext = '.i3.gz'

totaljobs = len(file_list)
filecutsuff = file_list[0].replace('.i3.zst', '')
filenamelist = filecutsuff.split("_",20)
filenameprefix = filenamelist[0]
processfilename = filenamelist[0]
for i in range(1,len(filenamelist)) :
	if "Subrun" in filenamelist[i] :
		if 'Level2pass2' in filecutsuff:
			filenameprefix = filenameprefix + "_" + filenamelist[i]+"_000"
			break
		else :
			filenameprefix = filenameprefix + "_" + "Subrun000"
			break
	filenameprefix = filenameprefix + "_" + filenamelist[i]
	processfilename = processfilename + "_" + filenamelist[i]

startnumber = 9999999999

job_string = '''#!/bin/bash 

eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
#RUNNUM=$( printf '%05d' $1 )
#tar -xvjf {}$RUNNUM.tar.bz2 -C /data/user/tmcelroy/domeff/datatemp
/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/metaprojects/combo/V00-00-04/env-shell.sh /home/tmcelroy/icecube/domeff/process_splineMPE_2015.py -g {} -d {} -r $1 -t {} -o {} -y {}

'''.format(files_dir+"/"+filenameprefix,opts["gcd"],files_dir+filenameprefix,ext,pathtosub+filenameprefix,opts["pre"])
processfilename = processfilename + '.sh'
with open(opts["out"] + '/jobscripts/' + processfilename + '.sh', 'w') as ofile:
	ofile.write(job_string)
	subprocess.Popen(['chmod','777',opts["out"] + '/jobscripts/' + processfilename + '.sh'])

submit_string = '''
executable = {}/jobscripts/{}

transfer_input_files = domanalysis.py,event.py,general.py,geometry.py,process_splineMPE_2015.py,writeEvent.py

Arguments = $(Process)
output = /home/tmcelroy/icecube/domeff/out/{}_$(Process).out
error = /home/tmcelroy/icecube/domeff/error/{}_$(Process).err
log = /scratch/tmcelroy/domeff/log/{}_$(Process).log

Universe = vanilla
request_memory = 4GB
request_cpus = 1

notification = never

+TransferOutput=""

queue {}
'''.format(opts["out"],processfilename + '.sh',processfilename,processfilename,processfilename,str(totaljobs))

submissionfilename = 'domeff_process_' + folder + '.submit'

with open(opts["out"] + '/submissionscripts/' + submissionfilename, 'w') as ofile:
	ofile.write(submit_string)

submit = subprocess.Popen(['condor_submit',opts["out"] + '/submissionscripts/' + submissionfilename])
