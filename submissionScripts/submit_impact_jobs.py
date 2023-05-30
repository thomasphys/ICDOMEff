#!/usr/bin/env python
import os, sys
import subprocess
import numpy as np

opts = {}
opts["out"] = "/data/user/tmcelroy/domeff/"

basefolderlist = ["/data/user/tmcelroy/domeff/hd5/RDE_hqe_0.963/"]

impactedges = [np.arccos(0.0)*180./np.pi,np.arccos(0.25)*180./np.pi,np.arccos(0.5)*180./np.pi,np.arccos(0.75)*180./np.pi,np.arccos(1.0)*180./np.pi]
impactname = ["low","med","high","vert"]

#impactedges = [0.75,1.0]
#impactname = ["vert"]

for i in range(len(impactname)) :
	for basefolder in basefolderlist :

		print(basefolder)
		all_files = os.listdir(basefolder)
		folderlist = []
		for folder in all_files:
			if os.path.isdir(basefolder+folder):
				folderlist.append(folder)

		for folder in folderlist :
			print(folder)
			job_string = '''#!/bin/bash 

eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/metaprojects/combo/V00-00-04/env-shell.sh python /home/tmcelroy/icecube/domeff/ProcessDomInfo.py -d {} -o {} -f GaisserH4a -i {} {}

'''.format(basefolder+folder,basefolder+"/impact"+impactname[i]+"/"+folder,impactedges[i],impactedges[i+1])
			procesfilename = 'impactjob_'+impactname[i]+folder+"_"+'.sh'
			with open(opts["out"] + '/jobscripts/' + procesfilename, 'w') as ofile:
				ofile.write(job_string)
			subprocess.Popen(['chmod','777',opts["out"] + 'jobscripts/' + procesfilename])

			submit_string = '''
executable = {}/jobscripts/{}

output = /home/tmcelroy/icecube/domeff/out/Impact.out
error = /home/tmcelroy/icecube/domeff/error/Impact.err
log = /scratch/tmcelroy/domeff/log/Impact.log

Universe = vanilla
request_memory = 4GB
request_cpus = 1

notification = never

+TransferOutput=""

queue
'''.format(opts["out"],procesfilename)

			submissionfilename = 'impact_process_'+impactname[i]+ folder + '.submit'

			with open(opts["out"] + '/submissionscripts/' + submissionfilename, 'w') as ofile:
				ofile.write(submit_string)
			submit = subprocess.Popen(['condor_submit','/data/user/tmcelroy/domeff/submissionscripts/' + submissionfilename])
