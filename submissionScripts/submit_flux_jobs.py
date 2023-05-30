#!/usr/bin/env python
import os, sys
import subprocess

opts = {}
opts["out"] = "/data/user/tmcelroy/domeff/"

basefolderlist = ["/data/user/tmcelroy/domeff/hd5/RDE_hqe_0.963/"]

#impactedges = [0.0,0.25,0.5,0.75]
#impactname = ["low","med","high"]

flux = ["GaisserH3a", "GaisserH4a_IT", "GaisserHillas", "Hoerandel", "Hoerandel5", "Hoerandel_IT"]
fluxname = ["GaisserH3a","GaisserH4aIT","GaisserHillas","Hoerandel","Hoerandel5","HoerandelIT"]

for i in range(len(fluxname)) :
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

/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/metaprojects/combo/V00-00-04/env-shell.sh python /home/tmcelroy/icecube/domeff/ProcessDomInfo.py -d {} -o {} -f {}

'''.format(basefolder+folder,basefolder+"/"+fluxname[i]+"/"+folder,flux[i])
			procesfilename = 'fluxjob_'+fluxname[i]+folder+"_"+'.sh'
			with open(opts["out"] + '/jobscripts/' + procesfilename, 'w') as ofile:
				ofile.write(job_string)
			subprocess.Popen(['chmod','777',opts["out"] + 'jobscripts/' + procesfilename])

			submit_string = '''
executable = {}/jobscripts/{}

output = /home/tmcelroy/icecube/domeff/out/Flux.out
error = /home/tmcelroy/icecube/domeff/error/Flux.err
log = /scratch/tmcelroy/domeff/log/Flux.log

Universe = vanilla
request_memory = 4GB
request_cpus = 1

notification = never

+TransferOutput=""

queue
'''.format(opts["out"],procesfilename)

			submissionfilename = 'flux_process_'+fluxname[i]+ folder + '.submit'

			with open(opts["out"] + '/submissionscripts/' + submissionfilename, 'w') as ofile:
				ofile.write(submit_string)
			submit = subprocess.Popen(['condor_submit','/data/user/tmcelroy/domeff/submissionscripts/' + submissionfilename])
