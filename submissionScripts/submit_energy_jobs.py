#!/usr/bin/env python
import os, sys
import subprocess

opts = {}
opts["out"] = "/data/user/tmcelroy/domeff/"

basefolderlist = ["/data/user/tmcelroy/domeff/hd5/RDE_hqe_0.963/"]

energyedges = [0.0,30.0,999999.0]
energyname = ["high","low"]

for i in range(len(energyname)) :
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

/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/metaprojects/combo/V00-00-04/env-shell.sh python /home/tmcelroy/icecube/domeff/ProcessDomInfo.py -d {} -o {} -f GaisserH4a -p {} {}

'''.format(basefolder+folder,basefolder+"/"+energyname[i]+"energy/"+folder,energyedges[i],energyedges[i+1])
			procesfilename = 'energyjob_'+energyname[i]+folder+"_"+'.sh'
			with open(opts["out"] + '/jobscripts/' + procesfilename, 'w') as ofile:
				ofile.write(job_string)
			subprocess.Popen(['chmod','777',opts["out"] + 'jobscripts/' + procesfilename])

			submit_string = '''
executable = {}/jobscripts/{}

output = /home/tmcelroy/icecube/domeff/out/Energy.out
error = /home/tmcelroy/icecube/domeff/error/Energy.err
log = /scratch/tmcelroy/domeff/log/Energy.log

Universe = vanilla
request_memory = 4GB
request_cpus = 1

notification = never

+TransferOutput=""

queue
'''.format(opts["out"],procesfilename)

			submissionfilename = 'energy_process_'+energyname[i]+ folder + '.submit'

			with open(opts["out"] + '/submissionscripts/' + submissionfilename, 'w') as ofile:
				ofile.write(submit_string)
			submit = subprocess.Popen(['condor_submit','/data/user/tmcelroy/domeff/submissionscripts/' + submissionfilename])
