#!/usr/bin/env python
import os, sys
import subprocess

opts = {}
opts["out"] = "/data/user/tmcelroy/domeff/"

#basefolderlist = ["/data/user/tmcelroy/domeff/hd5/RDE_hqe_0.963/"]
basefolderlist = ["/data/user/tmcelroy/domeff/hd5/eff100redo/","/data/user/tmcelroy/domeff/hd5/eff090redo/","/data/user/tmcelroy/domeff/hd5/eff110redo/","/data/user/tmcelroy/domeff/hd5/eff120redo/"]


#impactedges = [0.0,0.25,0.5,0.75]
#impactname = ["low","med","high"]

zenithedges = [40,55,70]
zenithname = ["high","low"]
extraname = ["eff100","eff090","eff110","eff120"]

for i in range(len(zenithname)) :
	for j,basefolder in enumerate(basefolderlist) :

		print(basefolder)
		all_files = os.listdir(basefolder)
		folderlist = []
		for folder in all_files:
			if os.path.isdir(basefolder+folder):
				folderlist.append(folder)

		for folder in folderlist :
			if "zenith" in folder :
				continue
			print(folder)
			job_string = '''#!/bin/bash 

eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/metaprojects/combo/V00-00-04/env-shell.sh python /home/tmcelroy/icecube/domeff/ProcessDomInfo.py -d {} -o {} -f GaisserH4a -z {} {}

'''.format(basefolder+folder,basefolder+"/"+zenithname[i]+"zenith/"+folder,zenithedges[i],zenithedges[i+1])
			procesfilename = 'zenithjob_'+extraname[j]+"_"+zenithname[i]+folder+"_"+'.sh'
			with open(opts["out"] + '/jobscripts/' + procesfilename, 'w') as ofile:
				ofile.write(job_string)
			subprocess.Popen(['chmod','777',opts["out"] + 'jobscripts/' + procesfilename])

			submit_string = '''
executable = {}/jobscripts/{}

output = /home/tmcelroy/icecube/domeff/out/Zenith.out
error = /home/tmcelroy/icecube/domeff/error/Zenith.err
log = /scratch/tmcelroy/domeff/log/Zenith.log

Universe = vanilla
request_memory = 4GB
request_cpus = 1

notification = never

+TransferOutput=""

queue
'''.format(opts["out"],procesfilename)

			submissionfilename = 'zenith_process_'+extraname[j]+"_"+zenithname[i]+ folder + '.submit'

			with open(opts["out"] + '/submissionscripts/' + submissionfilename, 'w') as ofile:
				ofile.write(submit_string)
			submit = subprocess.Popen(['condor_submit','/data/user/tmcelroy/domeff/submissionscripts/' + submissionfilename])
