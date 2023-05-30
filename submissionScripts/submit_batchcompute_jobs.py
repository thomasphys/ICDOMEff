#!/usr/bin/env python
import os, sys
import subprocess

opts = {}
opts["data"] = "/data/user/tmcelroy/domeff/" #sys.argv[1].replace("\n","")
opts["out"] = "/data/user/tmcelroy/domeff/" #sys.argv[2].replace("\n","")
opts["weight"] = 1.0 #sys.argv[3].replace("\n","")

#basefolderlist = ["/data/user/tmcelroy/domeff/hd5/eff090/","/data/user/tmcelroy/domeff/hd5/eff100/","/data/user/tmcelroy/domeff/hd5/eff110/","/data/user/tmcelroy/domeff/hd5/eff120/"]
#basefolderlist = ["/data/user/tmcelroy/domeff/hd5/RDE/"]
#basefolderlist = ["/data/user/tmcelroy/domeff/hd5/RDE_hqe_1.0269/","/data/user/tmcelroy/domeff/hd5/RDE_hqe_0.99/","/data/user/tmcelroy/domeff/hd5/RDE_hqe_1.1/"]
basefolderlist = ["/data/user/tmcelroy/domeff/hd5/RDE_hqe_0.963/"]

for basefolder in basefolderlist :

	eff = "IC86" #basefolder.split("/",100)[-2]

	all_files = os.listdir(basefolder)
	folderlist = []
	for folder in all_files:
		if os.path.isdir(basefolder+folder):
			folderlist.append(folder)

	for folder in folderlist :

		job_string = '''#!/bin/bash 

eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/metaprojects/combo/V00-00-04/env-shell.sh python /home/tmcelroy/icecube/domeff/ProcessDomInfoReweight.py -d {} -e {} -o {} -w {} -f GaisserH4a

'''.format(basefolder+folder,eff,basefolder+folder+"_reweight",opts["weight"])
		procesfilename = 'weightjob_' + eff+folder + '.sh'
		with open(opts["out"] + '/jobscripts/' + procesfilename, 'w') as ofile:
			ofile.write(job_string)
			subprocess.Popen(['chmod','777',opts["out"] + 'jobscripts/' + procesfilename])

		submit_string = '''
executable = {}/jobscripts/{}

output = /home/tmcelroy/icecube/domeff/out/ReweightPlot.out
error = /home/tmcelroy/icecube/domeff/error/ReweightPlot.err
log = /scratch/tmcelroy/domeff/log/ReweightPlot.log

Universe = vanilla
request_memory = 4GB
request_cpus = 1

notification = never

+TransferOutput=""

queue
'''.format(opts["out"],procesfilename)

		submissionfilename = 'weightplot_process_' +eff+ folder + '.submit'

		with open(opts["out"] + '/submissionscripts/' + submissionfilename, 'w') as ofile:
			ofile.write(submit_string)

		submit = subprocess.Popen(['condor_submit','/data/user/tmcelroy/domeff/submissionscripts/' + submissionfilename])
