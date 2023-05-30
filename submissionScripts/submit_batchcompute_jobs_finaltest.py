#!/usr/bin/env python
import os, sys
import subprocess

opts = {}
opts["out"] = "/data/user/tmcelroy/domeff/" #sys.argv[1].replace("\n","")

#basefolderlist = ["/data/user/tmcelroy/domeff/hd5/eff103/"]

#basefolderlist = ["/data/user/tmcelroy/domeff/hd5/eff090redo/","/data/user/tmcelroy/domeff/hd5/eff100redo/","/data/user/tmcelroy/domeff/hd5/eff110redo/"]
#basefolderlist = ["/data/user/tmcelroy/domeff/hd5/RDE/"]
#basefolderlist = ["/data/user/tmcelroy/domeff/hd5/RDE_hqe_0.99/","/data/user/tmcelroy/domeff/hd5/RDE_hqe_1.1/"]
#basefolderlist = ["/data/user/tmcelroy/domeff/hd5/RDE_hqe_0.963/"]
basefolderlist = ["/data/user/tmcelroy/domeff/hd5/RDE_hqe_0.963_shift2/"]

#efflist = ["eff090","eff100","eff110","eff120"]
#efflist = ["final"]
efflist = ["RDE_hqe_0.99","RDE_hqe_1.1"]
#efflist = ["RDE_hqe_0.963"]
for i, basefolder in enumerate(basefolderlist) :

	eff = efflist[i] #basefolder.split("/",100)[-2]

	all_files = os.listdir(basefolder)
	folderlist = []
	for folder in all_files:
		if os.path.isdir(basefolder+folder):
			folderlist.append(folder)

	for folder in folderlist :

		job_string = '''#!/bin/bash 

eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/metaprojects/combo/V00-00-04/env-shell.sh python /home/tmcelroy/icecube/domeff/ProcessDomInfo.py -d {} -e {} -o {} -f GaisserH4a

'''.format(basefolder+folder,eff,basefolder+folder)
		procesfilename = 'domjob_' + eff+folder + '.sh'
		with open(opts["out"] + '/jobscripts/' + procesfilename, 'w') as ofile:
			ofile.write(job_string)
			subprocess.Popen(['chmod','777',opts["out"] + 'jobscripts/' + procesfilename])

		submit_string = '''
executable = {}/jobscripts/{}

output = /home/tmcelroy/icecube/domeff/out/DOMredocheck.out
error = /home/tmcelroy/icecube/domeff/error/DOMredocheck.err
log = /scratch/tmcelroy/domeff/log/DOMredocheck.log

Universe = vanilla
request_memory = 4GB
request_cpus = 1

notification = never

+TransferOutput=""

queue
'''.format(opts["out"],procesfilename)

		submissionfilename = 'domfinal_process_' +eff+ folder + '.submit'

		with open(opts["out"] + '/submissionscripts/' + submissionfilename, 'w') as ofile:
			ofile.write(submit_string)

		submit = subprocess.Popen(['condor_submit','/data/user/tmcelroy/domeff/submissionscripts/' + submissionfilename])
