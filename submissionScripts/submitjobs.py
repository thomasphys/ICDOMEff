#!/usr/bin/env python
import os, sys
import subprocess
import time

totaljobs = 0
jobcountproc = subprocess.Popen(['condor_q'],stdout=subprocess.PIPE)
jobcountout, error = jobcountproc.communicate()
for line in jobcountout.splitlines() :
	sline = str(line)
	if "Total for tmcelroy" in sline :
		 jobarray = sline.split(",",100)
		 for job in jobarray :
		 	if "idle" in job :
				print(job.split(" ",2))	
				njobs = int(job.split(" ",2)[1])
		 		totaljobs = totaljobs +  njobs 
		 	elif "running" in job :
				print(job.split(" ",2))
				njobs = int(job.split(" ",2)[1])
		 		totaljobs = totaljobs + njobs 

avaliblejobs = 12000-totaljobs

print(avaliblejobs)

#jobcount = []
#submission = []
completejobs  = open("submittedjobs.txt", "a")    
completejobs_read  = open("submittedjobs.txt", "r")
cj = completejobs_read.readlines()

#jobfilelist = ["datarunlist/2015/runlist.txt","mcjobs.txt"]
jobfilelist = ["mcjobs.txt"]
#jobfilelist = ["chargejobs.txt","mcjobs.txt"]
#			   "",
#			   ""] 

totaljobs = 0
jobsrun = 0
for filename in jobfilelist :
	file_object  = open(filename, "r") 	
	fl = file_object.readlines()
	totaljobs += len(fl)
	for line in fl:
		data_array=line.split(":")
		if len(data_array) != 2 : continue
		if data_array[1] in cj : 
			jobsrun += 1
			continue
		#jobcount.append(int(data_array[0]))
		#submission.append(data_array[1])
		if avaliblejobs >= int(data_array[0]):
			print(data_array[1])
			avaliblejobs -= int(data_array[0])
			print(avaliblejobs)
			print(data_array[1].split(" ",1))
			job = data_array[1].replace("\n","")
			subprocess.Popen(job.split(" ",100),stdout=subprocess.PIPE)
			time.sleep(5)
			completejobs.write(data_array[1]+"\n")
	file_object.close()
completejobs_read.close()
completejobs.close()
print(str(jobsrun)+"/"+str(totaljobs))
