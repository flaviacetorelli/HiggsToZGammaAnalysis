#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='This script splits CMSSW tasks in multiple jobs and sends them on HERCULES')

parser.add_argument("-l", "--label",          required=True,     type=str,  help="job label")
parser.add_argument("-b", "--baseFolder",     required=True,     type=str,  help="base folder")
parser.add_argument("-e", "--exeName",        required=True,     type=str,  help="absolute path of executable")
parser.add_argument("-i", "--inputFolder",    required=True,     type=str,  help="folder where input files are stored")
parser.add_argument("-o", "--outputFolder",   required=True,     type=str,  help="folder where to store output files")
parser.add_argument("-f", "--outputFileName", required=True,     type=str,  help="base name of output files [outputFileName]_i.root")
parser.add_argument("-c", "--configFile",     required=True,     type=str,  help="CMSSW config file to be run")
parser.add_argument("-n", "--nJobs",          required=True,     type=int,  help="number of jobs")
parser.add_argument("-q", "--queue",          default="workday",     type=str,  help="condor queue to use")
parser.add_argument("-s", "--submit",                                       help="submit jobs", action='store_true')
parser.add_argument("-v", "--verbose",                                      help="increase output verbosity", action='store_true')
parser.add_argument("--htcondor",                                       help="use htcondor", action='store_true')

args = parser.parse_args()


#create input file list
command = '/bin/find '+args.inputFolder+' -type f | grep root | grep -v failed'
inputFiles = subprocess.check_output(command,shell=True).splitlines()
nInputFiles = len(inputFiles)
nFilesPerJob = int(math.ceil(1.*nInputFiles/args.nJobs))
if args.verbose:
   print inputFiles
   print nInputFiles
   print args.nJobs
   print nFilesPerJob

print 
print 'START'
print 

currDir = os.getcwd()

print

try:
   subprocess.check_output(['mkdir','jobs'])
except subprocess.CalledProcessError as e:
   print e.output
try:
   subprocess.check_output(['mkdir','jobs/'+args.label])
except subprocess.CalledProcessError as e:
   print e.output


##### loop for creating and sending jobs #####
for x in range(1, args.nJobs+1):
   
   ##### creates directory and file list for job #######
   jobDir = currDir+'/jobs/'+args.label+'/job_'+str(x)
   os.system('mkdir '+jobDir)
   os.chdir(jobDir)
   #os.system("sed '"+str(1+interval*(x-1))+","+str(interval*x)+"!d' ../../"+FileList+" > list.txt ")
   
   ##### creates input file list #######
   with open(jobDir+'/inputFileList_'+str(x)+'.txt', 'w') as list:
      for fileIt in range(0, nFilesPerJob):
         if (x-1)*nFilesPerJob+fileIt < len(inputFiles):
            list.write(inputFiles[(x-1)*nFilesPerJob+fileIt]+'\n')

   ###if args.verbose:
    ###  print inputFileList
   
   ##### creates config file #######
   with open(args.baseFolder+'/'+args.configFile) as fi:
      contents = fi.read()
      #replaced_contents = contents.replace('INPUTFILELIST', jobDir+'/inputFileList_'+str(x)+'.txt')
      #replaced_contents = replaced_contents.replace('OUTPUTFILE', args.outputFolder+'/'+args.outputFileName+'_'+args.label+'_'+str(x))
   with open(jobDir+"/config.cfg", "w") as fo:
      fo.write(contents)
   inputFileListName = jobDir+'/inputFileList_'+str(x)+'.txt'
   command = 'sed -i \"s%^inputFileList .*$%inputFileList '+jobDir+'/inputFileList_'+str(x)+'.txt'+'%\" '+jobDir+'/config.cfg'
   os.system(command)
   command = 'sed -i \"s%^outputFileName .*$%outputFileName '+args.outputFolder+'/'+args.outputFileName+'_'+args.label+'_'+str(x)+'%\" '+jobDir+'/config.cfg'
   os.system(command)
   
   ##### creates jobs #######
   with open('job_'+str(x)+'.sh', 'w') as fout:
      fout.write("#!/bin/sh\n")
      fout.write("echo\n")
      fout.write("echo 'START---------------'\n")
      fout.write("echo 'current dir: ' ${PWD}\n")
      fout.write("cd /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/CMSSW_10_5_0/src\n")
      fout.write("eval `scramv1 runtime -sh`\n")
      fout.write("cd "+str(args.baseFolder)+"\n")
      fout.write("echo 'current dir: ' ${PWD}\n")
      fout.write("source scripts/setup.sh\n")
      fout.write("unbuffer "+args.exeName+" "+jobDir+"/config.cfg > "+jobDir+"/out.log\n")
      fout.write("echo 'STOP---------------'\n")
      fout.write("echo\n")
      fout.write("echo\n")
   os.system("chmod 755 job_"+str(x)+".sh")
   
   ###### sends bjobs ######

   if args.htcondor:
      with open('htcondor_job_'+str(x)+'.sub', 'w') as fout:
         fout.write("executable            = %s/job_%s.sh\n"%(jobDir,str(x)))
         fout.write("arguments             = $(ClusterId) $(ProcId)\n")
         fout.write("output                = %s/job_%s.$(ClusterId).$(ProcId).out\n"%(jobDir,str(x)) )
         fout.write("error                 = %s/job_%s.$(ClusterId).$(ProcId).err\n"%(jobDir,str(x)) )
         fout.write("log                   = %s/job_%s.$(ClusterId).log\n"%(jobDir,str(x)) )
         fout.write("transfer_output_files = \"\"\n" )
         fout.write("+JobFlavour           = \"%s\"\n"%(args.queue) )
         fout.write("queue \n")
      if args.submit:
         command = "condor_submit "+jobDir+" htcondor_job_%s.sub"%(str(x))
         print command
         os.system(command)
         print "htcondor job nr. " + str(x) + " submitted"

      
   os.chdir("../..")
   
print
print "your jobs:"
os.system("condor_q")
print
print 'END'
print
