# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015

import numpy as np

try:
    from randomAtomBox import atombox
    from lammps import lammps
    from SHCPostProc import SHCPostProc 
    from prepareBatch import *
    from errors import *
except(ImportError):
    print "Failed in importing, try deleting .pyc files if you changed between Taito and Sisu."

import os, sys
from subprocess import call, Popen, PIPE

def main(filePrefix):

#    filePrefix='140415c'
    quenchFileInit='quench_asi_nofilename.lmp'
    # quenchFile='quench_asi_'+filePrefix+'.lmpc'
    nemdFileInit='nemd_asi_nofilename.lmp'
    # nemdFile='nemd_asi_'+filePrefix+'.lmpc'
    dataFile=filePrefix+'_Si.dat'
    restartFile=filePrefix+'.quenched.restart'
    LAMMPSParams='-var filename '+filePrefix

    length=200 # Length of the system in Angstroms
    width=70 # Width of the system in Angstroms
    # width=30 # Width of the system in Angstroms
    mass=28.0 # Mass to be used in writing the LAMMPS data file
    rho=2.291 # Density in g/cm^3
    # rho=2.33 # Larkin Si.sw
    Natoms=np.int(np.round(rho*1e-3/1e-6*length*width**2*1e-30/(mass*1.6605e-27)))
    print "The system contains %d atoms." % (Natoms)

    timeQuench='0-12:00:00'
    timeNemd='0-12:00:00'
    # partition='test'
    # partition='parallel'
    partition='small'
    nodes=4

    if os.path.isfile(restartFile):
        answer=raw_input('The restartfile '+restartFile+' exists, do you want to perform the quenching anyway ([y/n], default=y)?.\n')
        # print answer
        if answer=='n':
            runQuench=False
            print 'Skipping the quenching.'
        else:
            runQuench=True
            print 'Performing the quenching.'
    else:
        runQuench=True
        
    # sys.exit()
    # Determine the system from the hostname
    p=Popen(["hostname"],shell=True, stdout=PIPE, stderr=PIPE)
    out,err=p.communicate()
    host=out.split('-')[0]
    print "Running in "+host+"."

    system=host.upper()
    if system=='SISU':
        LammpsExec='/wrk/kisaaski/lammps_sisu/lmp_sisu'
        ncpus=nodes*24
    elif system=='TAITO':
        LammpsExec='/wrk/kisaaski/lammps_sisu/lmp_taito_040215'
        ncpus=nodes*16
    else:
        raise StandardError, "Could not determine the system."

    print "%d atoms per core." % (round(Natoms/ncpus))
    # sys.exit()

    if runQuench:
        # Create the box of silicon atoms, write to datafile
        ab=atombox(length,width,Natoms)
        ab.fillBox(seed=1234)
        ab.writeToFile(dataFile,mass)
        del ab

        batchQuenchScript='batch_quench_'+filePrefix+'.sh'
        # Prepare the batch script
        prepareBatch(batchQuenchScript,LammpsExec,quenchFileInit,execParams=LAMMPSParams,outputFile=filePrefix+'_quench.out',ntasks=ncpus,nodes=nodes,time=timeQuench,partition=partition,system=system,out=filePrefix+'_quench.slurm',jobname=filePrefix+'_quench.job')

        # sys.exit()
        # Submit the batch script
        try:
            command1="sbatch "+batchQuenchScript 
            print "Shell command: " + command1
            p=Popen([command1],shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            out,err=p.communicate()
            p.wait() 
            if err!='': # Failed the submission
                raise SubmissionError("Failed the submission of quench job.",err)
            if out!='':
                print "Output from the first submission: "+ out.strip('\n') 
                quenchJobId=str(out.split(' ')[-1]).strip(' ').strip('\n')
                print "Job id="+str(quenchJobId)
        except SubmissionError as S:
            # if err!='':
            print "Error: ", S.message
            print "Output: " + S.error
            sys.exit()
    else: # Quench not run, set quenchJobId to empty string so that no dependency is created for the NEMD job.
        quenchJobId=''
    # sys.exit()
    # Prepare the second batch script
    batchNemdScript='batch_nemd_'+filePrefix+'.sh'
    # Prepare the NEMD batch script
    prepareBatch(batchNemdScript,LammpsExec,nemdFileInit,execParams=LAMMPSParams,outputFile=filePrefix+'_nemd.out',ntasks=ncpus,time=timeNemd,nodes=nodes,partition=partition,system=system,dependency='afterany:'+quenchJobId,out=filePrefix+'_nemd.slurm',jobname=filePrefix+'_nemd.job')

    # Submit the batch script

    try:
        command2="sbatch "+batchNemdScript
        print "Shell command: " + command2
        # p=Popen([command2],shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        p=Popen([command2],shell=True, stdout=PIPE, stderr=PIPE)
        out,err=p.communicate()
        p.wait() 
        if err!='': # Failed the submission
            raise SubmissionError("Failed the submission of NEMD job.",err)
        if out!='':
            print "Output from the second submission: "+ out.strip('\n') 
        NemdJobId=str(out.split(' ')[-1]).strip(' ').strip('\n')
        # print "Job id="+str(NemdJobId)
    except SubmissionError as S:
        # if err!='':
        print "Error: ", S.message
        print "Error output from the submission: " + S.error
        sys.exit()

if __name__=="__main__":
    import argparse

    parser=argparse.ArgumentParser()
    parser.add_argument("filePrefix",help="The used file prefix.")
    args=parser.parse_args()
    filePrefix=args.filePrefix

    main(filePrefix)
