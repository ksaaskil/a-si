#!/usr/bin/python
# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015

from __future__ import division
import numpy as np
from SHCPostProc import SHCPostProc

def main(fileprefix,KijFilePrefix=None):

    if KijFilePrefix is None:
        KijFilePrefix=fileprefix
    # fileprefix='140415a'
    # Datafolder 
    dataFolder='DATA'
    outputFolder=dataFolder+'/'+fileprefix+'_tar'
    # Post-processor searches/saves file "KijFilePrefix.Kij.npy"
    # KijFilePrefix='130415a' 
    # KijFilePrefix=fileprefix
    # Use this restart file if the force constant file cannot be found
    restartFile=KijFilePrefix+'.quenched.restart' 

    # Create the data folder
    from subprocess import call
    command=["mkdir","-p",outputFolder]
    print " ".join(command)
    call(command)
    
    dt_md=2.5e-15 # Timestep used in MD, affects the frequency grid
    widthWin=0.1e12 # Width of the Daniell smoothing window in Hz
    # The velocity dump file from LAMMPS
    fileVels=fileprefix+'.vels.dat' 
    # The compactly formatted velocity file, produced using a C++ script if not found
    fileCompactVels=fileprefix+'.vels.dat.compact' 

    # Correct the units, this assumes the unit of eV/(A^2)*(A/ps)^2 for the v_iK_{ij}v_j product (LAMMPS metal units)
    scaleFactor=1.602e-19/(1e-20)*1e4
    
    # Prepare the post-processor
    pP=SHCPostProc(fileCompactVels,KijFilePrefix,
                   dt_md=dt_md,scaleFactor=scaleFactor,
                   LAMMPSDumpFile=fileVels,widthWin=widthWin,
                   LAMMPSRestartFile=restartFile,
                   NChunks=100,chunkSize=50000,
                   backupPrefix=fileprefix,
                   reCalcVels=False,
                   reCalcFC=False)
    # Post-process
    pP.postProcess() # All variables will be contained in the object pP
    
    # Various output options
    
    # Pickling the post-processing instance
    import cPickle as pickle
    with open(outputFolder+'/'+fileprefix+'_PP.pckl','w') as f:
        pickle.dump(pP,f)
    
    # Saving into numpy files
    np.save(outputFolder+'/'+fileprefix+'_oms.npy',pP.oms_fft)
    np.save(outputFolder+'/'+fileprefix+'_SHC.npy',pP.SHC_smooth)

    # Saving to file
    print "Saving to file "+outputFolder+'/'+fileprefix+'_SHC.txt'
    np.savetxt(outputFolder+'/'+fileprefix+'_SHC.txt',np.column_stack((pP.oms_fft,pP.SHC_smooth)))

    # Copy relevant files to the DATA folder
    command="cp "+fileprefix+'.aveinput_*.dat'+" "+outputFolder
    print command
    call(command,shell=True)

    command="cp "+fileprefix+'.Ti*.dat'+" "+outputFolder
    print command
    call(command,shell=True)  

    # Tar the output folder without the DATA folder 
    command=["tar","-czvf",fileprefix+'_tar.tgz','--directory='+dataFolder,outputFolder.strip(dataFolder+'/')]
    print " ".join(command)
    call(command)

    # Plotting if available
    # import matplotlib.pylab as plt
    # plt.plot(pP.oms_fft/(2*np.pi*1.0e12),pP.SHC_smooth)
    # plt.xlabel('Frequency (THz)')
    # plt.ylabel('Spectral current')
    # plt.savefig(fileprefix+'_SHC.eps')

def calcFC_fromFile(datafile):
    
    print "Reading the file %s." % (datafile)
    with open(datafile,'rb') as f:
        
        s1=f.readline().split()
        NL=int(s1[1])
        print 'N_left=%d' % (NL)
        s2=f.readline().split()
        NR=int(s2[1])
        print 'N_right=%d' % (NR)
        
        ids_L=[]
        ids_R=[]
        
        for k in range(NL+NR):
            s=f.readline().split()
            atom_id=int(s[0])
            atom_type=int(s[1])
            if atom_type==6:
                ids_L.append(k)
            elif atom_type==7:
                ids_R.append(k)
            
        print "len(ids_L)=%d, len(ids_R)=%d" % (len(ids_L), len(ids_R))
                
        HSTEP=f.readline().split()[1]
        HSTEP=float(HSTEP)
        print "HSTEP=%f" % (HSTEP)

        Kij=np.zeros((3*NL,3*NR))
        
        for k in range(3*NL):
        # for k in range(3*NL): # Loop over particles
            print "i=%d/%d" % (k,3*NL)
            
            # Positive direction shift
            # Skip ten lines
            for kk in range(9):
                f.readline()
            A=np.fromfile(f,count=4*(NL+NR),sep=" ")
            ids=A[0::4]
            fx1=A[1::4]
            fy1=A[2::4]
            fz1=A[3::4]
            
            # Negative direction shift
            for kk in range(9):
                f.readline()
            A=np.fromfile(f,count=4*(NL+NR),sep=" ")
            ids=A[0::4]
            fx2=A[1::4]
            fy2=A[2::4]
            fz2=A[3::4]
            
            dFx=fx1[ids_R]-fx2[ids_R]
            dFy=fy1[ids_R]-fy2[ids_R]
            dFz=fz1[ids_R]-fz2[ids_R]
            
            Kij[k,0::3]=dFx
            Kij[k,1::3]=dFy
            Kij[k,2::3]=dFz

    Kij/=(2*HSTEP)

    KijFile=filePrefix+'.Kij.npy'
    with open(KijFile,'wb') as f:
        print "Writing to %s." % (KijFile)
        np.save(f,Kij)

    with open(filePrefix+'.ids_L.npy','wb') as f:
        np.save(f,ids_L)

    with open(filePrefix+'.ids_R.npy','wb') as f:
        np.save(f,ids_R)
            

    # return Kij, ids_L, ids_R

if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser()
    parser.add_argument("filePrefix",help="The used file prefix.")
    parser.add_argument("KijFilePrefix",nargs='?',help="The prefix for the used force constant file, script searches for *.Kij.npy where *=argument (optional, default=fileprefix)")
    args=parser.parse_args()
    filePrefix=args.filePrefix
    if args.KijFilePrefix is None:
        KijFilePrefix=filePrefix
    else:
        KijFilePrefix=args.KijFilePrefix
    import os
    if not os.path.isfile(KijFilePrefix+'.Kij.npy'):    
        print "Reading Kij values from the Fij.dat file."
        calcFC_fromFile(filePrefix+'.Fij.dat')
    
    #import sys
    #sys.exit()
    main(fileprefix=filePrefix,KijFilePrefix=KijFilePrefix)
