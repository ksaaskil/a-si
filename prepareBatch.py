# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015

def prepareBatch(filename,executable,inputFile,execParams='',
                 outputFile='exec_output.out',
                 time='00:30:00',nodes=1,ntasks=1,
                 partition='test',
                 out='slurm_output.out',
                 system='TAITO',
                 dependency=None,
                 jobname=None,
                 mem_per_cpu='1000M'):
    ''' 
    Prepare a batch file for SLURM. 

    Arguments: Self-explanatory.
    '''
    if system=='TAITO':
        ntasks_per_node=16
    elif system=='SISU':
        ntasks_per_node=24
    else:
        raise ValueError, "System must be 'TAITO' or 'SISU'."

    print "Preparing the batch file "+filename+'.'
    with open(filename,'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH -t '+time+'\n')
        f.write('#SBATCH -p '+partition+'\n')
        f.write('#SBATCH --nodes='+str(nodes)+'\n')
        f.write('#SBATCH -n '+str(ntasks)+'\n')
        f.write('#SBATCH -o '+out+'\n')
        f.write('#SBATCH --mem-per-cpu='+str(mem_per_cpu)+'\n')
        f.write('#SBATCH --ntasks-per-node='+str(ntasks_per_node)+'\n')
        
        if dependency is not None and not dependency.endswith(':'):
            f.write('#SBATCH --dependency='+str(dependency)+'\n')
        if jobname is not None:
            f.write('#SBATCH --job-name='+str(jobname)+'\n')
        if system=='TAITO':
            f.write('#SBATCH --constraint=snb\n')
            f.write('srun '+executable+' '+execParams+' < '+inputFile+' > ' + outputFile +'\n')
        elif system=='SISU':
            f.write('aprun -n '+str(ntasks)+' -m '+ str(mem_per_cpu)+' '+executable+' '+execParams+ ' < '+inputFile+' > ' + outputFile+'\n')         
            
    
