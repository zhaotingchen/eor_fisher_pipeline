import os
import sys
import numpy as np
import shutil
import h5py
from config_template import *

pars_list = ['HII_EFF_FACTOR','ION_Tvir_MIN']

for rng in rng_list:
    file_name = 'rng_'+str(rng)+'/run_rest.sh'
    arglist = ()
    save_file = 'rng_'+str(rng)+'/'+save_file_f
    file = '_f'
    with h5py.File(save_file,'r') as f:
        for i in range(len(step_size)):
            for j in pars_list:
                key ='fisher'+file+'_'+str(i)+'_'+j
                if key not in list(f.keys()):
                    arglist += ('python run21cmfast.py faint '+str(i)+' '+j,)
                else:
                    if 'tb_lc' not in list(f[key].keys()):
                        arglist += ('python run21cmfast.py faint '+str(i)+' '+j,)
        key = 'fisher'+file+'_fid'
        if key not in list(f.keys()):
            arglist += ('python run21cmfast.py faint '+str(i+1),)
        else:
            if 'tb_lc' not in list(f[key].keys()):
                arglist += ('python run21cmfast.py faint '+str(i+1),)
    
    save_file = 'rng_'+str(rng)+'/'+save_file_b
    file = '_b'
    with h5py.File(save_file,'r') as f:
        for i in range(len(step_size)):
            for j in pars_list:
                key ='fisher'+file+'_'+str(i)+'_'+j
                if key not in list(f.keys()):
                    arglist += ('python run21cmfast.py bright '+str(i)+' '+j,)
                else:
                    if 'tb_lc' not in list(f[key].keys()):
                        arglist += ('python run21cmfast.py bright '+str(i)+' '+j,)
        key = 'fisher'+file+'_fid'
        if key not in list(f.keys()):
            arglist += ('python run21cmfast.py bright '+str(i+1)+' faint',)
        else:
            if 'tb_lc' not in list(f[key].keys()):
                arglist += ('python run21cmfast.py bright '+str(i+1)+' filler',)
    
    with open('run_pipeline.sh','r') as f:
        lines = [line.rstrip() for line in f]
    
    tidy_syscall = []
    for line in lines:
        if '--array' in line:
            continue
        if '--job-name' in line:
            tidy_syscall += ['#SBATCH --job-name=tidy_sim_eor_'+str(rng)]
            continue
        if 'cd ' in line:
            break
        tidy_syscall += [line,]
    
    for arg in arglist:
        tidy_syscall += [arg,]
    tidy_syscall += ['python signalcov.py faint','python signalcov.py bright']
    tidy_syscall += ['echo "****ELAPSED "$SECONDS"s tidy_sim_eor"']
    with open(file_name,'w+') as f:
        for line in tidy_syscall:
            f.write(line+'\n')

rng_start = rng_list[0]
rng_finish = rng_list[-1]
filetxt = '''
for i in {{{rng_start}..{rng_finish}..1}}
do
cd rng_$i
sbatch run_rest.sh
cd ../
done
'''.format(**locals())
with open('tidy.sh','w+') as f:
    f.write(filetxt)