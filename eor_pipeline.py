from config_template import rng_list,step_size
import shutil
import os
import sys

for rng in rng_list:
    if not os.path.exists('rng_'+str(rng)):
        os.mkdir('rng_'+str(rng))
    
for rng in rng_list:
    shutil.copyfile('config_template.py','rng_'+str(rng)+'/config.py')
    shutil.copyfile('initialize.py','rng_'+str(rng)+'/initialize.py')
    shutil.copyfile('run21cmfast.py','rng_'+str(rng)+'/run21cmfast.py')
    shutil.copyfile('signalcov.py','rng_'+str(rng)+'/signalcov.py')
    

for rng in rng_list:
    with open('rng_'+str(rng)+'/config.py','a') as f:
        f.write('rng = '+str(rng))

with open('template.sh') as file:
    lines = [line.rstrip() for line in file]

flag = True
for i,line in enumerate(lines):
    if not flag and '#SBATCH' not in line:
        break
    if flag and '#SBATCH' in line:
        flag=False

lines[i] = '#SBATCH --array='+str(rng_list[0])+'-'+str(rng_list[-1])

for i,line in enumerate(lines):
    if 'do' == line:
        break
lines[i-1] = 'for i in {0..'+str(len(step_size)-1)+'..1}'

flag = True
for i,line in enumerate(lines):
    if not flag and 'do' == line:
        break
    if flag and 'do' == line:
        flag=False
lines[i-1] = 'for i in {0..'+str(len(step_size)-1)+'..1}'

for i,line in enumerate(lines):
    if 'done' == line:
        break
lines[i+1] = 'python run21cmfast.py faint '+str(len(step_size))+' HII_EFF_FACTOR'

flag = True
for i,line in enumerate(lines):
    if not flag and 'done' == line:
        break
    if flag and 'done' == line:
        flag=False
lines[i+1] = 'python run21cmfast.py bright '+str(len(step_size))+' HII_EFF_FACTOR'

with open('run.sh', 'w') as f:
    for line in lines:
        f.write(f"{line}\n")