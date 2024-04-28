import h5py
from config import save_file_f,save_file_b,astro_pars_f,astro_pars_b,cache_dir,temp_dir,rng,step_size
import sys
import os

# create the save file
f = h5py.File(save_file_f, "w")
# store the parameters
f.create_group('fiducial_pars')
grp = f['fiducial_pars']
grp.attrs.update(astro_pars_f)
grp.attrs.update({'rng':rng})
f.close()

f = h5py.File(save_file_b, "w")
f.create_group('fiducial_pars')
grp = f['fiducial_pars']
grp.attrs.update(astro_pars_b)
grp.attrs.update({'rng':rng})
f.close()

if not os.path.isdir(cache_dir):
    os.mkdir(cache_dir)

if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)