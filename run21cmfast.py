import py21cmfast as p21c
from py21cmfast import plotting
import os
import h5py
import time
from multiprocessing import cpu_count
from astropy.cosmology import Planck18
import numpy as np
from hiimtool.basic_util import Specs,f_21
import scipy
import sys
from config import *

run_setting = sys.argv[-3] #faint or bright
i = int(sys.argv[-2]) # which step size
j = sys.argv[-1] # which parameter
if run_setting == 'faint':
    file = '_f'
elif run_setting == 'bright':
    file = '_b'
else:
    raise ValueError('run_setting has to be faint or bright')

save_file = locals()['save_file'+file]

max_num_thread = cpu_count()
#print(max_num_thread)

#fiducial
astro_pars = locals()['astro_pars'+file]

#step_size

if i>= len(locals()['step_size']):
    #fiducial
    grp_name = 'fisher'+file+'_fid'
else:
    step_size = locals()['step_size'][i]
    astro_pars[j] = astro_pars[j]*(1+step_size)
    grp_name = 'fisher'+file+'_'+str(i)+'_'+j

#if i == 0:
#    astro_pars['HII_EFF_FACTOR'] = astro_pars['HII_EFF_FACTOR']*(1+step_size)
#elif i == 1:
#    astro_pars['HII_EFF_FACTOR'] = astro_pars['HII_EFF_FACTOR']*(1-step_size)
#elif i == 2:
#    astro_pars['ION_Tvir_MIN'] = astro_pars['ION_Tvir_MIN']*(1+step_size)
#elif i == 3:
#    astro_pars['ION_Tvir_MIN'] = astro_pars['ION_Tvir_MIN']*(1-step_size)

user_pars["N_THREADS"]= max_num_thread

st_time = time.time()
lightcone = p21c.run_lightcone(
    redshift = z_min,
    flag_options = flag_op,
    user_params = user_pars,
    lightcone_quantities=quantities,
    global_quantities=quantities,
    astro_params = astro_pars,
    random_seed = rng,
    direc=cache_dir,
    write=True,
)
print(time.time()-st_time)

fname = lightcone.save(fname='lightcone.h5', direc=temp_dir,)

lightcone.gather(
    fname=fname,
    direc=cache_dir,
    kinds=("brightness_temp", "init","spin_temp","ionized_box","perturb_field"),
    clean=True,
);

lc = p21c.LightCone.read(fname) # read lightcone
lc_z = lc.lightcone_redshifts


box_len = user_pars['BOX_LEN'] #in Mpc
hii_dim = user_pars['HII_DIM'] #number of cells

z_range = (f_21/freq_range[-1]-1,f_21/freq_range[0]-1)
sel_node = ((lc_z>z_range[0])*(lc_z<z_range[1]))
z_range = (lc_z[sel_node].min(),lc_z[sel_node].max())
freq_range = (f_21/(1+lc_z[sel_node].max()),f_21/(1+lc_z[sel_node].min()))
num_node = sel_node.sum()
delta_dc = box_len/hii_dim
comov_len = delta_dc*(num_node-1)

freqarr = np.diff(freqarr)/2+freqarr[:-1]
num_ch = len(freqarr)
sp = Specs(cosmo=cosmo,
           freq_start_hz = freq_range[0],
           num_channels = num_ch,
           deltav_ch = np.diff(freqarr)[0],
           FWHM_ref=1, FWHM_freq_ref=1
)

# read the temperature lightcone
simfile = h5py.File(fname,"r")
Tb_lc = np.array(simfile['lightcones']['brightness_temp'])[:,:,sel_node]

# apply a smoothing function to match the frequency resolution setting
Tb_lc_smooth = np.zeros((Tb_lc.shape+(3,)))
Tb_lc_smooth = np.broadcast_to(Tb_lc[:,:,:,None],Tb_lc_smooth.shape)
Tb_lc_smooth = Tb_lc_smooth.reshape((Tb_lc.shape[:2]+(-1,)))
smooth_window = np.exp(-(np.linspace(-2,2,5)**2/2/2.5**2))
smooth_window /= smooth_window.sum()
Tb_lc_smooth = scipy.ndimage.convolve1d(Tb_lc_smooth, smooth_window, axis=-1,)

z_coor_super = np.linspace(0,comov_len,Tb_lc_smooth.shape[-1])
z_coor_sub = cosmo.comoving_distance(sp.z_arr()).value
z_coor_sub -= z_coor_sub.min()
z_corr_indx = np.argmin(np.abs(z_coor_sub[:,None]-z_coor_super[None,:]),axis=-1)

Tb_lc_smooth = Tb_lc_smooth[:,:,z_corr_indx]

f = h5py.File(save_file,'a')
f.create_group(grp_name)
grp = f[grp_name]
grp.attrs.update(astro_pars)
dset = grp.create_dataset(
    "tb_lc", data = Tb_lc_smooth)
f.close()

os.remove(fname)


