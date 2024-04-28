from astropy.cosmology import Planck18
import numpy as np

cosmo=Planck18

cache_dir = '_cache'
temp_dir = 'temp'

save_file_f = "sim_result_faint.hdf5"
save_file_b = "sim_result_bright.hdf5"

#fiducial
astro_pars_f = {
    'HII_EFF_FACTOR': 65.0,
    'R_BUBBLE_MAX ': 50.0,
    'ION_Tvir_MIN': 4.70,
    'L_X': 40.0,
    'NU_X_THRESH':500.0,
    'X_RAY_SPEC_INDEX':1.0,
}

astro_pars_b = {
    'HII_EFF_FACTOR': 150.0,
    'R_BUBBLE_MAX ': 50.0,
    'ION_Tvir_MIN': 5.10,
    'L_X': 40.0,
    'NU_X_THRESH':500.0,
    'X_RAY_SPEC_INDEX':1.0,
}

step_size = [0.001,0.005,0.01,0.02,0.05,0.1]

quantities = ("brightness_temp", 'density', 'xH_box','velocity')
flag_op = {"USE_TS_FLUCT":True,"SUBCELL_RSD":True,'INHOMO_RECO':True}
z_min = 7.0
user_pars = {"HII_DIM":400, "BOX_LEN": 800,"N_THREADS":None,
            "MINIMIZE_MEMORY":False,}

freqarr = np.linspace(145*1e6,175*1e6,151)

freq_range = (145*1e6,175.1*1e6) # in Hz

kx_sample_file = '/home/zchen/inter_sim/vis_power/kx_i.npy'
ky_sample_file = '/home/zchen/inter_sim/vis_power/ky_i.npy'
count_i_file =  '/home/zchen/inter_sim/vis_power/count_i.npy'

#for cylindrical power
kperpedges = np.logspace(-2,0,21)

# number of slices per axis along transverse plane for jackknife
nsub_axis = 5

#for 1D power
k1dedges = np.logspace(np.log10(0.05),0,21)

mu_beam = 0.26755280609040966

n_wedge = 3

#rng
