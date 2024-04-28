import scipy
from scipy.signal import blackmanharris
import torch
import py21cmfast as p21c
from py21cmfast import plotting
import os
import h5py
import time
from multiprocessing import cpu_count
from astropy.cosmology import Planck18
import numpy as np
from hiimtool.basic_util import Specs,p2dim,f_21,get_mask_renorm_simple,get_corr_mat
from hiimtool.hiimimage import get_power_3d,grid_im
from hiimtool.util import bin_3d_to_1d_lowmem,bin_3d_to_cy_lowmem,p2dim
from scipy.special import legendre
from torch.fft import fft,fftn,fftfreq
from functools import cached_property
import matplotlib.colors as colors
from hiimtool.power_util import MultiPower
from scipy.special import legendre
from astropy.io import fits
from astropy.wcs import WCS
from config import *
import sys

if torch.cuda.is_available():
    device='cuda'
else:
    device='cpu'

def write_into_hdf5(file,data,grp,key_name):
    with h5py.File(file,"a") as f:
        if key_name in list(f[grp].keys()):
            del f[grp][key_name]
        f[grp][key_name] = data

run_setting = sys.argv[-1] #faint or bright
if run_setting == 'faint':
    file = '_f'
elif run_setting == 'bright':
    file = '_b'
else:
    raise ValueError('run_setting has to be faint or bright')
save_file = locals()['save_file'+file]
box_len = user_pars['BOX_LEN'] #in Mpc
hii_dim = user_pars['HII_DIM'] #number of cells

freqarr = np.diff(freqarr)/2+freqarr[:-1]
num_ch = len(freqarr)
sp = Specs(cosmo=cosmo,
           freq_start_hz = freq_range[0],
           num_channels = num_ch,
           deltav_ch = np.diff(freqarr)[0],
           FWHM_ref=1, FWHM_freq_ref=1
)

comov_len = (cosmo.comoving_distance(sp.z_arr()).max()-cosmo.comoving_distance(sp.z_arr()).min()).to('Mpc').value

pos_lc = np.zeros((hii_dim,hii_dim,num_ch,3))
pos_lc[:,:,:,0] = (np.linspace(0,hii_dim,hii_dim)*box_len/(hii_dim))[:,None,None]
pos_lc[:,:,:,1] = (np.linspace(0,hii_dim,hii_dim)*box_len/(hii_dim))[None,:,None]
pos_lc[:,:,:,2] = (np.linspace(0,num_ch,num_ch)*(comov_len)/num_ch)[None,None,:]

def get_signal_cov(grp_name):
    with h5py.File(save_file,"r") as f:
        Tb_lc_smooth = np.array(f[grp_name]['tb_lc'])
        astro_pars = dict(f[grp_name].attrs)
    mp = MultiPower(
        pos = pos_lc,
        len_side = [box_len,box_len,pos_lc[:,:,:,2].max()],
        N_side = [hii_dim,hii_dim,num_ch],
        cosmo=cosmo,
        field=Tb_lc_smooth,
    )
    kpara_lc = mp.k_para
    kz_edges = np.sort(np.unique(np.abs(kpara_lc)))
    kz_edges = np.append(kz_edges[:-1]-np.diff(kz_edges)/2,kz_edges[-2:]+np.diff(kz_edges)[-1]/2)
    mp.kperpedges = kz_edges[::2]
    
    k_lc = mp.k_mode
    kperp = mp.k_perp
    kpara = mp.k_para
    kz = kpara[None,None,:]/k_lc
    pl_2 = legendre(2)(kz)
    pl_4 = legendre(4)(kz)
    power_3d_0 = mp.power_3d_0
    kx_sample = np.load(kx_sample_file)
    ky_sample = np.load(ky_sample_file)
    count_i =  np.load(count_i_file)
    kx_vec = mp.kxvec
    ky_vec = mp.kyvec
    kx_vec = np.fft.fftshift(kx_vec)
    ky_vec = np.fft.fftshift(ky_vec)
    kx_vec = np.append(kx_vec-np.diff(kx_vec)[0]/2,kx_vec[-1]+np.diff(kx_vec)[0]/2)
    ky_vec = np.append(ky_vec-np.diff(ky_vec)[0]/2,ky_vec[-1]+np.diff(ky_vec)[0]/2)
    
    count_grid,_,_ = np.histogram2d(kx_sample,ky_sample,bins=[kx_vec,ky_vec])
    assert count_grid.max()==1 # Grid resolution is high enough
    power_3d_0_grid = power_3d_0.copy()
    power_3d_0_grid = np.fft.fftshift(power_3d_0_grid,axes=(0,1))
    power_sample = power_3d_0_grid[count_grid>0]
    write_into_hdf5(save_file,power_sample,grp_name,'power_sample')
    #with h5py.File(save_file,"a") as f:
    #    f[grp_name]['power_sample'] = power_sample
    kperp_sp = np.sqrt(kx_sample**2+ky_sample**2)
    ksp_mode = np.sqrt(kperp_sp[:,None]**2+kpara[None,:]**2)
    kz_sp = kpara[None,:]/ksp_mode
    pl2_sp = legendre(2)(kz_sp)
    pl4_sp = legendre(4)(kz_sp)
    weight_grid,_,_ = np.histogram2d(kx_sample,ky_sample,bins=[kx_vec,ky_vec],weights=count_i)
    weight_grid = np.broadcast_to(weight_grid[weight_grid>0][:,None],ksp_mode.shape)
    ptbcy_0 = bin_3d_to_cy_lowmem(power_3d_0.T,kperp.T,kperpedges)
    ptbcy_0 = np.fft.fftshift(ptbcy_0,axes=0)
    write_into_hdf5(save_file,kperpedges,grp_name,'kperpedges')
    write_into_hdf5(save_file,ptbcy_0,grp_name,'ptbcy_0')
    #with h5py.File(save_file,"a") as f:
    #    f[grp_name]['kperpedges'] = kperpedges
    #    f[grp_name]['ptbcy_0'] = ptbcy_0
    i_indx = 0
    j_indx = 0
    nx,ny,nz = Tb_lc_smooth.shape
    nx_sub = nx//nsub_axis
    ny_sub = ny//nsub_axis
    
    mask = np.ones_like(Tb_lc_smooth)
    mask[i_indx*nx_sub:(i_indx+1)*nx_sub,
    j_indx*ny_sub:(j_indx+1)*ny_sub] = 0.0
    
    mask_renorm = get_mask_renorm_simple(mask)
    
    k1dcen = (k1dedges[1:]+k1dedges[:-1])/2
    
    ptb1d_0_sub = np.zeros((nsub_axis**2,len(k1dedges)-1))
    ptb1d_2_sub = np.zeros((nsub_axis**2,len(k1dedges)-1))
    ptb1d_4_sub = np.zeros((nsub_axis**2,len(k1dedges)-1))
    
    i = 0
    for i_indx in range(nsub_axis):
        for j_indx in range(nsub_axis):
            #print(i)
            mask = np.ones_like(Tb_lc_smooth)
            mask[i_indx*nx_sub:(i_indx+1)*nx_sub,
            j_indx*ny_sub:(j_indx+1)*ny_sub] = 0.0
            mp_sub = MultiPower(
                     pos = pos_lc,
                     len_side = [box_len,box_len,pos_lc[:,:,:,2].max()],
                     N_side = [hii_dim,hii_dim,num_ch],
                     cosmo=cosmo,
                     field=Tb_lc_smooth*mask,
                     window=blackmanharris(num_ch),
                     device=device,
                     )
            #sub_renorm = get_mask_renorm_simple(mask)
            power_3d_0_sub = mp_sub.power_3d_0*mask_renorm
            power_3d_0_grid = power_3d_0_sub.copy()
            power_3d_0_grid = np.fft.fftshift(power_3d_0_grid,axes=(0,1))
            power_sample = power_3d_0_grid[count_grid>0]
            ptb1d_0_sub[i] = bin_3d_to_1d_lowmem(power_sample,ksp_mode,k1dedges,weights=weight_grid)
            ptb1d_2_sub[i] = bin_3d_to_1d_lowmem(power_sample*pl2_sp*5,ksp_mode,k1dedges,weights=weight_grid)
            ptb1d_4_sub[i] = bin_3d_to_1d_lowmem(power_sample*pl4_sp*9,ksp_mode,k1dedges,weights=weight_grid)
            i += 1
    ptbsub_arr = np.array([ptb1d_0_sub,ptb1d_2_sub,ptb1d_4_sub])
    diff_arr = ptbsub_arr - ptbsub_arr.mean(axis=1)[:,None,:]
    diff_arr = np.transpose(diff_arr,axes=(0,-1,1)).reshape((-1,nsub_axis**2))
    cov_ps = (diff_arr[:,None,:]*diff_arr[None,:,:]).mean(axis=-1)*(nsub_axis**2-1)/(nsub_axis**2)
    corr_ps = cov_ps/np.sqrt(np.diagonal(cov_ps))[:,None]/np.sqrt(np.diagonal(cov_ps))[None,:]
    write_into_hdf5(save_file,k1dedges,grp_name,'k1dedges')
    write_into_hdf5(save_file,cov_ps,grp_name,'cov_ps')
    write_into_hdf5(save_file,ptbsub_arr,grp_name,'ptbsub_arr')
    #with h5py.File(save_file,"a") as f:
    #    f[grp_name]['k1dedges'] = k1dedges
    #    f[grp_name]['cov_ps'] = cov_ps
    #    f[grp_name]['ptbsub_arr'] = ptbsub_arr
    mu_list = np.linspace(mu_beam,1.0,n_wedge+1)
    mu_list = np.append(0.0,mu_list)
    mu_sample = np.abs(kz_sp)
    i_indx = 0
    j_indx = 0
    nx,ny,nz = Tb_lc_smooth.shape
    nx_sub = nx//nsub_axis
    ny_sub = ny//nsub_axis
    
    mask = np.ones_like(Tb_lc_smooth)
    mask[i_indx*nx_sub:(i_indx+1)*nx_sub,
    j_indx*ny_sub:(j_indx+1)*ny_sub] = 0.0
    
    mask_renorm = get_mask_renorm_simple(mask)
    
    
    ptb1d_0_sub_ad = np.zeros((nsub_axis**2,len(k1dedges)-1))
    ptb1d_2_sub_ad = np.zeros((nsub_axis**2,len(k1dedges)-1))
    ptb1d_4_sub_ad = np.zeros((nsub_axis**2,len(k1dedges)-1))
    
    
    i = 0
    for i_indx in range(nsub_axis):
        for j_indx in range(nsub_axis):
            #print(i)
            mask = np.ones_like(Tb_lc_smooth)
            mask[i_indx*nx_sub:(i_indx+1)*nx_sub,
            j_indx*ny_sub:(j_indx+1)*ny_sub] = 0.0
            mp_sub = MultiPower(
                     pos = pos_lc,
                     len_side = [box_len,box_len,pos_lc[:,:,:,2].max()],
                     N_side = [hii_dim,hii_dim,num_ch],
                     cosmo=cosmo,
                     field=Tb_lc_smooth*mask,
                     window=blackmanharris(num_ch),
                     device=device,
                     )
            #sub_renorm = get_mask_renorm_simple(mask)
            power_3d_0_sub = mp_sub.power_3d_0*mask_renorm
            power_3d_0_grid = power_3d_0_sub.copy()
            power_3d_0_grid = np.fft.fftshift(power_3d_0_grid,axes=(0,1))
            power_sample = power_3d_0_grid[count_grid>0]
            #ptb1d_0_sub[i] = bin_3d_to_1d_lowmem(power_sample,ksp_mode,k1dedges)
            #ptb1d_2_sub[i] = bin_3d_to_1d_lowmem(power_sample*pl2_sp*5,ksp_mode,k1dedges)
            #ptb1d_4_sub[i] = bin_3d_to_1d_lowmem(power_sample*pl4_sp*9,ksp_mode,k1dedges)
            sel_indx = (np.abs(mu_sample)>=mu_beam)
            ptb1d_0_sub_ad[i] = bin_3d_to_1d_lowmem(power_sample,ksp_mode,k1dedges,weights=weight_grid*sel_indx)
            ptb1d_2_sub_ad[i] = bin_3d_to_1d_lowmem(power_sample*pl2_sp*5,ksp_mode,k1dedges,weights=weight_grid*sel_indx)
            ptb1d_4_sub_ad[i] = bin_3d_to_1d_lowmem(power_sample*pl4_sp*9,ksp_mode,k1dedges,weights=weight_grid*sel_indx)
            i += 1

    ptbsub_arr_ad = np.array([ptb1d_0_sub_ad,ptb1d_2_sub_ad,ptb1d_4_sub_ad])
    diff_arr = ptbsub_arr_ad - ptbsub_arr_ad.mean(axis=1)[:,None,:]
    diff_arr = np.transpose(diff_arr,axes=(0,-1,1)).reshape((-1,nsub_axis**2))
    cov_ps_ad = (diff_arr[:,None,:]*diff_arr[None,:,:]).mean(axis=-1)*(nsub_axis**2-1)/(nsub_axis**2)
    corr_ps_ad = cov_ps_ad/np.sqrt(np.diagonal(cov_ps_ad))[:,None]/np.sqrt(np.diagonal(cov_ps_ad))[None,:]
    write_into_hdf5(save_file,cov_ps_ad,grp_name,'cov_ps_ad')
    write_into_hdf5(save_file,ptbsub_arr_ad,grp_name,'ptbsub_arr_ad')
    
    #with h5py.File(save_file,"a") as f:
    #    data = f[grp_name]
    #    data['cov_ps_ad'] = cov_ps_ad
    #    data['ptbsub_arr_ad'] = ptbsub_arr_ad

    i_indx = 0
    j_indx = 0
    nx,ny,nz = Tb_lc_smooth.shape
    nx_sub = nx//nsub_axis
    ny_sub = ny//nsub_axis
    
    mask = np.ones_like(Tb_lc_smooth)
    mask[i_indx*nx_sub:(i_indx+1)*nx_sub,
    j_indx*ny_sub:(j_indx+1)*ny_sub] = 0.0
    
    mask_renorm = get_mask_renorm_simple(mask)
    
    
    ptb1d_0_sub_w = np.zeros((len(mu_list)-1,nsub_axis**2,len(k1dedges)-1))
    ptb1d_2_sub_w = np.zeros((len(mu_list)-1,nsub_axis**2,len(k1dedges)-1))
    ptb1d_4_sub_w = np.zeros((len(mu_list)-1,nsub_axis**2,len(k1dedges)-1))
    
    
    i = 0
    for i_indx in range(nsub_axis):
        for j_indx in range(nsub_axis):
            #print(i)
            mask = np.ones_like(Tb_lc_smooth)
            mask[i_indx*nx_sub:(i_indx+1)*nx_sub,
            j_indx*ny_sub:(j_indx+1)*ny_sub] = 0.0
            mp_sub = MultiPower(
                     pos = pos_lc,
                     len_side = [box_len,box_len,pos_lc[:,:,:,2].max()],
                     N_side = [hii_dim,hii_dim,num_ch],
                     cosmo=cosmo,
                     field=Tb_lc_smooth*mask,
                     window=blackmanharris(num_ch),
                     device=device,
                     )
            #sub_renorm = get_mask_renorm_simple(mask)
            power_3d_0_sub = mp_sub.power_3d_0*mask_renorm
            power_3d_0_grid = power_3d_0_sub.copy()
            power_3d_0_grid = np.fft.fftshift(power_3d_0_grid,axes=(0,1))
            power_sample = power_3d_0_grid[count_grid>0]
            for wedge_indx in range(len(mu_list)-1):
                sel_indx = (mu_sample>=mu_list[wedge_indx])*(mu_sample<mu_list[wedge_indx+1])
                ptb1d_0_sub_w[wedge_indx,i] = bin_3d_to_1d_lowmem(power_sample,ksp_mode,k1dedges,weights=weight_grid*sel_indx)
                ptb1d_2_sub_w[wedge_indx,i] = bin_3d_to_1d_lowmem(power_sample*pl2_sp*5,ksp_mode,k1dedges,weights=weight_grid*sel_indx)
                ptb1d_4_sub_w[wedge_indx,i] = bin_3d_to_1d_lowmem(power_sample*pl4_sp*9,ksp_mode,k1dedges,weights=weight_grid*sel_indx)
            i += 1
    cov_ps_wedge = []
    corr_ps_wedge = []
    for wedge_indx in range(len(mu_list)-1):
        ptbsub_arr_w = np.array([ptb1d_0_sub_w[wedge_indx],ptb1d_2_sub_w[wedge_indx],ptb1d_4_sub_w[wedge_indx]])
        diff_arr = ptbsub_arr_w - ptbsub_arr_w.mean(axis=1)[:,None,:]
        diff_arr = np.transpose(diff_arr,axes=(0,-1,1)).reshape((-1,nsub_axis**2))
        cov_ps_w = (diff_arr[:,None,:]*diff_arr[None,:,:]).mean(axis=-1)*(nsub_axis**2-1)/(nsub_axis**2)
        corr_ps_w = cov_ps_w/np.sqrt(np.diagonal(cov_ps_w))[:,None]/np.sqrt(np.diagonal(cov_ps_w))[None,:]
        cov_ps_wedge += [cov_ps_w,]
        corr_ps_wedge += [corr_ps_w,]
    
    corr_ps_wedge = np.array(corr_ps_wedge)
    cov_ps_wedge = np.array(cov_ps_wedge)
    write_into_hdf5(save_file,np.array([ptb1d_0_sub_w,ptb1d_2_sub_w,ptb1d_4_sub_w]),grp_name,'ptbsub_arr_wedge')
    write_into_hdf5(save_file,cov_ps_wedge,grp_name,'cov_ps_wedge')
    
    #with h5py.File(save_file,"a") as f:
    #    f[grp_name]['ptbsub_arr_wedge'] = np.array([ptb1d_0_sub_w,ptb1d_2_sub_w,ptb1d_4_sub_w])
    #    f[grp_name]['cov_ps_wedge'] = cov_ps_wedge
    return 1

with h5py.File(save_file,"r") as f:
    keys_list = list(f.keys())

for grp_n in keys_list:
    if 'fisher' in grp_n:
        with h5py.File(save_file,"r") as f:
            data_list = list(f[grp_n].keys())
        if 'cov_ps_wedge' not in data_list:
            get_signal_cov(grp_n)
    

