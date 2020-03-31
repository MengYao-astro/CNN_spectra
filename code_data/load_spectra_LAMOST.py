#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.interpolate import Akima1DInterpolator


# In[2]:


# Read file list text
list_LAMOST_txt = open('../../data/LAMOST/dr5_SNR100/list.txt').readlines()
# complete fits' names with address
list_LAMOST = ['../../data/LAMOST/dr5_SNR100/' + i.replace('\n','') for i in list_LAMOST_txt[1:]]
print('dataset size = ',len(list_LAMOST))


# ## Check out the wavelength range

# In[3]:


wv_min = np.zeros(len(list_LAMOST))
wv_max = np.zeros(len(list_LAMOST))
for i in range(len(list_LAMOST)):
    hdul = fits.open(list_LAMOST[i])
    wavelength = hdul[0].data[2]
    wv_min[i] = min(wavelength)
    wv_max[i] = max(wavelength)
inter_start = int(max(wv_min))+1
inter_end = int(min(wv_max))
print('interpolation should start at wv = ',inter_start)
print('interpolation should stop at wv = ',inter_end)


# ### Run a test on one spectrum

# In[22]:


print('\ntest:')
#test_on_num = 31505
test_on_num = 20
print(list_LAMOST[test_on_num])
# read fits
hdul_test = fits.open(list_LAMOST[test_on_num])
# print subclass
subclass_test = hdul_test[0].header['SUBCLASS']
print('subclass = ',subclass_test)

"""
try:
    teff_test = hdul_test[0].header['teff']
    logg_test = hdul_test[0].header['logg']
    print('T_eff = ', teff_test, ' K')
    print("log_g = ", logg_test)
except:
    print('Teff and logg are Not found')
"""

#read flux and wavelength
flux_test = hdul_test[0].data[0]
wl_test = hdul_test[0].data[2]
# count data points
print('the number of data points = ', len(wl_test))
fig0, axs0 = plt.subplots(3, 1, dpi = 150, figsize=(6,8))
#fig0.suptitle('{}'.format(subclass_test))
axs0[0].plot(wl_test,flux_test,c='b',linewidth=0.7)
axs0[0].set_ylabel('flux')
axs0[0].set_title('Data')



cs_test = CubicSpline(wl_test, flux_test)
wv_inter_test = np.linspace(inter_start,9000,4096)
flux_inter_test_cs = cs_test(wv_inter_test)
axs0[1].plot(wv_inter_test,flux_inter_test_cs,linewidth=0.7)
axs0[1].set_ylabel('flux')
axs0[1].set_title('cubic spline')


akima_test = Akima1DInterpolator(wl_test, flux_test)
flux_inter_test_akima = akima_test(wv_inter_test)
axs0[2].plot(wv_inter_test,flux_inter_test_akima,linewidth=0.7)
axs0[2].set_ylabel('flux')
axs0[2].set_xlabel('wavelength')
axs0[2].set_title('Akima')

fig0.tight_layout()
fig0.show()

print('relative diff average = ',sum(abs(flux_inter_test_cs - flux_inter_test_akima)/flux_inter_test_akima)/4096)


# In[5]:


hdul_test[0].data[2]


# ## Load LAMOST data

# In[6]:


def load_spec(list_spec,inter_start,inter_end):
    '''
    Define a function that can generate arrays gives flux and wavelength of spectra, as well as their subclass.
    list_spec : a 'list'(data type) of spectra's names, shoulds contain complete folder address.
                For example, 'data1/yao/2nd_research/data/LAMOST/spec***********.fits' for each row.
    '''
    num_spec = len(list_spec)   # length of the list, shows the number of spectra
    n_datapoints = 4096
    fits_obsid = np.zeros((num_spec), dtype=int)  
    spectra_flux = np.zeros((num_spec,n_datapoints),dtype = np.float32) # prepare an array to store flux 
    #spectra_wl = np.zeros((num_spec,n_datapoints))   # prepare an array to store wavelength
    label = np.zeros((num_spec), dtype=int)             # prepare an array to store class label
    interpolate_wv = np.linspace(inter_start,inter_end,n_datapoints)
    for i in range (num_spec):   # for every fits
        hdul = fits.open(list_spec[i])  # open it 
        fits_obsid[i] = hdul[0].header["OBSID"]
        flux = hdul[0].data[0] # get the flux 
        wavelength = hdul[0].data[2]   # get the wavelength
        akima = Akima1DInterpolator(wavelength,flux)  # Using akima interpolation
        spectra_flux[i] = akima(interpolate_wv) # write flux into spectta
        """
        subclass = hdul[0].header['SUBCLASS']  # read the subclass
        if subclass[0] == 'O':     # if type O, label = 0
            label[i] = 0
        elif subclass[0] == 'B':   # if      B,         1
            label[i] = 1
        elif subclass[0] == 'A':   # if      A,         2
            label[i] = 2
        elif subclass[0] == 'F':   # if      F,         3
            label[i] = 3
        elif subclass[0] == 'G':   # if      G,         4
            label[i] = 4
        elif subclass[0] == 'K':   # if      K,         5
            label[i] = 5
        elif subclass[0] == 'M':   # if      M,         6
            label[i] = 6
        else:                      # else, save -1 and will delete those fits data afterwards.
            label[i] = -1
        """
            
    return interpolate_wv, spectra_flux, fits_obsid        


# In[22]:


wv_inter, spectra_flux, fits_obsid= load_spec(list_LAMOST,inter_start, inter_end)
stellar_spectra = spectra_flux

np.save('../../data/npy/LAMOST_dr5_SNR100_obsid', fits_obsid)
np.save('../../data/npy/LAMOST_dr5_SNR100_spectra',stellar_spectra)
#np.save('../../data/npy/LAMOST_label_SNR100_ab',stellar_label)



