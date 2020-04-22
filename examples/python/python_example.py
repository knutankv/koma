#!/usr/bin/env python
# coding: utf-8

# # Python example

# In[2]:


import koma, koma.oma

import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import resample, detrend

import pandas as pd


# ## Load data and define input
# 

# In[13]:


# Import and specify properties of input data
df = pd.read_csv('response_data.csv', sep=',',header=None)
data = df.values      # data including noise
levels = 3
fs = 35
t = np.arange(0,1/fs*(data.shape[0]),1/fs)

# Cov-SSI settings
i = 20
s = 4
orders = np.arange(2, 50+2, 2)
stabcrit = {'freq':0.2, 'damping': 0.2, 'mac': 0.3}

# Noise specification
noise_factor = 1.0


# ## Add artificial noise and plot response

# In[4]:


noise = np.std(data) * noise_factor
data_noised = data + noise*np.random.randn(data.shape[0], data.shape[1])

fig, ax = plt.subplots(nrows=3, ncols=1, num=1)
ax[0].plot(t, data_noised[:,0])
ax[0].plot(t, data[:,0], color='IndianRed', alpha=1)
ax[1].plot(t, data_noised[:,3])
ax[1].plot(t, data[:,3], color='IndianRed', alpha=1)
ax[2].plot(t, data_noised[:,6])
ax[2].plot(t, data[:,6], color='IndianRed', alpha=1)

__ = [a.set_xticks([]) for a in ax[0:2]]
__ = ax[2].set_xlabel('Time [s]')

__ = ax[0].set_ylabel('Level 3')
__ = ax[1].set_ylabel('Level 2')
__ = ax[2].set_ylabel('Level 1')


# ## Preprocess

# In[5]:


fs_rs = 3  # resample frequency
data_rs = resample(detrend(data), int(np.ceil(data.shape[0]*fs_rs/fs)), axis=0) # detrend and resample
fs_rs = data_rs.shape[0]/data.shape[0]*fs


# ## Cov-SSI call

# In[6]:


lambd, phi = koma.oma.covssi(data_rs, fs_rs, i, orders)


# ## Postprocessing and visualization

# In[14]:


# Establishing stable poles
s = 2
lambd_stab, phi_stab, orders_stab, ix_stab = koma.oma.find_stable_poles(lambd, phi, orders, s, stabcrit=stabcrit)


# In[15]:


# 
print(lambd_stab)

