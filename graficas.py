# -*- coding: utf-8 -*-
"""
Created on Sun Dec 30 13:46:28 2018

@author: malopez
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Reading data from files
datos_theta001 = pd.read_pickle('./alpha09_beta0_theta001_mu5.pkl')
datos_theta10 = pd.read_pickle('./alpha09_beta0_theta10_mu4.pkl')
datos_theta100 = pd.read_pickle('./alpha09_beta0_theta100_mu3.pkl')

antonio_001 = pd.read_csv('./datosAntonio/sim_alpha09_beta0_T5.dat', sep=' ', names=['time', 'temp'])
antonio_10 = pd.read_csv('./datosAntonio/sim_alpha09_beta0_T4.dat', sep=' ', names=['time', 'temp'])
antonio_100 = pd.read_csv('./datosAntonio/sim_alpha09_beta0_T3.dat', sep=' ', names=['time', 'temp'])

teoricos_001 = pd.read_csv('./datosTeoricos/alpha09_beta0_T5.dat', sep='\t', names=['time', 'temp'])
teoricos_10 = pd.read_csv('./datosTeoricos/alpha09_beta0_T4.dat', sep='\t', names=['time', 'temp'])
teoricos_100 = pd.read_csv('./datosTeoricos/alpha09_beta0_T3.dat', sep='\t', names=['time', 'temp'])

# From the last 2000 datapoints of each simulation, I calculate the mean stationary temp.
temp_stationary = pd.concat((datos_theta001['temp'][-2000:], 
                             datos_theta10['temp'][-2000:],
                             datos_theta100['temp'][-2000:]), axis=0).mean()
temp_tras_stationary = pd.concat((datos_theta001['temp_tras'][-2000:], 
                             datos_theta10['temp_tras'][-2000:],
                             datos_theta100['temp_tras'][-2000:]), axis=0).mean()

# =============================================================================
# Intento de corregir el escalamiento en el tiempo
# =============================================================================
# =============================================================================
# datos_theta001['time'] = datos_theta001['time']*temp_tras_stationary / (np.sqrt(np.pi) * temp_stationary)
# datos_theta10['time'] = datos_theta10['time']*temp_tras_stationary / (np.sqrt(np.pi) * temp_stationary)
# datos_theta100['time'] = datos_theta100['time']*temp_tras_stationary / (np.sqrt(np.pi) * temp_stationary)
# =============================================================================

# =============================================================================
# datos_theta001['time'] = datos_theta001['time']*np.sqrt(temp_tras_stationary / (np.pi * temp_stationary))
# datos_theta10['time'] = datos_theta10['time']*np.sqrt(temp_tras_stationary / (np.pi * temp_stationary))
# datos_theta100['time'] = datos_theta100['time']*np.sqrt(temp_tras_stationary / (np.pi * temp_stationary))
# =============================================================================

datos_theta001['time'] = datos_theta001['time']*np.sqrt(temp_tras_stationary / (np.pi))
datos_theta10['time'] = datos_theta10['time']*np.sqrt(temp_tras_stationary / (np.pi))
datos_theta100['time'] = datos_theta100['time']*np.sqrt(temp_tras_stationary / (np.pi))


# Now, we scale the temperature with respect to the stationary value
# so that the stationary value is 1
datos_theta001.temp /= temp_stationary
datos_theta10.temp /= temp_stationary
datos_theta100.temp /= temp_stationary


# Plot
plt.style.use('seaborn-whitegrid')

fig, ax = plt.subplots(figsize=(5, 4), dpi=500) #(7,6)
plt.yticks(np.arange(0, 7, 1))
plt.xticks(np.arange(0, 14, 1))
ax.set_xlabel('t')
ax.set_ylabel('T/Ts')
ax.set_xlim([0,13])
ax.set_ylim([0,6])
ax2 = ax
ax3 = ax

ax4 = ax
ax5 = ax
ax6 = ax

ax3.scatter(datos_theta100['time'], datos_theta100['temp'], marker='.', s=3)
ax2.scatter(datos_theta10['time'], datos_theta10['temp'], marker='.', s=3)
ax.scatter(datos_theta001['time'], datos_theta001['temp'], marker='.', s=3)

# Plot theoretical data
ax4.plot(teoricos_100['time'], teoricos_100['temp'], linewidth=0.5)
ax4.plot(teoricos_10['time'], teoricos_10['temp'], linewidth=0.5)
ax4.plot(teoricos_001['time'], teoricos_001['temp'], linewidth=0.5)
# Plot MD simulations 
# =============================================================================
# ax4.scatter(antonio_100['time'], antonio_100['temp'], marker='x', s=6, color='b')
# ax5.scatter(antonio_10['time'], antonio_10['temp'], marker='x', s=6, color='orange')
# ax6.scatter(antonio_001['time'], antonio_001['temp'], marker='x', s=6, color='g')
# =============================================================================

