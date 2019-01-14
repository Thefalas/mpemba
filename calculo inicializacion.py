# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 19:03:56 2018

@author: malopez
"""

import pandas as pd

filePath = './Datos/sim0023/RgradsHCS.dat'

columnas = ['foto', 'cols_per_particle', 'time', 'temp_tras', 
            'temp_rot', 'theta', 'temp']
# theta = temp_rot/temp_tras
# temp = 0.5*(temp_tras + temp_rot)

datos = pd.read_csv(filePath, sep='\t', names=columnas)
datos.drop('foto', axis=1, inplace=True)

# =============================================================================
# fig, ax = plt.subplots(figsize=(3, 3), dpi=300)
# plt.scatter(x=datos['cols_per_particle'], y=datos['temp'])
# =============================================================================

# Stationary values
temp_stationary = datos['temp'][-2000:].mean()
temp_tras_stationary = datos['temp_tras'][-2000:].mean()
temp_rot_stationary = datos['temp_rot'][-2000:].mean()
theta_stationary = datos['theta'][-2000:].mean()

# Initial values
temp_0 = datos['temp'][0]
temp_tras_0 = datos['temp_tras'][0]
temp_rot_0 = datos['temp_rot'][0]
theta_0 = datos['theta'][0]

# =============================================================================
# multiples = 3
# theta_ratio = 100
# =============================================================================

# Calculation of future initialization values
def printInitValues(multiples, theta_ratio):
    future_t0 = multiples * temp_stationary
    future_tr0 = 2*theta_ratio*future_t0/(1+theta_ratio)
    future_tt0 = 2*future_t0 - future_tr0
    
    print('theta = ', theta_ratio, ' multiples = ', multiples)
    print('Stationary temperature = ', temp_stationary)
    print('Future TR0 = ', future_tr0, ' Future TT0 = ', future_tt0)
    print('')
    
printInitValues(5, 0.01)
printInitValues(4, 10)
printInitValues(3, 100)

# =============================================================================
# # Assuming theta = temp_rot/temp_tras
# b_future_t0 = multiples * temp_stationary
# b_future_tr0 = 2*theta_ratio*b_future_t0/(1+theta_ratio)
# #b_future_tr0 = (2*multiples*theta_ratio/(1+theta_ratio)) * temp_stationary
# b_future_tt0 = 2*b_future_t0 - b_future_tr0
# =============================================================================


# =============================================================================
# # Assuming theta = temp_rot/temp
# future_t0 = multiples * temp_stationary
# future_tr0 = theta_ratio * future_t0 
# future_tt0 = 2*future_t0 - future_tr0
# =============================================================================


# =============================================================================
# datos2 = pd.read_pickle('C:/Users/malopez/Desktop/mpemba/theta001mult5.pkl')
# plt.scatter(datos2['time'], datos2['temp'], color='r')
# plt.scatter(datos['time'], datos['temp'], color='g')
# =============================================================================


# =============================================================================
# pd.to_pickle(datos, './datosthetax.pkl')
# =============================================================================
