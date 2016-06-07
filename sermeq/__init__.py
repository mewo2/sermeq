"""Useful functions for glaciology.
"""
from __future__ import division

import numpy as np
from scipy.special import erf

#: Number of seconds in a day
day = 86400.
#: Number of seconds in a year
year = 365.25 * day

#: Specific heat capacity of ice (:math:`\mathrm{J}\,\mathrm{kg}^{-1}\,\mathrm{K}^{-1}`)
specific_heat_capacity = 2097
#: Latent heat of fusion of ice (:math:`\mathrm{J}\,\mathrm{kg}^{-1}`)
latent_heat = 333500.
#: Thermal diffusivity of ice (:math:`\mathrm{m}^2\,\mathrm{s}^{-1}`)
thermal_diffusivity = 1.09e-6
#: Thermal conductivity of ice (:math:`\mathrm{W}\,\mathrm{m}^{-1}\,\mathrm{K}^{-1}`)
thermal_conductivity = 2.1


zero_celsius = 273.15

R = 8.314
As = (3.615e-13, 1.733e3)
Qs = (6e4, 13.9e4)

def kelvin(celsius):
    """Convert temperature from Celsius to Kelvin
    """
    return celsius + zero_celsius

def celsius(kelvin):
    """Convert temperature from Kelvin to Celsius
    """
    return kelvin - zero_celsius

def flow_param(temperature):
    """Calculate the flow parameter :math:`A` for ice at a given temperature
    """
    A = As[temperature > -10]
    Q = Qs[temperature > -10]
    return A * np.exp(-Q/(R * kelvin(temperature)))

def temp_from_flow(A=None, B=None):
    """Calculate the ice temperature, given the flow parameter :math:`A` or :math:`B`
    """
    if A is not None and B is not None:
        raise Exception
    if A is None:
        A = B ** -3
    A_ = As[A > 4.3e-25]
    Q = Qs[A > 4.3e-25]
    return celsius(-Q/(R*np.log(A / A_)))

def water_frac_from_flow(A=None, B=None, T=0):
    if A is None:
        A = B ** -3
    ratio = A / flow_param(T)
    return (ratio - 1) / 181.25

def robin_temperature(height, thickness, surfacetemp, heatflux, accumulation):
    """Compute the steady-state temperature according to Robin (1955)
    """
    q = np.sqrt(accumulation / (2 * thermal_diffusivity * thickness))
    T = surfacetemp - (heatflux * np.sqrt(np.pi) / (2 * thermal_conductivity * q)
            * (erf(height * q) - erf(thickness * q)))
    return T
