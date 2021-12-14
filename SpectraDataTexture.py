# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


#!/usr/bin/env python
# vim:set ft=python fileencoding=utf-8 sr et ts=4 sw=4 : See help 'modeline'
import numpy as np
'''
    == A few notes about color ==

    Color   Wavelength(nm) Frequency(THz)
    Red     620-750        484-400
    Orange  590-620        508-484
    Yellow  570-590        526-508
    Green   495-570        606-526
    Blue    450-495        668-606
    Violet  380-450        789-668

    f is frequency (cycles per second)
    l (lambda) is wavelength (meters per cycle)
    e is energy (Joules)
    h (Plank's constant) = 6.6260695729 x 10^-34 Joule*seconds
                         = 6.6260695729 x 10^-34 m^2*kg/seconds
    c = 299792458 meters per second
    f = c/l
    l = c/f
    e = h*f
    e = c*h/l

    List of peak frequency responses for each type of
    photoreceptor cell in the human eye:
        S cone: 437 nm
        M cone: 533 nm
        L cone: 564 nm
        rod:    550 nm in bright daylight, 498 nm when dark adapted.
                Rods adapt to low light conditions by becoming more sensitive.
                Peak frequency response shifts to 498 nm.

'''

import sys
import os
import traceback
import optparse
import time
import logging
import math
from matplotlib import pyplot as plt
from astropy.io import fits
from scipy import interpolate
from PIL import Image
def wavelength_to_rgb(wavelength, gamma=0.8):

    '''This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''

    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    # R *= 255
    # G *= 255
    # B *= 255
    return ((R), (G), (B))



# Press the green button in the gutter to run the script.
pi = np.pi
h = 6.626e-34
c = 3.0e+8
k = 1.38e-23

def planck(wav, T):
    a = 2.0*h*pi*c**2
    b = h*c/(wav*k*T)
    intensity = a/ ( (wav**5)*(math.e**b - 1.0) )
    return intensity

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
range = 750-380

wav = np.linspace(380,750,1024)
# flux = planck((wav/10**9), 5700)
# flux = flux*(1/flux.max())
#
# spec = []
# for i,lam in enumerate(wav):
#     R,G,B = wavelength_to_rgb(lam,0.8)
#     R = R*flux[i]
#
#     G = G*flux[i]
#     B = B*flux[i]
#     spec.append([R,G,B])
# spec = np.array(spec)
#
# spec=  np.rint(spec).astype('int')
# Image = np.repeat(spec[ np.newaxis,:,:],1024,axis= 0)

# plt.imshow(Image, interpolation='nearest')
# plt.show()


with open('SpectraFits/extractlist.txt') as f:
    lines = f.readlines()
wav = np.linspace(380,780,4096)
SpecIm = []
for line in lines:
    try:
        file = line.removesuffix('\n')
        print(file)
        data = np.loadtxt('SpectraFits/'+file)
        wave = data[:,0]
        Flux = data[:,1]
        spec = interpolate.interp1d(wave, Flux)
        NewFlux = spec(wav*10)
        NewFlux = NewFlux/NewFlux.max()
        SpecIm.append(NewFlux)
        SpecIm.append(NewFlux)
        SpecIm.append(NewFlux)

    except:
        print("An exception occurred")

Input = np.array(SpecIm)


im = Image.fromarray((Input*255).astype(np.uint8))
im.save('TestData.png')
im.show()