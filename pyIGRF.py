#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pyIGRF: code to synthesise magnetic field values from the 13th generation of the
        International Geomagnetic Reference Field, released in December 2020

 @author: Ciaran Beggan (British Geological Survey)
 
 See https://www.ngdc.noaa.gov/IAGA/vmod/ for information on the IGRF
 
 Based on existing codes: igrf13.f (FORTRAN) and chaosmagpy (Python3)
 
 With acknowledgements to: Clemens Kloss (DTU Space), David Kerridge (BGS) and
     Ashley Smith (Univ. of Edinburgh), william Brown and Grace Cox.
 
     This is a program for synthesising geomagnetic field values from the 
     International Geomagnetic Reference Field series of models as agreed
     in December 2019 by IAGA Working Group V-MOD. 
     
     It is the 13th generation IGRF, ie the 12th revision. 
     The main-field models for 1900.0, 1905.0,..1940.0 and 2020.0 are 
     non-definitive, those for 1945.0, 1950.0,...2015.0 are definitive and
     the secular-variation model for 2020.0 to 2025.0 is non-definitive.

     Main-field models are to degree and order 10 (ie 120 coefficients)
     for 1900.0-1995.0 and to 13 (ie 195 coefficients) for 2000.0 onwards. 
     The predictive secular-variation model is to degree and order 8 (ie 80
     coefficients).

     Options include values at different locations at different
     times (spot), values at same location at one year intervals
     (time series), grid of values at one time (grid); geodetic or
     geocentric coordinates, latitude & longitude entered as decimal
     degrees or degrees & minutes (not in grid), and outputs main field 
     or secular variation or both (grid only).

Dependencies: numpy, scipy

 Recent history of code:
     Initial release: April 2020 (Ciaran Beggan, BGS)
 
"""
from scipy import interpolate
import igrf_utils as iut
import io_options as ioo

IGRF_FILE = r'./IGRF13.shc'
igrf = iut.load_shcfile(IGRF_FILE, None)


if __name__ == '__main__':
    # Introduction text and initial option selection
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    print(' ')
    print('******************************************************')
    print('*              IGRF SYNTHESIS PROGRAM                *')
    print('*                                                    *')
    print('* A program for the computation of geomagnetic       *')
    print('* field elements from the International Geomagnetic  *')
    print('* Reference Field (13th generation) as revised in    *')
    print('* December 2019 by the IAGA Working Group V-MOD.     *')
    print('*                                                    *')
    print('* It is valid for dates from 1900.0 to 2025.0,       *')
    print('* values up to 2030.0 will be computed with          *')
    print('* reduced accuracy. Values for dates before 1945.0   *')
    print('* and after 2015.0 are non-definitive, otherwise the *')
    print('* values are definitive.                             *')
    print('*                                                    *')
    print('*                                                    *')
    print('*            (on behalf of) IAGA Working Group V-MOD *')
    print('******************************************************')
    print(' ')
    print('Enter name of output file')
    print('or press "Return" for output to screen')
      
    name = input("Enter filename: ")
    
    if not name:         # is empty 
        print('Printing to screen')
    else:
        print(name)
     
    
    while 1:
        print( 'Choose an option:')
        print( '1 - values at one or more locations & dates')
        print( '2 - values at yearly intervals at one location')
        print( '3 - values on a latitude/longitude grid at one date')
        iopt = input('->')
        iopt = iut.check_int(iopt)
        if iopt < 1 or iopt > 3: continue
        else:
            break
        
    # Parse the inputs for computing the main field and SV values. 
    # Convert geodetic to geocentric coordinates if required 
    if iopt == 1:
        date, alt, lat, colat, lon, itype, sd, cd = ioo.option1()
    elif iopt ==2:
         date, alt, lat, colat, lon, itype, sd, cd = ioo.option2()
    else:
         date, alt, lat, colat, lon, itype, sd, cd = ioo.option3()
        
        
    # Interpolate the geomagnetic coefficients to the desired date(s)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    f = interpolate.interp1d(igrf.time, igrf.coeffs, fill_value='extrapolate')
    coeffs = f(date)    
    
    # Compute the main field B_r, B_theta and B_phi value for this location 
    Br, Bt, Bp = iut.synth_values(coeffs.T, alt, colat, lon,
                              igrf.parameters['nmax'])
    
    # For the SV, find the 5 year period in which the date lies and compute
    # the SV within that period. IGRF has constant SV between each 5 year period
    # We don't need to subtract 1900 but it makes it clearer:
    epoch = (date-1900)//5    
    epoch_start = epoch*5
    # Add 1900 back on plus 1 year to account for SV in nT per year (nT/yr):
    coeffs_sv = f(1900+epoch_start+1) - f(1900+epoch_start)   
    Brs, Bts, Bps = iut.synth_values(coeffs_sv.T, alt, colat, lon,
                              igrf.parameters['nmax'])
    
    # Use the main field coefficients from the start of each five epch
    # to compute the SV for Dec,Inc,Hor and Total Field (F) 
    # [Note: these are non-linear components of X, Y and Z so treat separately]
    coeffsm = f(1900+epoch_start);
    Brm, Btm, Bpm = iut.synth_values(coeffsm.T, alt, colat, lon,
                              igrf.parameters['nmax'])
    
    
    # Rearrange to X, Y, Z components 
    X = -Bt; Y = Bp; Z = -Br
    # For the SV
    dX = -Bts; dY = Bps; dZ = -Brs 
    Xm = -Btm; Ym = Bpm; Zm = -Brm
    # Rotate back to geodetic coords if needed
    if (itype == 1):
        t = X; X = X*cd + Z*sd;  Z = Z*cd - t*sd
        t = dX; dX = dX*cd + dZ*sd;  dZ = dZ*cd - t*sd
        t = Xm; Xm = Xm*cd + Zm*sd;  Zm = Zm*cd - t*sd
        
    # Compute the fourn non-linear components 
    dec, hoz, inc, eff = iut.xyz2dhif(X,Y,Z)
    # The IGRF SV coefficients are relative to the main field components 
    # at the start of each five year epoch e.g. 2010, 2015, 2020
    decs, hozs, incs, effs = iut.xyz2dhif_sv(Xm, Ym, Zm, dX, dY, dZ)
    
    
    # Finallt, parse the outputs for writing to screen or file
    if iopt == 1:
        ioo.write1(name, date, alt, lat, colat, lon, X, Y, Z, dX, dY, dZ, \
                  dec, hoz, inc, eff, decs, hozs, incs, effs, itype)
        if name:
            print('Written to file: ' + name )
    elif iopt == 2:
        ioo.write2(name, date, alt, lat, colat, lon, X, Y, Z, dX, dY, dZ, \
                  dec, hoz, inc, eff, decs, hozs, incs, effs, itype)
        if name:
            print('Written to file: ' + name )
    else:
        ioo.write3(name, date, alt, lat, colat, lon, X, Y, Z, dX, dY, dZ, \
                  dec, hoz, inc, eff, decs, hozs, incs, effs, itype)
        if name:
            print('Written to file: ' + name )
    