#ifndef HEADERS_H_INCLUDED
#define HEADERS_H_INCLUDED

//Libraries
#include <iostream>

//Project files
#include "edwa_simulation.h"




/*
Conventions:
------------
=Coding
 * All constants are upper case
 * Classes and structs and namespaces begin with the first letter in upper case
 * New types are all in lower case
 * Variables are going to be in the same case as in the book
 * n5, n6 are Ytterbium levels, while n1, n2, n3, n4 are erbium levels
 * The numbers 12 in a switch statement means the change from level 1 to 2 whereas 21 is from 2 to 1
 of erbium levels
 * Frequency will be abbreviated "freq", wavelength "wl" 
 * Braced initialization will be used where possible and "=" assignment will be used when changing
 an already existing variable
 * If a function call has many variables and can't fit in a single line, each variable will be 
 split in a line with 1 tab separation
 E.G.:
    double const flux {Utility::return_photon_flux(
        wl1,
        pump_pow1,
        p.A,
        56
    )};
 
=Simulation
 * Steady state assumed
 * Erbium ions absorb 976 but not emit it
 * Ytterbium ions absorb and emit 976 but not 1500-1600
 * No ASE signals in the 900-1000 nm region
 * 1 Forward and 1 Backward pump will be assumed
 * Single amplifier configuration
 * Signal enters from the Forward position
 * Forward ASE is 0 at the beginning, Backwards ASE is 0 at the End
 * Units are in mks (metres, watts, seconds, etc)
 * ASE signals going to be used from 1500 to 1600, number specified at the start
 * Fixed number of Signals specified at the start
 * All quantities will be in these units: Hertz, meters, seconds, Watts, Joules with the exception
 of map indexes which will be picometres for convenience (wl_map and overlap)
 * All Erbium and Yttebium concentration related scattering losses are assumed to be the same for
 different wavelengths (i.e signal scattering loss  due to Er at 1520 nm is the same as 1600 nm)
 * 

=Issues/Notes
 * Incorporate cross relaxation (Ccr) and Upconversion (Cup) functions to calculate those values
based on the erbium/ytterbium concentration
 * Refractive index is assumed to be 1 in wl to freq conversions, this is not true
 * Incorporate refractive index and cross section calculation before the main calculation
*/


#endif  /* HEADERS_H_INCLUDED */