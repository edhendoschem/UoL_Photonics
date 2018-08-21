#ifndef CROSS_SECTION_H_INCLUDED
#define CROSS_SECTION_H_INCLUDED

#include <iostream>
#include <array>
#include <cmath>


namespace Simulation
{

    //Returns erbium absorption cross section from 1450 to 1600 nm in m^2. Input in picometres
    double erbium_abs(int const x) noexcept;

    //Erbium emission cross section from 1450 to 1600 nm in m^2. Input in pm
    double erbium_emi(int const x) noexcept;

    //Ytterbium absortpion cross section from 850 to 1150 nm in m^2. Input in pm
    double ytterbium_abs(int const x) noexcept;

    //Ytterbium emission cross section from 850 to 1150 nm in m^2. Input in pm
    double ytterbium_emi(int const x) noexcept;

} //End of namespace simulation

#endif