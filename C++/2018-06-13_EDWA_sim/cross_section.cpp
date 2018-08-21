#include "cross_section.h"   


//Returns erbium absorption cross section from 1450 to 1600 and 962 to 990 nm
double Simulation::erbium_abs(int const x_) noexcept 
{
    double const x {static_cast<double>(x_) * 1.0e-3};
    if (x >= 1450.0 && x < 1524.59) {
        constexpr std::array<double, 5> coef {
                                    9.735520850960E-08,
                                    -5.733291173420E-04,
                                    1.266032348160E+00,
                                    -1.242399633250E+03,
                                    4.571533145180E+05
                                    };

        return (coef[0] * pow(x, 4) + coef[1] * pow(x, 3) + coef[2] * pow(x, 2) + coef[3] * x
                                    + coef[4]) * 1.0e-25;

    }

    if (x >= 1524.59 && x < 1535.50) {
        constexpr std::array<double, 2> coef {
                                    3.254002865100E-01,
                                    -4.936677289960E+02
                                    };

        return (coef[0] * x + coef[1]) * 1.0e-25;
    }

    if (x >= 1535.50 && x < 1545.09) {
        constexpr std::array<double, 3> coef {
                                    2.481273818570E-02,
                                    -7.680909314830E+01,
                                    5.944391474140E+04
                                    };

        return (coef[0] * pow(x, 2) + coef[1] * x + coef[2]) * 1.0e-25;
    }

    if (x >= 1545.09 && x < 1549.14) {
        constexpr std::array<double, 2> coef {
                                    3.703703703700E-02,
                                    -5.474555555560E+01,
                                    };

        return (coef[0] * x + coef[1]) * 1.0e-25;
    }

    if (x >= 1549.14 && x <= 1600.5) {
        constexpr std::array<double, 4> coef {
                                    -2.104690961080E-05,
                                    1.008731048010E-01,
                                    -1.611558226700E+02,
                                    8.582251098830E+04,
                                    };

        return (coef[0] * pow(x, 3) + coef[1] * pow(x, 2) + coef[2] * x + coef[3]) * 1.0e-25;
    }

    if (x >= 962.0 && x < 976.0) {
        constexpr std::array<double, 2> coef {
                                    0.152027027025,
                                    -145.384864865
                                    };

        return (coef[0] * x + coef[1]) * 1.0e-25;
    }

    if (x >= 976.0 && x <= 990.0) {
        constexpr std::array<double, 2> coef {
                                    -0.133461538464,
                                    133.132564103
                                    };

        return (coef[0] * x + coef[1]) * 1.0e-25;
    }

    std::cout<<"Error invalid wavelength range in erbium_abs(): x = "<<x<<", returning -1.0\n";
    return -1.0;
}

//Erbium emission cross section from 1450 nm to 1600 nm in m^2//aquiaqui
double Simulation::erbium_emi(int const x_) noexcept 
{
    double const x {static_cast<double>(x_) * 1.0e-3};
    if (x >= 1450.0 && x < 1524.59) {
        constexpr std::array<double, 5> coef {
                                    2.160448531710E-07,
                                    -1.274295973020E-03,
                                    2.818613327600E+00,
                                    -2.770917365500E+03,
                                    1.021518706090E+06,
                                    };

        return (coef[0] * pow(x, 4) + coef[1] * pow(x, 3) + coef[2] * pow(x, 2) + coef[3] * x
                                    + coef[4]) * 1.0e-25;
    }

    if (x >= 1524.59 && x < 1535.50) {
        constexpr std::array<double, 3> coef {
                                    2.764411482320E-02,
                                    -8.411410261070E+01,
                                    6.398673948220E+04
                                    };

        return (coef[0] * pow(x, 2) + coef[1] * x + coef[2]) * 1.0e-25;
    }

    if (x >= 1535.50 && x < 1543.61) {
        constexpr std::array<double, 3> coef {
                                    3.700996665420E-02,
                                    -1.144608410390E+02,
                                    8.850172318370E+04
                                    };

        return (coef[0] * pow(x, 2) + coef[1] * x + coef[2]) * 1.0e-25;
    }

    if (x >= 1543.61 && x < 1550.98) {
        constexpr std::array<double, 2> coef {
                                    1.071913161470E-01,
                                    -1.618115875170E+02
                                    };

        return (coef[0] * x + coef[1]) * 1.0e-25;
    }

    if (x >= 1550.98 && x <= 1600.5) {
        constexpr std::array<double, 4> coef {
                                    -3.226607677830E-05,
                                    1.546350413260E-01,
                                    -2.470411464540E+02,
                                    1.315626807900E+05
                                    };

        return (coef[0] * pow(x, 3) + coef[1] * pow(x, 2) + coef[2] * x + coef[3]) * 1.0e-25;
    }

    std::cout<<"Error invalid wavelength range in erbium_emi(): x = "<<x<<", returning -1.0\n";
    return -1.0;
}


double Simulation::ytterbium_abs(int const x_) noexcept 
{
    double const x {static_cast<double>(x_) * 1.0e-3};
    if (x >= 850.0 && x < 856.98) {
        return 0.0;
    }

    if (x >= 856.98 && x < 966.70) {
        constexpr std::array<double, 5> coef {
                                    1.171073941740E-04,
                                    -4.249066925600E-01,
                                    5.775754324150E+02,
                                    -3.485875637330E+05,
                                    7.881627669790E+07
                                    };

        return (coef[0] * pow(x, 4) + coef[1] * pow(x,3) + coef[2] * pow(x,2) + coef[3] * x + coef[4]) * 1.0e-27;
    }

    if (x >= 966.70 && x < 973.99) {
        constexpr std::array<double, 2> coef {
                                    3.042441700960E+02,
                                    -2.937153992320E+05
                                    };

        return (coef[0] * x + coef[1]) * 1.0e-27;
    }

    if (x >= 973.99 && x < 983.17) {
        constexpr std::array<double, 2> coef {
                                    -2.998477366260E+02,
                                    2.946640769960E+05
                                    };

        return (coef[0] * x + coef[1]) * 1.0e-27;
    }

    if (x >= 983.17 && x < 1045.05) {
        constexpr std::array<double, 5> coef {
                                    9.568284931450E-05,
                                    -3.915462900430E-01,
                                    6.007627689730E+02,
                                    -4.096207838530E+05,
                                    1.047215458420E+08
                                    };

        return (coef[0] * pow(x, 4) + coef[1] * pow(x, 3) + coef[2] * pow(x, 2)  + coef[3] * x  + coef[4]) * 1.0e-27;
    }

    if (x >= 1045.05 && x <= 1150.5) {
        return 0.0;
    }

    std::cout<<"Error invalid wavelength range in ytterbium_abs(): x = "<<x<<", returning -1.0\n";
    return -1.0;
}


double Simulation::ytterbium_emi(int const x_) noexcept 
{
    double const x {static_cast<double>(x_) * 1.0e-3};
    
    if (x >= 850.0 && x < 886.46) {

        return 0.0;
    }

    if (x >= 886.46 && x < 967.33) {
        constexpr std::array<double, 5> coef {
                                    6.554178215310E-05,
                                    -2.408830756610E-01,
                                    3.319142623400E+02,
                                    -2.032175908870E+05,
                                    4.664708696060E+07
                                    };

        return (coef[0] * pow(x, 4) + coef[1] * pow(x,3) + coef[2] * pow(x, 2) + coef[3] * x + coef[4]) * 1.0e-27;
    }

    if (x >= 967.33 && x < 973.7) {
        constexpr std::array<double, 2> coef {
                                    3.539097978230E+02,
                                    -3.420655147280E+05
                                    };

        return (coef[0] * x + coef[1]) * 1.0e-27;
    }

    if (x >= 973.7 && x < 981.68) {
        constexpr std::array<double, 2> coef {
                                    -2.646654040400E+02,
                                    2.602782738380E+05
                                    };

        return (coef[0] * x + coef[1]) * 1.0e-27;
    }

    if (x >= 981.68 && x < 1053.56) {
        constexpr std::array<double, 5> coef {
                                    3.455662757690E-04,
                                    -1.415741483130E+00,
                                    2.174361844750E+03,
                                    -1.483749312730E+06,
                                    3.795629871630E+08
                                    };

        return (coef[0] * pow(x, 4) + coef[1] * pow(x, 3) + coef[2] * pow(x,2) + coef[3] * x + coef[4]) * 1.0e-27;
    }

    if (x >= 1053.56 && x <= 1150.5) {
        constexpr std::array<double, 2> coef {
                                    -3.559654146570E+00,
                                    4.088236107100E+03
                                    };

        return (coef[0] * x + coef[1]) * 1.0e-27;
    }

    std::cout<<"Error invalid wavelength range in ytterbium_emi(): x = "<<x<<", returning -1.0\n";
    return -1.0;
}