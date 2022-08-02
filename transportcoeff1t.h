#ifndef TRANSPORTCOEFF1TBINARY_H
#define TRANSPORTCOEFF1TBINARY_H

#include <o2omixturedata.h>

struct TranspCoeff_o2{
    double shear_viscosity;
    double bulk_viscosity;
    double thermal_conductivity;
};

struct TranspCoeff_o2o{
    double shear_viscosity_o2;
    double bulk_viscosity_o2;
    double thermal_conductivity_o2;
    double shear_viscosity_o;
    double bulk_viscosity_o;
    double thermal_conductivity_o;
    double shear_viscosity;
    double bulk_viscosity;
    double thermal_conductivity;
    double diffusion_o2o2;
    double diffusion_o2o;
    double diffusion_oo;
};

TranspCoeff_o2 calculate_o2transportcoeff(double temperature, o2mixturedata data, OmegaInt_o2 omegaint);
TranspCoeff_o2o calculate_o2otransportcoeff(double temperature, double mass_fraction_o2, o2mixturedata data, OmegaInt_o2 omegaint);

#endif // TRANSPORTCOEFF1TBINARY_H
