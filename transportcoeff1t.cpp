#include "transportcoeff1t.h"

TranspCoeff_o2 calculate_o2transportcoeff(double temperature, o2mixturedata data, OmegaInt_o2 omegaint){
    TranspCoeff_o2 result;
    result.shear_viscosity = 5*kB*temperature/(8*omegaint.omega22_o2o2_LJ);
    result.thermal_conductivity = 15*kB*result.shear_viscosity/(4*massO2_particle) + 3*kB*temperature/(8*omegaint.omega11_o2o2_LJ)*(data.specific_heat_rot + data.specific_heat_vibr);
    result.bulk_viscosity = 0; // not calculated now
    return result;
}

//not finished, look for exceptions

TranspCoeff_o2o calculate_o2otransportcoeff(double temperature, double density, double mass_fraction_o2, o2omixturedata data, OmegaInt_o2o omegaint){
    TranspCoeff_o2o result;
    result.shear_viscosity_o2 = 5*kB*temperature/(8*omegaint.omega22_o2o2_LJ);
    result.thermal_conductivity_o2 = 15*kB*result.shear_viscosity_o2/(4*massO2_particle) + 3*kB*temperature/(8*omegaint.omega11_o2o2_LJ)*(data.specific_heat_rot + data.specific_heat_vibr);
    result.bulk_viscosity_o2 = 0; // not calculated now

    result.shear_viscosity_o = 5*kB*temperature/(8*omegaint.omega22_oo_LJ);
    result.thermal_conductivity_o = 15*kB*result.shear_viscosity_o/(4*massO_particle);
    result.bulk_viscosity_o2 = 0; // not calculated now

    double d_o2o = 3*kB*kB*data.molar_mass_total*temperature*(massO_particle + massO2_particle)/(16*density*UniversalGasConstant*massO_particle*massO2_particle*omegaint.omega11_o2o_LJ);
    result.diffusion_o2o2 = d_o2o*std::pow(UniversalGasConstant/(kB*data.molar_mass_total),2)*massO_particle*massO2_particle*(1-mass_fraction_o2)/mass_fraction_o2;
    result.diffusion_o2o = -d_o2o*std::pow(UniversalGasConstant/(kB*data.molar_mass_total),2)*massO_particle*massO2_particle;
    result.diffusion_oo = d_o2o*std::pow(UniversalGasConstant/(kB*data.molar_mass_total),2)*massO_particle*massO2_particle*mass_fraction_o2/(1-mass_fraction_o2);

    double phi_o2o = sqrt(2)*pow( 1 + sqrt(result.shear_viscosity_o2/result.shear_viscosity_o)*pow(massO_molecular/massO2_molecular,0.25) ,2)/(4*sqrt(1 + massO2_molecular/massO_molecular));
    double phi_oo2 = sqrt(2)*pow( 1 + sqrt(result.shear_viscosity_o/result.shear_viscosity_o2)*pow(massO2_molecular/massO_molecular,0.25) ,2)/(4*sqrt(1 + massO_molecular/massO2_molecular));
    double g_o2o = 0.3765343/sqrt( 1 + massO2_molecular/massO_molecular)*pow(1 + sqrt(result.shear_viscosity_o2/result.shear_viscosity_o)*pow(massO_molecular/massO2_molecular,0.25) , 2);
    double g_oo2 = 0.3765343/sqrt( 1 + massO_molecular/massO2_molecular)*pow(1 + sqrt(result.shear_viscosity_o/result.shear_viscosity_o2)*pow(massO2_molecular/massO_molecular,0.25) , 2);

    result.shear_viscosity = result.bulk_viscosity_o2/(1 + phi_o2o*massO2_particle*(1-mass_fraction_o2)/massO_particle/mass_fraction_o2) + result.bulk_viscosity_o/(1 + phi_oo2*massO_particle*mass_fraction_o2/massO2_particle/(1 - mass_fraction_o2));

    result.thermal_conductivity = result.thermal_conductivity_o2/(1 + g_o2o*massO_particle*mass_fraction_o2/massO2_particle/(1 - mass_fraction_o2)) + result.thermal_conductivity_o/(1 + g_oo2*massO2_particle*(1 - mass_fraction_o2)/massO_particle/mass_fraction_o2);

    if (mass_fraction_o2 == 0.0){
        result.bulk_viscosity = result.bulk_viscosity_o;
    }
    if (mass_fraction_o2 == 1.0){
        result.bulk_viscosity = result.bulk_viscosity_o2;
    }
    else {
        result.bulk_viscosity = pow(pow(result.bulk_viscosity_o2,3/4) + pow(result.bulk_viscosity_o,3/4)*massO2_particle*(1-mass_fraction_o2)/massO_particle/mass_fraction_o2,4/3);
    }

    return result;
}
