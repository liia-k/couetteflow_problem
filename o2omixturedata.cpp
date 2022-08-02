#include "o2omixturedata.h"

std::array<double,num_vibr_level_o2> vibr_exp_index(double temperature){
    std::array<double,num_vibr_level_o2> result;
    for (size_t i = 0; i < result.size(); i++){
        result[i] = exp(-vibr_energy_o2[i]/(UniversalGasConstant*temperature));
    }
    return result;
}

o2mixturedata calculate_o2mixturedata(double temperature){
    o2mixturedata result;
    std::array<double,num_vibr_level_o2> vibr_exp = vibr_exp_index(temperature);
    double s1 = 0;
    double s2 = 0;
    double s3 = 0;
    for (size_t i = 0; i < vibr_exp.size(); i++){
        s1 += vibr_energy_o2[i]*vibr_exp[i];
        s2 += vibr_exp[i];
        s3 += (vibr_energy_o2[i])*(vibr_energy_o2[i])*vibr_exp[i];
    }

    result.energy_vibr_o2 = (s1/s2 + vibr_enrgy_o2_0_level)/massO2_molecular;
    result.energy_int_o2 = result.energy_vibr_o2 + kB*temperature/massO2_particle;
    result.energy_mixture = result.energy_int_o2 + 3*kB*temperature/massO2_particle/2;

    result.specific_heat_rot = kB/massO2_particle;
    result.specific_heat_vibr = (s3/s2 - (s1/s2)*(s1/s2))/(massO2_molecular*UniversalGasConstant*temperature);
    result.specific_heat_mixture = result.specific_heat_rot + result.specific_heat_vibr + 3*kB/(2*massO2_particle);

    result.Gamma = 1 + UniversalGasConstant/(massO2_molecular*result.specific_heat_mixture);
    return result;
}

o2omixturedata calculate_o2omixturedata(double temperature, double mass_fraction_o2)
{
    o2omixturedata result;
    std::array<double,num_vibr_level_o2> vibr_exp = vibr_exp_index(temperature);
    double s1 = 0;
    double s2 = 0;
    double s3 = 0;
    for (size_t i = 0; i < vibr_exp.size(); i++){
        s1 += vibr_energy_o2[i]*vibr_exp[i];
        s2 += vibr_exp[i];
        s3 += (vibr_energy_o2[i])*(vibr_energy_o2[i])*vibr_exp[i];
    }

    result.energy_vibr_o2 = mass_fraction_o2*(s1/s2 + vibr_enrgy_o2_0_level)/massO2_molecular;
    result.energy_int_o2 = result.energy_vibr_o2 + kB*temperature/massO2_particle;
    result.energy_mixture = result.energy_int_o2 + 3*kB*temperature/massO2_particle/2 + 3*kB*temperature/massO_particle/2 + energy_formation_o() + energy_formation_o2();

    result.enthalpy_o = 5*kB*temperature/massO_particle/2 + energy_formation_o();
    result.enthalpy_o2 = result.energy_int_o2 + energy_formation_o2() + 5*kB*temperature/massO2_particle/2;

    result.specific_heat_rot = kB/massO2_particle;
    result.specific_heat_vibr = mass_fraction_o2*(s3/s2 - (s1/s2)*(s1/s2))/(massO2_molecular*UniversalGasConstant*temperature);
    result.specific_heat_o = 3*kB/(2*massO_particle);
    result.specific_heat_o2 = result.specific_heat_rot + result.specific_heat_vibr + 3*kB/(2*massO2_particle);
    result.specific_heat_mixture = result.specific_heat_o + result.specific_heat_o2;

    result.molar_mass_total = 1/(mass_fraction_o2/massO2_molecular + (1-mass_fraction_o2)/massO_molecular);
    result.Gamma = 1 + UniversalGasConstant/(result.molar_mass_total*result.specific_heat_mixture);
    return result;
}

OmegaInt_o2 calculate_omega_int_o2(double temperature){
    OmegaInt_o2 result;
    double omega_rs11 = sqrt(M_PI*kB*temperature/(massO2_particle))*molecule_dist_o2o2*molecule_dist_o2o2;
    double omega_rs22 = sqrt(4*M_PI*kB*temperature/(massO2_particle))*molecule_dist_o2o2*molecule_dist_o2o2;
    double x11 = log(kB*temperature/min_value_potential_o2o2) + a11_LJ;
    double x22 = log(kB*temperature/min_value_potential_o2o2) + a22_LJ;
    result.omega11_o2o2_LJ = omega_rs11/(f11[0] + f11[1]/(x11*x11) + f11[2]/x11 + f11[3]*x11 + f11[4]*x11*x11 + f11[5]*std::pow(x11,3));
    result.omega22_o2o2_LJ = omega_rs22/(f22[0] + f22[1]/(x22*x22) + f22[2]/x22 + f22[3]*x22 + f22[4]*x22*x22 + f22[5]*std::pow(x22,3));
    return result;
}

OmegaInt_o2o calculate_omega_int_o2o(double temperature){
    OmegaInt_o2o result;
    double omega_rs11_o2o2 = sqrt(M_PI*kB*temperature/(massO2_particle))*molecule_dist_o2o2*molecule_dist_o2o2;
    double omega_rs22_o2o2 = sqrt(4*M_PI*kB*temperature/(massO2_particle))*molecule_dist_o2o2*molecule_dist_o2o2;
    double omega_rs11_oo = sqrt(M_PI*kB*temperature/(massO_particle))*molecule_dist_oo*molecule_dist_oo;
    double omega_rs22_oo = sqrt(4*M_PI*kB*temperature/(massO_particle))*molecule_dist_oo*molecule_dist_oo;
    double omega_rs11_o2o = sqrt(M_PI*kB*temperature*(massO2_particle + massO_particle)/(2*massO2_particle*massO_particle))*molecule_dist_o2o*molecule_dist_o2o;

    double x11_o2o2 = log(kB*temperature/min_value_potential_o2o2) + a11_LJ;
    double x22_o2o2 = log(kB*temperature/min_value_potential_o2o2) + a22_LJ;
    double x11_oo = log(kB*temperature/min_value_potential_oo) + a11_LJ;
    double x22_oo = log(kB*temperature/min_value_potential_oo) + a22_LJ;
    double x11_o2o = log(kB*temperature/min_value_potential_o2o) + a11_LJ;

    result.omega11_o2o2_LJ = omega_rs11_o2o2/(f11[0] + f11[1]/(x11_o2o2*x11_o2o2) + f11[2]/x11_o2o2 + f11[3]*x11_o2o2 + f11[4]*x11_o2o2*x11_o2o2 + f11[5]*std::pow(x11_o2o2,3));
    result.omega22_o2o2_LJ = omega_rs22_o2o2/(f22[0] + f22[1]/(x22_o2o2*x22_o2o2) + f22[2]/x22_o2o2 + f22[3]*x22_o2o2 + f22[4]*x22_o2o2*x22_o2o2 + f22[5]*std::pow(x22_o2o2,3));
    result.omega11_oo_LJ = omega_rs11_oo/(f11[0] + f11[1]/(x11_oo*x11_oo) + f11[2]/x11_oo + f11[3]*x11_oo + f11[4]*x11_oo*x11_oo + f11[5]*std::pow(x11_oo,3));
    result.omega22_oo_LJ = omega_rs22_oo/(f22[0] + f22[1]/(x22_oo*x22_oo) + f22[2]/x22_oo + f22[3]*x22_oo + f22[4]*x22_oo*x22_oo + f22[5]*std::pow(x22_oo,3));
    result.omega11_o2o_LJ = omega_rs11_o2o/(f11[0] + f11[1]/(x11_o2o*x11_o2o) + f11[2]/x11_o2o + f11[3]*x11_o2o + f11[4]*x11_o2o*x11_o2o + f11[5]*std::pow(x11_o2o,3));

    return result;
}
