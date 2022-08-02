#ifndef O2OMIXTUREDATA_H
#define O2OMIXTUREDATA_H

#include "global.h"
#include <array>

static const double diameterO2 = 1.32e-10; // molecular diameter in m
static const double diameterO = 0.6e-10; // atomic diameter in m
static const double massO2_molecular = 31.9988; //amu = gram/mol
static const double massO2_particle = massO2_molecular/Nav/1.0e3; // kg
static const double massO_molecular = 15.9994; // amu
static const double massO_particle = massO_molecular/Nav/1.0e3; // kg

struct macroParam
{
    double density             = 0;
    double massfraction_o2     = 0;
    double pressure            = 0;
    double velocity_tau        = 0;
    double velocity_n          = 0;
    double temp                = 0;
    double soundSpeed          = 0;
    bool isLeftContact         = false;
    QString gas                = "O2,O";
};

constexpr int num_vibr_level_o2 = 47; // including 0 level
constexpr double vibr_enrgy_o2_0_level = 0.09745;
constexpr double vibr_energy_o2[num_vibr_level_o2] = {0.09745, 0.29054, 0.48125, 0.66945, 0.85495, 1.03785, 1.21795, 1.39525, 1.56965, 1.74105, 1.90935, 2.07465, 2.23675, 2.39555, 2.55105, 2.70315, 2.85175, 2.99675, 3.13815, 3.2758, 3.40975, 3.53985, 3.66585, 3.78785, 3.90565, 4.01925, 4.12855, 4.2334, 4.33362, 4.42921, 4.52001, 4.6059, 4.68674, 4.76239, 4.83271, 4.89753, 4.95668, 5.00998, 5.05725, 5.0983, 5.132932, 5.161, 5.182427, 5.197353, 5.2063752, 5.21083, 5.21246}; // eV = J/mol
//constexpr double vibr_energy_o2[num_vibr_level_o2] = {-5.1153, -4.9221, -4.7315, -4.5433, -4.3578, -4.1749, -3.9948, -3.8175, -3.6431, -3.4717, -3.3034, -3.1381, -2.9760, -2.8172, -2.6617, -2.5096, -2.3610, -2.2160, -2.0746, -1.9369, -1.8030, -1.6729, -1.5469, -1.4249, -1.3071, -1.1935, -1.0842, -0.97939, -0.87913, -0.78354, -0.69274, -0.60685, -0.52601, -0.45036, -0.38004, -0.31522, -0.25607, -0.20277, -0.15550, -0.11445, -0.079818, -0.051751, -0.030323, -0.015397, -0.0063748, -0.0019261, -0.00029275}; // eV = J/mol

std::array<double,num_vibr_level_o2> vibr_exp_index(double temperature);

constexpr double energy_formation_o2(){
    return 0; // chemical reactions are not included yet
}
constexpr double energy_formation_o(){
    return 0;
}

struct o2mixturedata
{
    double energy_vibr_o2;
    double energy_int_o2;
    double energy_mixture;
    double specific_heat_vibr;
    double specific_heat_rot;
    double specific_heat_mixture;
    double Gamma;
};

struct o2omixturedata
{
    double energy_vibr_o2;
    double energy_int_o2;
    double energy_mixture;
    double enthalpy_o2;
    double enthalpy_o;
    double molar_mass_total;
    double specific_heat_vibr;
    double specific_heat_rot;
    double specific_heat_o2;
    double specific_heat_o;
    double specific_heat_mixture;
    double Gamma;
};

o2mixturedata calculate_o2mixturedata(double temperature);
o2omixturedata calculate_o2omixturedata(double temperature, double mass_fraction_o2);

constexpr double molecule_dist_o2o2 = 3.621e-10; // m
constexpr double molecule_dist_o2o = 3.185e-10;
constexpr double molecule_dist_oo = 2.75e-10;
constexpr double min_value_potential_o2o2 = 97.5*kB; // J
constexpr double min_value_potential_oo = 80.0*kB;
constexpr double min_value_potential_o2o = 88.075*kB; // wrong, find ~correct value
constexpr double a11_LJ = 1.4;
constexpr double a22_LJ = 1.5;
constexpr double f11[6] = {-0.16845, -0.02258, 0.19779, 0.64373, -0.09267, 0.00711};
constexpr double f22[6] = {-0.40811, -0.05086, 0.34010, 0.70375, -0.10699, 0.00763};

struct OmegaInt_o2{
    double omega11_o2o2_LJ;
    double omega22_o2o2_LJ;
};

struct OmegaInt_o2o{
    double omega11_o2o2_LJ;
    double omega11_o2o_LJ;
    double omega11_oo_LJ;
    double omega22_o2o2_LJ;
    double omega22_oo_LJ;
};

OmegaInt_o2 calculate_omega_int_o2(double temperature);
OmegaInt_o2o calculate_omega_int_o2o(double temperature);

#endif // O2OMIXTUREDATA_H
