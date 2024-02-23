#include "SiProperties.h"

#include <iostream>
using namespace std;

SiProperties::SiProperties()
    : Properties()
{
    // https://pdg.lbl.gov/2020/AtomicNuclearProperties/index.html

    // Following parameters are for use in Bethe-Bloch formula for dE/dx.
        fZ = 14;                ///< atomic number
        fA = 28.0855;            ///< atomic mass (g/mol)
        fI = 173.0;            ///< mean excitation energy (eV)

    // https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_silicon_Si.pdf
       fSa = 0.14921;            ///< Sternheimer parameter a
       fSk = 3.2546;            ///< Sternheimer parameter k
      fSx0 = 0.2015;            ///< Sternheimer parameter x0
      fSx1 = 2.8716;            ///< Sternheimer parameter x1
    fScbar = 4.4355;            ///< Sternheimer parameter Cbar
    fSdelta0 = 0.14;            ///< Sternheimer parameter delta0

    fTemperature = 293; // assume room temperature
}

SiProperties::~SiProperties()
{}


double SiProperties::Density()
{
    return 2.329;
}

// //----------------------------------------------------------------------------------
// void SiProperties::PrintInfo()
// {
//     double mass = 0.10565837;
//     cout << "        Density [g/cc]: " << Density() << endl;
//     cout << "dE/dx MIP mu [MeV/cm]: " << Eloss(0.3633, 0.10565837, 0) << endl;
//     cout << "dE Fluctu. [MeV^2/cm]: " << ElossVar(0.3633, 0.10565837) << endl;
//     cout << "4*zeta [MeV]: " << FWHM(0.3633, 0.10565837) << endl;
//     cout << "MPV for MIP mu at 35 um thickness [MeV/cm]: " << MPV(0.3633, 0.10565837, 35e-4) << endl;
//     cout << "MPV for MIP mu at 120 um thickness [MeV/cm]: " << MPV(0.3633, 0.10565837, 120e-4) << endl;
//     cout << "MPV for MIP mu at 200 um thickness [MeV/cm]: " << MPV(0.3633, 0.10565837, 200e-4) << endl;

//     // cout << MOM(4e-3, mass) << endl;
//     // cout << KE(0.3633, mass) << endl;
//     cout << Eloss(MOM(4e-3, mass), mass, 0) << endl;
//     cout << MPV(MOM(4e-3, mass), mass, 120e-4) << endl;
//     cout << FWHM(MOM(4e-3, mass), mass) << endl;

//     cout << Eloss(MOM(1e-3, mass), mass, 0) << endl;
//     cout << MPV(MOM(1e-3, mass), mass, 120e-4) << endl;
//     cout << FWHM(MOM(1e-3, mass), mass) << endl;

//     // cout << Eloss(0.3633, mass, 0) << endl;
//     // cout << MPV(0.3633, mass, 120e-4) << endl;
//     // cout << FWHM(0.3633, mass) << endl;

// }
