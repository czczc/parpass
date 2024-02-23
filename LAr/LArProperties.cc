#include "LArProperties.h"

#include <iostream>
using namespace std;

LArProperties::LArProperties()
    : Properties()
{
    // https://pdg.lbl.gov/2020/AtomicNuclearProperties/index.html

    // Following parameters are for use in Bethe-Bloch formula for dE/dx.
        fZ = 18;                ///< Ar atomic number
        fA = 39.948;            ///< Ar atomic mass (g/mol)
        fI = 188.00;            ///< Ar mean excitation energy (eV)

       fSa = 0.1956;            ///< Ar Sternheimer parameter a
       fSk = 3.0000;            ///< Ar Sternheimer parameter k
      fSx0 = 0.2000;            ///< Ar Sternheimer parameter x0
      fSx1 = 3.0000;            ///< Ar Sternheimer parameter x1
    fScbar = 5.2146;            ///< Ar Sternheimer parameter Cbar
    fSdelta0 = 0.;            ///< Sternheimer parameter delta0

    fTemperature = 87.8;     // K
}

LArProperties::~LArProperties()
{}

double LArProperties::Density()
{
    // temperature is assumed to be in degrees Kelvin
    // density is nearly a linear function of temperature.  see the NIST tables for details
    // slope is between -6.2 and -6.1, intercept is 1928 kg/m^3
    // this parameterization will be good to better than 0.5%.
    // density is returned in g/cm^3
    return -0.00615*Temperature() + 1.928;
}

//----------------------------------------------------------------------------------
void LArProperties::PrintInfo()
{
    double mass = 0.10565837;
    cout << "       Temperature [K]: " << Temperature() << endl;
    cout << "        Density [g/cc]: " << Density() << endl;
    cout << "dE/dx MIP mu [MeV/cm]: " << Eloss(MOM(0.25, mass), mass, 0) << endl;
    cout << "  MPV 5GeV mu [MeV/cm]: " << MPV(5, mass, 1.) << endl;
    cout << " dE Fluctu. [MeV^2/cm]: " << ElossVar(1, mass) << endl;

}
