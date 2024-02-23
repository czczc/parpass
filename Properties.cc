#include "Properties.h"

#include <iostream>
using namespace std;

Properties::Properties()
{
    // https://pdg.lbl.gov/2020/AtomicNuclearProperties/index.html

    // use Si as default
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

Properties::~Properties()
{}


double Properties::Density()
{
    // use Si as default
    return 2.329 + 0.*Temperature();
}


//----------------------------------------------------------------------------------
// Restricted mean energy loss (dE/dx) in units of MeV/cm.
//
// For unrestricted mean energy loss, set tcut = 0, or tcut large.
//
// Arguments:
//
// mom  - Momentum of incident particle in GeV/c.
// mass - Mass of incident particle in GeV/c^2.
// tcut - Maximum kinetic energy of delta rays (MeV).
//
// Returned value is positive.
//
// Based on Bethe-Bloch formula as contained in particle data book.
// Material parameters (stored in Properties.fcl) are taken from
// pdg web site http://pdg.lbl.gov/AtomicNuclearProperties/
//----------------------------------------------------------------------------------
double Properties::Eloss(double mom, double mass, double tcut)
{
    // Some constants.
    double K = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
    double me = 0.510998918; // Electron mass (MeV/c^2).

    // Calculate kinematic quantities.
    double bg = mom / mass;           // beta*gamma.
    double gamma = sqrt(1. + bg*bg);  // gamma.
    double beta = bg / gamma;         // beta (velocity).
    double mer = 0.001 * me / mass;   // electron mass / mass of incident particle.
    double tmax = 2.*me* bg*bg / (1. + 2.*gamma*mer + mer*mer);  // Maximum delta ray energy (MeV).
    // Make sure tcut does not exceed tmax.
    if(tcut == 0. || tcut > tmax) { tcut = tmax; }

    // Calculate density effect correction (delta).
    double x = log10(bg);
    double delta = 0.;
    if(x >= fSx0) {
        delta = 2. * log(10.) * x - fScbar;
        if(x < fSx1) {
            delta += fSa * pow(fSx1 - x, fSk);
        }
    }

    // Calculate stopping number.
    double B = 0.5 * log(2.*me*bg*bg*tcut / (1.e-12 * fI*fI))
      - 0.5 * beta*beta * (1. + tcut / tmax) - 0.5 * delta;
    // Don't let the stopping number become negative.
    if(B < 1.) B = 1.;

    // Calculate dE/dx.
    double dedx = Density() * K*fZ*B / (fA * beta*beta);
    // cout << fZ << ", " << fA << ", " << Density() << endl;
    return dedx;   // MeV/cm.
}

//----------------------------------------------------------------------------------
// Energy loss fluctuation (sigma_E^2 / length in MeV^2/cm).
//
// Arguments:
//
// mom  - Momentum of incident particle in GeV/c.
//
// Based on Bichsel formula referred to but not given in pdg.
//----------------------------------------------------------------------------------
double Properties::ElossVar(double mom, double mass)
{
    // Some constants.
    double K = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
    double me = 0.510998918; // Electron mass (MeV/c^2).

    // Calculate kinematic quantities.
    double bg = mom / mass;          // beta*gamma.
    double gamma2 = 1. + bg*bg;      // gamma^2.
    double beta2 = bg*bg / gamma2;   // beta^2.

    // Calculate final result.
    double result = gamma2 * (1. - 0.5 * beta2) * me * (fZ / fA) * K * Density();
    return result;
}


//----------------------------------------------------------------------------------
// Most probable energy loss (dE/dx) in units of MeV/cm.
//
// Arguments:
//
// mom  - Momentum of incident particle in GeV/c.
// mass - Mass of incident particle in GeV/c^2.
// thickness - thickness of material in cm.
//----------------------------------------------------------------------------------
double Properties::MPV(double mom, double mass, double thickness)
{
    // Some constants.
    double K = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
    double me = 0.510998918; // Electron mass (MeV/c^2).

    // Calculate kinematic quantities.
    double bg = mom / mass;           // beta*gamma.
    double gamma = sqrt(1. + bg*bg);  // gamma.
    double beta = bg / gamma;         // beta (velocity).
    double mer = 0.001 * me / mass;   // electron mass / mass of incident particle.

    // Calculate density effect correction (delta).
    double x = log10(bg);
    double delta = 0.;
    if(x >= fSx0) {
        delta = 2. * log(10.) * x - fScbar;
        if(x < fSx1) {
            delta += fSa * pow(fSx1 - x, fSk);
        }
    }

    // Calculate MPV
    double zeta = K/2 * fZ/fA * thickness*Density()/beta/beta;
    double mpv = zeta * (
        log(2.*me*bg*bg/ (1.e-6 * fI))
        + log(zeta/(1.e-6 * fI))
        + 0.2 - beta*beta - delta
    ) / thickness;
    return mpv; // MeV/cm

}

double Properties::FWHM(double mom, double mass)
{
    // Some constants.
    double K = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
    double me = 0.510998918; // Electron mass (MeV/c^2).

    // Calculate kinematic quantities.
    double bg = mom / mass;           // beta*gamma.
    double gamma = sqrt(1. + bg*bg);  // gamma.
    double beta = bg / gamma;         // beta (velocity).
    // double mer = 0.001 * me / mass;   // electron mass / mass of incident particle.

    // Calculate density effect correction (delta).
    // double x = log10(bg);
    // double delta = 0.;
    // if(x >= fSx0) {
    //     delta = 2. * log(10.) * x - fScbar;
    //     if(x < fSx1) {
    //         delta += fSa * pow(fSx1 - x, fSk);
    //     }
    // }

    // Calculate MPV
    double zeta = K/2 * fZ/fA *Density()/beta/beta; // MeV/cm
    return 4*zeta;
}


//----------------------------------------------------------------------------------
void Properties::PrintInfo()
{
    double mass = 0.10565837;
    cout << "        Density [g/cc]: " << Density() << endl;
    cout << "dE/dx MIP mu [MeV/cm]: " << Eloss(0.3633, 0.10565837, 0) << endl;
    cout << "dE Fluctu. [MeV^2/cm]: " << ElossVar(0.3633, 0.10565837) << endl;
    cout << "4*zeta [MeV]: " << FWHM(0.3633, 0.10565837) << endl;
    cout << "MPV for MIP mu at 35 um thickness [MeV/cm]: " << MPV(0.3633, 0.10565837, 35e-4) << endl;
    cout << "MPV for MIP mu at 120 um thickness [MeV/cm]: " << MPV(0.3633, 0.10565837, 120e-4) << endl;
    cout << "MPV for MIP mu at 200 um thickness [MeV/cm]: " << MPV(0.3633, 0.10565837, 200e-4) << endl;

    // cout << MOM(4e-3, mass) << endl;
    // cout << KE(0.3633, mass) << endl;
    cout << Eloss(MOM(4e-3, mass), mass, 0) << endl;
    cout << MPV(MOM(4e-3, mass), mass, 120e-4) << endl;
    cout << FWHM(MOM(4e-3, mass), mass) << endl;

    cout << Eloss(MOM(1e-3, mass), mass, 0) << endl;
    cout << MPV(MOM(1e-3, mass), mass, 120e-4) << endl;
    cout << FWHM(MOM(1e-3, mass), mass) << endl;

    // cout << Eloss(0.3633, mass, 0) << endl;
    // cout << MPV(0.3633, mass, 120e-4) << endl;
    // cout << FWHM(0.3633, mass) << endl;

}

double Properties::KE(double mom, double mass)
{
    double bg = mom / mass;           // beta*gamma.
    double gamma = sqrt(1. + bg*bg);  // gamma.
    return (gamma-1) * mass;  // KE in GeV
}

double Properties::MOM(double ke, double mass)
{
    return sqrt(ke*ke+2*ke*mass);
}
