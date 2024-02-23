#ifndef PROPERTIES_H
#define PROPERTIES_H

class Properties {
public:

    Properties();
    virtual ~Properties();

    //  methods
    virtual double Density();
    virtual void PrintInfo();
    
    double Temperature()  { return fTemperature; } 
    double Eloss(double mom, double mass, double tcut);
    double ElossVar(double mom, double mass);
    double MPV(double mom, double mass, double thickness);
    double FWHM(double mom, double mass);

    double KE(double mom, double mass);
    double MOM(double ke, double mass);



    // Following parameters are for use in Bethe-Bloch formula for dE/dx.
    double fZ;                ///< Ar atomic number
    double fA;                ///< Ar atomic mass (g/mol)
    double fI;                ///< Ar mean excitation energy (eV)
    double fSa;               ///< Sternheimer parameter a
    double fSk;               ///< Sternheimer parameter k
    double fSx0;              ///< Sternheimer parameter x0
    double fSx1;              ///< Sternheimer parameter x1
    double fScbar;            ///< Sternheimer parameter Cbar
    double fSdelta0;          ///< Sternheimer parameter delta0

    double fTemperature; // density coulbe be temperature dependent
};

#endif