#ifndef SIPROPERTIES_H
#define SIPROPERTIES_H

class SiProperties {
public:

    SiProperties();
    virtual ~SiProperties();

    //  methods
    void Init(); 
    void PrintInfo();

    double Density();
    
    double Eloss(double mom, double mass, double tcut);
    double ElossVar(double mom, double mass);


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
};

#endif