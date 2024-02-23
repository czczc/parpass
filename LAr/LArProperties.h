#ifndef LARPROPERTIES_H
#define LARPROPERTIES_H

#include "../Properties.h"

class LArProperties: public Properties {
public:

    LArProperties();
    virtual ~LArProperties();

    //  methods
    virtual double Density();

    virtual void PrintInfo();

};

#endif