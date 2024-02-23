#ifndef SIPROPERTIES_H
#define SIPROPERTIES_H

#include "../Properties.h"

class SiProperties: public Properties {
public:

    SiProperties();
    virtual ~SiProperties();

    //  methods
    virtual double Density();
};

#endif