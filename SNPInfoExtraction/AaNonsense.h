#pragma once
#include <string>
#include "ModelSNP.h"

using namespace std;

class AaNonsense : public virtual ModelSNP {
protected:
    char wt_aa = 0;
    virtual void makeModel(string line) = 0;
public:
    AaNonsense();
    virtual ~AaNonsense();
    virtual string condensedInfo() = 0;
};

AaNonsense::AaNonsense() {}
AaNonsense::~AaNonsense() {}