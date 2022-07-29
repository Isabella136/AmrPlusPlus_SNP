#pragma once
#include <string>
#include "ModelSNP.h"

using namespace std;

class AaNonstop : public virtual ModelSNP {
protected:
    char mt_aa = 0;
    virtual void makeModel(string line) = 0;
public:
    AaNonstop();
    virtual ~AaNonstop();
    virtual string condensedInfo() = 0;
};

AaNonstop::AaNonstop() {}
AaNonstop::~AaNonstop() {}