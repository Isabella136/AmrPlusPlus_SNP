#pragma once
#include <string>
#include "ModelSNP.h"

using namespace std;

class AaNonstop : public virtual ModelSNP {
protected:
    virtual void makeModel(string line) = 0;
public:
    AaNonstop();
    virtual ~AaNonstop();
    virtual string condensedInfo() = 0;
};

AaNonstop::AaNonstop() {}
AaNonstop::~AaNonstop() {}