#pragma once
#include <string>
#include "ModelInDel.h"

using namespace std;

class AaInsertion : public virtual ModelInDel {
protected:
    char mt_aa = 0;
    virtual void makeModel(string line) = 0;
public:
    AaInsertion();
    virtual ~AaInsertion();
    virtual string condensedInfo() = 0;
};

AaInsertion::AaInsertion() {}
AaInsertion::~AaInsertion() {}