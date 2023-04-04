#pragma once
#include <string>
#include "ModelInDel.h"

using namespace std;

class AaDeletion : public virtual ModelInDel {
protected:
    char wt_aa = 0;
    virtual void makeModel(string line) = 0;
public:
    AaDeletion();
    virtual ~AaDeletion();
    virtual string condensedInfo() = 0;
};

AaDeletion::AaDeletion() {}
AaDeletion::~AaDeletion() {}