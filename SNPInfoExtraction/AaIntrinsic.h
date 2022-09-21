#pragma once
#include <string>
#include "ModelSNP.h"

using namespace std;

class AaIntrinsic : public virtual ModelSNP {
protected:
    char wt_aa = 0;
    virtual void makeModel(string line) = 0;
public:
    AaIntrinsic();
    virtual ~AaIntrinsic();
    virtual string condensedInfo() = 0;
};

AaIntrinsic::AaIntrinsic() {}
AaIntrinsic::~AaIntrinsic() {}