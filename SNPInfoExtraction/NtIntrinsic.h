#pragma once
#include <string>
#include "ModelSNP.h"

using namespace std;

class NtIntrinsic : public virtual ModelSNP {
protected:
    char wt_nuc = 0;
    virtual void makeModel(string line) = 0;
public:
    NtIntrinsic();
    virtual ~NtIntrinsic();
    virtual string condensedInfo() = 0;
};

NtIntrinsic::NtIntrinsic() {}
NtIntrinsic::~NtIntrinsic() {}