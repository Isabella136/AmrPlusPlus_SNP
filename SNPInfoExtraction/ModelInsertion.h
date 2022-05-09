/*#pragma once
#include <string>
#include <vector>
#include "ModelSingle.h"

using namespace std;

class ModelInsertion : public virtual ModelSingle {
protected:
    vector<char> mutant_aa;
    virtual void makeModel(string line) = 0;
public:
    ModelInsertion();
    virtual string condensedSNPinfo() = 0;
};

ModelInsertion::ModelInsertion() {}*/