#pragma once
#include <string>
#include <list>
#include "ModelSingle.h"

using namespace std;

class ModelReg : public virtual ModelSingle {
    protected:
        char wt_aa = 0;
        list<char> mutant_aa;
        virtual void makeModel(string line) = 0;
    public:
        ModelReg();
        virtual string condensedSNPinfo() = 0;
};

ModelReg::ModelReg() {}