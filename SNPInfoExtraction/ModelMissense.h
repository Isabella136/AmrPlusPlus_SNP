#pragma once
#include <string>
#include <vector>
#include "ModelSingle.h"

using namespace std;

class ModelMissense : public virtual ModelSingle {
    protected:
        char wt_aa = 0;
        vector<char> mutant_aa;
        virtual void makeModel(string line) = 0;
    public:
        ModelMissense();
        virtual ~ModelMissense();
        virtual string condensedSNPinfo() = 0;
};

ModelMissense::ModelMissense() {}
ModelMissense::~ModelMissense() {}