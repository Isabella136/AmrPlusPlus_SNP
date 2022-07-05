#pragma once
#include <string>
#include "ModelSingle.h"

using namespace std;

class ModelMust : public virtual ModelSingle {
protected:
    char wt_aa = 0;
    virtual void makeModel(string line) = 0;
public:
    ModelMust();
    virtual ~ModelMust();
    virtual string condensedSNPinfo() = 0;
};

ModelMust::ModelMust() {}
ModelMust::~ModelMust() {}