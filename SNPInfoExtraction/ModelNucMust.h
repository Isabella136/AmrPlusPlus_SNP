#pragma once
#include <string>
#include "ModelSingle.h"

using namespace std;

class ModelNucMust : public virtual ModelSingle {
protected:
    char wt_nuc = 0;
    virtual void makeModel(string line) = 0;
public:
    ModelNucMust();
    virtual ~ModelNucMust();
    virtual string condensedSNPinfo() = 0;
};

ModelNucMust::ModelNucMust() {}
ModelNucMust::~ModelNucMust() {}