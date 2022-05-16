#pragma once
#include <string>
#include "ModelSingle.h"

using namespace std;

class ModelDeletion : public virtual ModelSingle {
protected:
    char wt_aa = 0;
    virtual void makeModel(string line) = 0;
public:
    ModelDeletion();
    virtual ~ModelDeletion();
    virtual string condensedSNPinfo() = 0;
};

ModelDeletion::ModelDeletion() {}
ModelDeletion::~ModelDeletion() {}