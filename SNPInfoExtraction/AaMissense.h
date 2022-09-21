#pragma once
#include <string>
#include <vector>
#include "ModelSNP.h"

using namespace std;

class AaMissense : public virtual ModelSNP {
    protected:
        char wt_aa = 0;
        vector<char> mutant_aa;
        virtual void makeModel(string line) = 0;
    public:
        AaMissense();
        virtual ~AaMissense();
        virtual string condensedInfo() = 0;
};

AaMissense::AaMissense() {}
AaMissense::~AaMissense() {}