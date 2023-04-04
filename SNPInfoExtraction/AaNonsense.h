#pragma once
#include <string>
#include "ModelSNP.h"

using namespace std;

class AaNonsense : public virtual ModelSNP {
protected:
    char wt_aa = 0;
    virtual void makeModel(string line) = 0;
public:
    AaNonsense();
    list<int> getFirstPos();
    virtual ~AaNonsense();
    virtual string condensedInfo() = 0;
};

AaNonsense::AaNonsense() {}
AaNonsense::~AaNonsense() {}
list<int> AaNonsense::getFirstPos() {
    return {5, pos};				//nonsense come after insertion(4) and before multiple(6)
}