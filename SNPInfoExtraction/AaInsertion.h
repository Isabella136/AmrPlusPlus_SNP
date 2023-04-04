#pragma once
#include <string>
#include "ModelInDel.h"

using namespace std;

class AaInsertion : public virtual ModelInDel {
protected:
    string mt_aa = "";
    virtual void makeModel(string line) = 0;
public:
    AaInsertion();
    virtual ~AaInsertion();
    list<int> getFirstPos();
    virtual string condensedInfo() = 0;
};

AaInsertion::AaInsertion() {}
AaInsertion::~AaInsertion() {}
list<int> AaInsertion::getFirstPos() {
    return {4, *(pos.begin())};  //insertion come after deletion(3) and before nonsense(5)
}