#pragma once
#include <string>
#include <list>
#include "../ModelNonsense.h"
#include "KargvaModel.h"

using namespace std;

class KargvaModelNonsense : public ModelNonsense, public KargvaModel {
private:
    void makeModel(string line);
    string addContext();
public:
    KargvaModelNonsense(string line);
    void addToModel(string line);
    bool includes(string line);
    string condensedSNPinfo();
};

KargvaModelNonsense::KargvaModelNonsense(string line) {
    makeModel(line);
}
void KargvaModelNonsense::addToModel(string line) {
    throw std::exception("should not have been called: model type is nonsense");
}
bool KargvaModelNonsense::includes(string line) {
    return false;
}
string KargvaModelNonsense::addContext() :KargvaModel.addContext(pos) {}
string KargvaModelNonsense::condensedSNPinfo() {
    return make_pair(make_pair(wt_aa, pos), mutant_aa);
}
void KargvaModelNonsense::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1, line.size() - 2));
}