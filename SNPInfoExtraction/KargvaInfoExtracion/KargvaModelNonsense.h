#pragma once
#include "../ModelNonsense.h"
#include "KargvaModel.h"

using namespace std;

class KargvaModelNonsense : public ModelNonsense, public KargvaModel {
private:
    void makeModel(string line);
public:
    KargvaModelNonsense(string line, string id);
    void addToModel(string line);
    bool includes(string line);
    string condensedSNPinfo();
};

KargvaModelNonsense::KargvaModelNonsense(string line, string id) :KargvaModel(id) {
    makeModel(line);
}
void KargvaModelNonsense::addToModel(string line) {
    throw std::exception("should not have been called: model type is nonsense");
}
bool KargvaModelNonsense::includes(string line) {
    return false;
}
string KargvaModelNonsense::condensedSNPinfo() {
    string toReturn = "Non:";
    toReturn += wt_aa;
    toReturn += pos;
    toReturn += '*';
    toReturn += addContext(pos - 2, pos);
    return toReturn;
}
void KargvaModelNonsense::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1, line.size() - 5));
}