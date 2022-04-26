#pragma once
#include <string>
#include <list>
#include "../ModelReg.h"
#include "KargvaModel.h"

using namespace std;

class KargvaModelReg : public ModelReg, public KargvaModel {
private:
    void makeModel(string line);
    string addContext();
public:
    KargvaModelReg(string line);
    void addToModel(string line);
    bool includes(string line);
    string condensedSNPinfo();
};

KargvaModelReg::KargvaModelReg(string line) {
    makeModel(line);
}
void KargvaModelReg::addToModel(string line) {
    mutant_aa.push_back(line.at(line.size() - 1));
}
bool KargvaModelReg::includes(string line) {
    string temp = "";
    temp += wt_aa;
    temp += pos;
    return temp == line.substr(0, line.size() - 1);
}
string KargvaModelReg::condensedSNPinfo() {
    return make_pair(make_pair(wt_aa, pos), mutant_aa);
}
void KargvaModelReg::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1, line.size() - 2));
    mutant_aa.push_back(line.at(line.size() - 1));
}