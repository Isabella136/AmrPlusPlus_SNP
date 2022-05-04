#pragma once
#include "../ModelReg.h"
#include "KargvaModel.h"

using namespace std;

class KargvaModelReg : public ModelReg, public KargvaModel {
private:
    void makeModel(string line);
public:
    KargvaModelReg(string line, string id, CARD_database* dbSeq);
    void addToModel(string line);
    bool includes(string line);
    string condensedSNPinfo();
};

KargvaModelReg::KargvaModelReg(string line, string id, CARD_database* dbSeq):KargvaModel(id, dbSeq) {
    makeModel(line);
}
void KargvaModelReg::addToModel(string line) {
    mutant_aa.push_back(line.at(line.size() - 1));
}
bool KargvaModelReg::includes(string line) {
    string temp = "";
    temp += wt_aa;
    temp += to_string(pos);
    return temp == line.substr(0, line.size() - 1);
}
string KargvaModelReg::condensedSNPinfo() {
    string toReturn = "Reg:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    for (int i = 0; i < mutant_aa.size(); i++)
        toReturn += mutant_aa[i];
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
void KargvaModelReg::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1, line.size() - 2));
    mutant_aa.push_back(line.at(line.size() - 1));
}