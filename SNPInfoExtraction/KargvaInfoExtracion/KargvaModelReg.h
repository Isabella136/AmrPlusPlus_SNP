#pragma once
#include "../ModelReg.h"
#include "KargvaModel.h"

using namespace std;

class KargvaModelReg : public ModelReg, public KargvaModel {
private:
    void makeModel(string line);
    bool includesAll(string line);
    bool includesMt(string snp);
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
bool KargvaModelReg::includesMt(string snp) {
    for (auto iter = this->mutant_aa.begin(); iter != this->mutant_aa.end(); ++iter) {
        if (snp.at(snp.size() - 1) == *iter)
            return true;
    }
    return false;
}
bool KargvaModelReg::includesAll(string line) {
    string temp = line;
    vector<string> snp;
    if (line.find(";") != -1) {
        snp.push_back(temp.substr(0, temp.find(";")));
        while (temp.find(";") != -1) {
            temp = temp.substr(temp.find(";") + 1);
            snp.push_back(temp.substr(0, temp.find(";")));
        }
    }
    else {
        snp.push_back(temp.substr(0, temp.find(",")));
        while (temp.find(",") != -1) {
            temp = temp.substr(temp.find(",") + 1);
            snp.push_back(temp.substr(0, temp.find(",")));
        }
    }
    for (auto iter = snp.begin(); iter != snp.end(); ++iter) {
        if ((*iter).find(to_string(this->pos)) != -1) {
            if ((*iter).find('-') != -1)
                return false;
            else if (this->includesMt(*iter))
                return true;
            else
                return false;
        }
    }
    return false;

}
bool KargvaModelReg::includes(string line) {
    if (line.find(";") != -1 || line.find(",") != -1)
        return this->includesAll(line);
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