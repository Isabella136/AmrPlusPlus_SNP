#pragma once
#include "../ModelMissense.h"
#include "KargvaModel.h"

using namespace std;

class KargvaModelMissense : public ModelMissense, public KargvaModel {
private:
    void makeModel(string line);
    bool includesAll(string line);
    bool includesMt(string snp);
public:
    KargvaModelMissense(string line, string id, shared_ptr<CARD_database> dbSeq);
    ~KargvaModelMissense();
    KargvaModelMissense(const KargvaModelMissense& other);
    Model* Clone();
    void addToModel(string line);
    bool includes(string line);
    string condensedSNPinfo();
};

KargvaModelMissense::KargvaModelMissense(string line, string id, shared_ptr<CARD_database> dbSeq):KargvaModel(id, dbSeq) {
    makeModel(line);
}
KargvaModelMissense::KargvaModelMissense(const KargvaModelMissense& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->mutant_aa = other.mutant_aa;
    this->cardID = other.cardID;
    this->databaseSequences = other.databaseSequences;
}
Model* KargvaModelMissense::Clone() {
    return new KargvaModelMissense(*this);
}
KargvaModelMissense::~KargvaModelMissense() {}
void KargvaModelMissense::addToModel(string line) {
    mutant_aa.push_back(line.at(line.size() - 1));
}
bool KargvaModelMissense::includesMt(string snp) {
    for (auto iter = this->mutant_aa.begin(); iter != this->mutant_aa.end(); ++iter) {
        if (snp.at(snp.size() - 1) == *iter)
            return true;
    }
    return false;
}
bool KargvaModelMissense::includesAll(string line) {
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
bool KargvaModelMissense::includes(string line) {
    if (line.find(";") != -1 || line.find(",") != -1)
        return this->includesAll(line);
    string temp = "";
    temp += wt_aa;
    temp += to_string(pos);
    return temp == line.substr(0, line.size() - 1);
}
string KargvaModelMissense::condensedSNPinfo() {
    string toReturn = "Mis:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    for (int i = 0; i < mutant_aa.size(); i++)
        toReturn += mutant_aa[i];
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
void KargvaModelMissense::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1, line.size() - 2));
    mutant_aa.push_back(line.at(line.size() - 1));
}