#pragma once
#include "../ModelDeletion.h"
#include "KargvaModel.h"

using namespace std;

class KargvaModelDeletion : public ModelDeletion, public KargvaModel {
private:
    void makeModel(string line);
public:
    KargvaModelDeletion(string line, string id, shared_ptr<CARD_database> dbSeq);
    ~KargvaModelDeletion();
    KargvaModelDeletion(const KargvaModelDeletion& other);
    Model* Clone();
    void addToModel(string line);
    bool includes(string line);
    string condensedSNPinfo();
};

KargvaModelDeletion::KargvaModelDeletion(string line, string id, shared_ptr<CARD_database> dbSeq) :KargvaModel(id, dbSeq) {
    makeModel(line);
}
KargvaModelDeletion::KargvaModelDeletion(const KargvaModelDeletion& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->cardID = other.cardID;
    this->databaseSequences = other.databaseSequences;
}
Model* KargvaModelDeletion::Clone() {
    return new KargvaModelDeletion(*this);
}
KargvaModelDeletion::~KargvaModelDeletion() {}
void KargvaModelDeletion::addToModel(string line) {
    throw std::exception("should not have been called: model type is deletion");
}
bool KargvaModelDeletion::includes(string line) {
    return false;
}
string KargvaModelDeletion::condensedSNPinfo() {
    string toReturn = "Del:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    toReturn += '-';
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
void KargvaModelDeletion::makeModel(string line) {
    if (isalpha(line.at(0))) {
        wt_aa = line.at(0);
        pos = stoi(line.substr(1, line.size() - 2));
    }
    else if (isalpha(line.at(1))) {
        wt_aa = line.at(1);
        pos = stoi(line.substr(2));
    }
    else {
        pos = stoi(line.substr(1, line.size() - 2));
        wt_aa = line.at(line.size() - 1);
    }
}