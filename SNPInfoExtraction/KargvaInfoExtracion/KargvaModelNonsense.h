#pragma once
#include "../ModelNonsense.h"
#include "KargvaModel.h"

using namespace std;

class KargvaModelNonsense : public ModelNonsense, public KargvaModel {
private:
    void makeModel(string line);
public:
    KargvaModelNonsense(string line, string id, shared_ptr<CARD_database> dbSeq);
    ~KargvaModelNonsense();
    KargvaModelNonsense(const KargvaModelNonsense& other);
    Model* Clone();
    void addToModel(string line);
    bool includes(string line);
    string condensedSNPinfo();
};

KargvaModelNonsense::KargvaModelNonsense(string line, string id, shared_ptr<CARD_database> dbSeq) :KargvaModel(id, dbSeq) {
    makeModel(line);
}
KargvaModelNonsense::KargvaModelNonsense(const KargvaModelNonsense& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->cardID = other.cardID;
    this->databaseSequences = other.databaseSequences;
}
Model* KargvaModelNonsense::Clone() {
    return new KargvaModelNonsense(*this);
}
KargvaModelNonsense::~KargvaModelNonsense() {}

void KargvaModelNonsense::addToModel(string line) {
    throw std::exception("should not have been called: model type is nonsense");
}
bool KargvaModelNonsense::includes(string line) {
    return false;
}
string KargvaModelNonsense::condensedSNPinfo() {
    string toReturn = "Non:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    toReturn += '*';
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
void KargvaModelNonsense::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1, line.size() - 5));
}