#pragma once
#include "../AaNonsense.h"
#include "KargvaModel.h"

using namespace std;

class KargvaAaNonsense : public AaNonsense, public KargvaModel {
private:
    void makeModel(string line);
public:
    KargvaAaNonsense(string line, string id, shared_ptr<CARD_database> dbSeq);
    ~KargvaAaNonsense();
    KargvaAaNonsense(const KargvaAaNonsense& other);
    InfoPipe* Clone();
    void addToModel(string line);
    bool includes(string line);
    string condensedInfo();
    string infoType();
};

KargvaAaNonsense::KargvaAaNonsense(string line, string id, shared_ptr<CARD_database> dbSeq) :KargvaModel(id, dbSeq) {
    makeModel(line);
}
KargvaAaNonsense::KargvaAaNonsense(const KargvaAaNonsense& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->cardID = other.cardID;
    this->databaseSequences = other.databaseSequences;
}
InfoPipe* KargvaAaNonsense::Clone() {
    return new KargvaAaNonsense(*this);
}
KargvaAaNonsense::~KargvaAaNonsense() {}

void KargvaAaNonsense::addToModel(string line) {
    throw std::exception("should not have been called: model type is nonsense");
}
bool KargvaAaNonsense::includes(string line) {
    return false;
}
string KargvaAaNonsense::condensedInfo() {
    string toReturn = "Non:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    toReturn += '*';
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
void KargvaAaNonsense::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1, line.size() - 5));
}
string KargvaAaNonsense::infoType() {
    return "Model";
}