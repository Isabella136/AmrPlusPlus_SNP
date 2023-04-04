#pragma once
#include "KargvaModel.h"
#include "../AaIntrinsic.h"

using namespace std;

class KargvaAaIntrinsic : public KargvaModel, public AaIntrinsic {
    private:
        void makeModel(string line);
    public:
        KargvaAaIntrinsic(string line, string id, shared_ptr<CARD_database> dbSeq);
        ~KargvaAaIntrinsic();
        KargvaAaIntrinsic(const KargvaAaIntrinsic& other);
        InfoPipe* Clone();
        void addToModel(string line);
        bool includes(string line);
        int getPos();
        string infoType();
        string condensedInfo();
};
KargvaAaIntrinsic::KargvaAaIntrinsic(string line, string id, shared_ptr<CARD_database> dbSeq) : KargvaModel(id, dbSeq) {
    makeModel(line);
}
KargvaAaIntrinsic::~KargvaAaIntrinsic() {}
KargvaAaIntrinsic::KargvaAaIntrinsic(const KargvaAaIntrinsic& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->cardID = other.cardID;
    this->databaseSequences = other.databaseSequences;
}
InfoPipe* KargvaAaIntrinsic::Clone() {
    return new KargvaAaIntrinsic(*this);
}
void KargvaAaIntrinsic::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1));
}
void KargvaAaIntrinsic::addToModel(string line) {
    throw std::exception("should not have been called: model type is intrinsic");
}
bool KargvaAaIntrinsic::includes(string line) {
    return false;
}
int KargvaAaIntrinsic::getPos() {
    return pos;
}
string KargvaAaIntrinsic::infoType() {
    return "Model";
}
string KargvaAaIntrinsic::condensedInfo() {
    string toReturn = "Amino:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}