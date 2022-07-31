#pragma once
#include "KargvaModel.h"
#include "../AaInsertion.h"

using namespace std;

class KargvaAaInsertion : public KargvaModel, public AaInsertion {
    private:
        void makeModel(string line);
    public:
        KargvaAaInsertion(string line, string id, shared_ptr<CARD_database> dbSeq);
        ~KargvaAaInsertion();
        KargvaAaInsertion(const KargvaAaInsertion& other);
        InfoPipe* Clone();
        void addToModel(string line);
        bool includes(string line);
        string condensedInfo();
        string infoType();
};

KargvaAaInsertion::KargvaAaInsertion(string line, string id, shared_ptr<CARD_database> dbSeq) : KargvaModel(id, dbSeq) {
    makeModel(line);
}
KargvaAaInsertion::~KargvaAaInsertion() {}
KargvaAaInsertion::KargvaAaInsertion(const KargvaAaInsertion& other) {
    this->cardID = other.cardID;
    this->databaseSequences = other.databaseSequences;
    this->mt_aa = other.mt_aa;
    this->pos = other.pos;
}
InfoPipe* KargvaAaInsertion::Clone() {
    return new KargvaAaInsertion(*this);
}
void KargvaAaInsertion::makeModel(string line) {
    int i = 0;
    for (i = 0; isalpha(line[i]); i++) {
        mt_aa += line[i];
    }
    line = line.substr(i, line.length()-i-1);
    while (line.find("/") != -1) {
        pos.push_back(stoi(line.substr(0, line.find("/"))));
        line = line.substr(line.find("/")+1);
    }
    pos.push_back(stoi(line));
}
void KargvaAaInsertion::addToModel(string line) {
    throw std::exception("should not have been called: model type is insertion");
}
bool KargvaAaInsertion::includes(string line) {
    return false;
}
string KargvaAaInsertion::condensedInfo() {
    string toReturn = "Ins:";
    toReturn += mt_aa;
    int firstPos = *(pos.begin());
    int lastPos = 0;
    for (auto iter = pos.begin(); iter != pos.end(); ++iter) {
        lastPos = *iter;
        if (lastPos != firstPos)
            toReturn += "/";
        toReturn += *iter;
    }
    toReturn += addContext(firstPos - 2, lastPos - 1, 0);
    return toReturn;
}
string KargvaAaInsertion::infoType() {
    return "Model";
}