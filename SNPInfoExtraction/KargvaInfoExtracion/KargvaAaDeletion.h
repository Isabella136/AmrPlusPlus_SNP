#pragma once
#include "../AaDeletion.h"
#include "KargvaModel.h"

using namespace std;

class KargvaAaDeletion : public AaDeletion, public KargvaModel {
private:
    void makeModel(string line);
public:
    KargvaAaDeletion(string line, string id, shared_ptr<CARD_database> dbSeq);
    ~KargvaAaDeletion();
    KargvaAaDeletion(const KargvaAaDeletion& other);
    InfoPipe* Clone();
    void addToModel(string line);
    bool includes(string line);
    string condensedInfo();
    string infoType();
};

KargvaAaDeletion::KargvaAaDeletion(string line, string id, shared_ptr<CARD_database> dbSeq) :KargvaModel(id, dbSeq) {
    makeModel(line);
}
KargvaAaDeletion::KargvaAaDeletion(const KargvaAaDeletion& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->cardID = other.cardID;
    this->databaseSequences = other.databaseSequences;
}
InfoPipe* KargvaAaDeletion::Clone() {
    return new KargvaAaDeletion(*this);
}
KargvaAaDeletion::~KargvaAaDeletion() {}
void KargvaAaDeletion::addToModel(string line) {
    throw std::exception("should not have been called: model type is deletion");
}
bool KargvaAaDeletion::includes(string line) {
    return false;
}
string KargvaAaDeletion::condensedInfo() {
    string toReturn = "Del:";
    toReturn += wt_aa;
    int firstPos = *(pos.begin());
    int lastPos = 0;
    for (auto iter = pos.begin(); iter != pos.end(); ++iter) {
        lastPos = *iter;
        if (lastPos != firstPos)
            toReturn += "/";
        toReturn += to_string(*iter);
    }
    toReturn += addContext(firstPos - 2, lastPos, wt_aa);
    return toReturn;
}
void KargvaAaDeletion::makeModel(string line) {
    if (isalpha(line.at(0))) {
        wt_aa = line.at(0);
        string posString = line.substr(1, line.size() - 2);
        while (posString.find("/") != -1) {
            pos.push_back(stoi(posString.substr(0, posString.find("/"))));
            posString = posString.substr(posString.find("/") + 1);
        }
        pos.push_back(stoi(posString.substr(0)));
    }
    else if (isalpha(line.at(1))) {
        wt_aa = line.at(1);
        pos.push_back(stoi(line.substr(2)));
    }
    else {
        pos.push_back(stoi(line.substr(1, line.size() - 2)));
        wt_aa = line.at(line.size() - 1);
    }
}
string KargvaAaDeletion::infoType() {
    return "Model";
}