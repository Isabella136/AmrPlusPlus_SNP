#pragma once
#include "../ModelDeletion.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureModelDeletion : public ModelDeletion, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureModelDeletion(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureModelDeletion();
    LiteratureModelDeletion(const LiteratureModelDeletion& other);
    Model* Clone();
    string condensedSNPinfo();
    int getPos();
};

LiteratureModelDeletion::LiteratureModelDeletion(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, false) {
    makeModel(line);
}
LiteratureModelDeletion::LiteratureModelDeletion(const LiteratureModelDeletion& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
}
Model* LiteratureModelDeletion::Clone() {
    return new LiteratureModelDeletion(*this);
}

int LiteratureModelDeletion::getPos() {
    return pos;
}

LiteratureModelDeletion::~LiteratureModelDeletion() {}

string LiteratureModelDeletion::condensedSNPinfo() {
    string toReturn = "Del:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    toReturn += '-';
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
void LiteratureModelDeletion::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1, line.size() - 2));
}