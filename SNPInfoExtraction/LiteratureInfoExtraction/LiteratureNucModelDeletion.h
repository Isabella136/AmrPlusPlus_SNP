#pragma once
#include "../ModelNucDeletion.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureNucModelDeletion : public ModelNucDeletion, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureNucModelDeletion(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureNucModelDeletion();
    LiteratureNucModelDeletion(const LiteratureNucModelDeletion& other);
    Model* Clone();
    string condensedSNPinfo();
    int getPos();
};

LiteratureNucModelDeletion::LiteratureNucModelDeletion(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, true) {
    makeModel(line);
}
LiteratureNucModelDeletion::LiteratureNucModelDeletion(const LiteratureNucModelDeletion& other) {
    this->wt_nuc = other.wt_nuc;
    this->pos = other.pos;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
}
Model* LiteratureNucModelDeletion::Clone() {
    return new LiteratureNucModelDeletion(*this);
}

int LiteratureNucModelDeletion::getPos() {
    return pos;
}

LiteratureNucModelDeletion::~LiteratureNucModelDeletion() {}

string LiteratureNucModelDeletion::condensedSNPinfo() {
    string toReturn = "NucDel:";
    toReturn += wt_nuc;
    toReturn += to_string(pos);
    toReturn += '-';
    toReturn += addContext(pos - 2, pos, wt_nuc);
    return toReturn;
}
void LiteratureNucModelDeletion::makeModel(string line) {
    wt_nuc = line.at(0);
    pos = stoi(line.substr(1, line.size() - 2));
}
