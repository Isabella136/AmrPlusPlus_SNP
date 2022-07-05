#pragma once
#include <string>
#include "../ModelNucMust.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureNucMustModel : public ModelNucMust, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureNucMustModel(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureNucMustModel();
    LiteratureNucMustModel(const LiteratureNucMustModel& other);
    Model* Clone();
    string condensedSNPinfo();
    int getPos();
};

LiteratureNucMustModel::LiteratureNucMustModel(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, true) {
    makeModel(line);
}
LiteratureNucMustModel::LiteratureNucMustModel(const LiteratureNucMustModel& other) {
    this->wt_nuc = other.wt_nuc;
    this->pos = other.pos;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
}
Model* LiteratureNucMustModel::Clone() {
    return new LiteratureNucMustModel(*this);
}

int LiteratureNucMustModel::getPos() {
    return pos;
}

LiteratureNucMustModel::~LiteratureNucMustModel() {}

string LiteratureNucMustModel::condensedSNPinfo() {
    string toReturn = "Nuc:";
    toReturn += wt_nuc;
    toReturn += to_string(pos);
    toReturn += addContext(pos - 2, pos, wt_nuc);
    return toReturn;
}
void LiteratureNucMustModel::makeModel(string line) {
    wt_nuc = line.at(0);
    pos = stoi(line.substr(1));
}