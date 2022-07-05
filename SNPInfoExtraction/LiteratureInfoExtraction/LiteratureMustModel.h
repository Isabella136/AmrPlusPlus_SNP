#pragma once
#include <string>
#include "../ModelMust.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureMustModel : public ModelMust, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureMustModel(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureMustModel();
    LiteratureMustModel(const LiteratureMustModel& other);
    Model* Clone();
    string condensedSNPinfo();
    int getPos();
};

LiteratureMustModel::LiteratureMustModel(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, false) {
    makeModel(line);
}
LiteratureMustModel::LiteratureMustModel(const LiteratureMustModel& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
}
Model* LiteratureMustModel::Clone() {
    return new LiteratureMustModel(*this);
}

int LiteratureMustModel::getPos() {
    return pos;
}

LiteratureMustModel::~LiteratureMustModel() {}

string LiteratureMustModel::condensedSNPinfo() {
    string toReturn = "Mis:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
void LiteratureMustModel::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1));
}