#pragma once
#include "../ModelMissense.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureModelMissense : public ModelMissense, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureModelMissense(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureModelMissense();
    LiteratureModelMissense(const LiteratureModelMissense& other);
    Model* Clone();
    string condensedSNPinfo();
    int getPos();
};

LiteratureModelMissense::LiteratureModelMissense(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, false) {
    makeModel(line);
}
LiteratureModelMissense::LiteratureModelMissense(const LiteratureModelMissense& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->mutant_aa = other.mutant_aa;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
}
Model* LiteratureModelMissense::Clone() {
    return new LiteratureModelMissense(*this);
}

int LiteratureModelMissense::getPos() {
    return pos;
}

LiteratureModelMissense::~LiteratureModelMissense() {}

string LiteratureModelMissense::condensedSNPinfo() {
    string toReturn = "Mis:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    for (int i = 0; i < mutant_aa.size(); i++)
        toReturn += mutant_aa[i];
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
void LiteratureModelMissense::makeModel(string line) {
    wt_aa = line.at(0);
    line = line.substr(1);
    int i = 0;
    while (isdigit(line.at(i))) {
        i += 1;
    }
    pos = stoi(line.substr(0, i));
    line = line.substr(i);
    for (int j = 0; j < line.length(); j++) {
        mutant_aa.push_back(line.at(j));
    }

}
