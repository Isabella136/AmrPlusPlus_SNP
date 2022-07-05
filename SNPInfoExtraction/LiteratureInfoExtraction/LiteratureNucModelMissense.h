#pragma once
#include "../ModelNucMissense.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureNucModelMissense : public ModelNucMissense, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureNucModelMissense(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureNucModelMissense();
    LiteratureNucModelMissense(const LiteratureNucModelMissense& other);
    Model* Clone();
    string condensedSNPinfo();
    int getPos();
};

LiteratureNucModelMissense::LiteratureNucModelMissense(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, true) {
    makeModel(line);
}
LiteratureNucModelMissense::LiteratureNucModelMissense(const LiteratureNucModelMissense& other) {
    this->wt_nuc = other.wt_nuc;
    this->pos = other.pos;
    this->mutant_nuc = other.mutant_nuc;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
}
Model* LiteratureNucModelMissense::Clone() {
    return new LiteratureNucModelMissense(*this);
}

int LiteratureNucModelMissense::getPos() {
    return pos;
}

LiteratureNucModelMissense::~LiteratureNucModelMissense() {}

string LiteratureNucModelMissense::condensedSNPinfo() {
    string toReturn = "Nuc:";
    toReturn += wt_nuc;
    toReturn += to_string(pos);
    for (int i = 0; i < mutant_nuc.size(); i++)
        toReturn += mutant_nuc[i];
    toReturn += addContext(pos - 2, pos, wt_nuc);
    return toReturn;
}
void LiteratureNucModelMissense::makeModel(string line) {
    wt_nuc = line.at(0);
    line = line.substr(1);
    int i = 0;
    while (isdigit(line.at(i))) {
        i += 1;
    }
    pos = stoi(line.substr(0, i));
    line = line.substr(i);
    for (int j = 0; j < line.length(); j++) {
        mutant_nuc.push_back(line.at(j));
    }

}