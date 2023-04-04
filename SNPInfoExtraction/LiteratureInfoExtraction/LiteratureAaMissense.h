#pragma once
#include "../AaMissense.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureAaMissense : public AaMissense, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureAaMissense(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureAaMissense();
    LiteratureAaMissense(const LiteratureAaMissense& other);
    InfoPipe* Clone();
    string condensedInfo();
    int getPos();
    string infoType();
};

LiteratureAaMissense::LiteratureAaMissense(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, false) {
    makeModel(line);
}
LiteratureAaMissense::LiteratureAaMissense(const LiteratureAaMissense& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->mutant_aa = other.mutant_aa;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
    this->nucleic = other.nucleic;
    this->source = other.source;
}
InfoPipe* LiteratureAaMissense::Clone() {
    return new LiteratureAaMissense(*this);
}

int LiteratureAaMissense::getPos() {
    return pos;
}

LiteratureAaMissense::~LiteratureAaMissense() {}

string LiteratureAaMissense::condensedInfo() {
    this->mutant_aa = sortMutant(this->mutant_aa);
    string toReturn = "Mis:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    for (int i = 0; i < mutant_aa.size(); i++)
        toReturn += mutant_aa[i];
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
void LiteratureAaMissense::makeModel(string line) {
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
    source = "Literature";
}
string LiteratureAaMissense::infoType() {
    return "Model";
}