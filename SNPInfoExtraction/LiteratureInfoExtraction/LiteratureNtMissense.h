#pragma once
#include "../NtMissense.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureNtMissense : public NtMissense, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureNtMissense(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureNtMissense();
    LiteratureNtMissense(const LiteratureNtMissense& other);
    InfoPipe* Clone();
    string condensedInfo();
    int getPos();
    string infoType();
};

LiteratureNtMissense::LiteratureNtMissense(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, true) {
    makeModel(line);
}
LiteratureNtMissense::LiteratureNtMissense(const LiteratureNtMissense& other) {
    this->wt_nuc = other.wt_nuc;
    this->pos = other.pos;
    this->mutant_nuc = other.mutant_nuc;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
    this->nucleic = other.nucleic;
    this->source = other.source;
}
InfoPipe* LiteratureNtMissense::Clone() {
    return new LiteratureNtMissense(*this);
}

int LiteratureNtMissense::getPos() {
    return pos;
}

LiteratureNtMissense::~LiteratureNtMissense() {}

string LiteratureNtMissense::condensedInfo() {
    string toReturn = "Nuc:";
    toReturn += wt_nuc;
    toReturn += to_string(pos);
    for (int i = 0; i < mutant_nuc.size(); i++)
        toReturn += mutant_nuc[i];
    toReturn += addContext(pos - 2, pos, wt_nuc);
    return toReturn;
}
void LiteratureNtMissense::makeModel(string line) {
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
    source = "Literature";
}
string LiteratureNtMissense::infoType() {
    return "Model";
}