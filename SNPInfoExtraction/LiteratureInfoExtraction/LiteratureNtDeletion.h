#pragma once
#include "../NtDeletion.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureNtDeletion : public NtDeletion, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureNtDeletion(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureNtDeletion();
    LiteratureNtDeletion(const LiteratureNtDeletion& other);
    InfoPipe* Clone();
    string condensedInfo();
    int getPos();
    string infoType();
};

LiteratureNtDeletion::LiteratureNtDeletion(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, true) {
    makeModel(line);
}
LiteratureNtDeletion::LiteratureNtDeletion(const LiteratureNtDeletion& other) {
    this->wt_nuc = other.wt_nuc;
    this->pos = other.pos;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
}
InfoPipe* LiteratureNtDeletion::Clone() {
    return new LiteratureNtDeletion(*this);
}

int LiteratureNtDeletion::getPos() {
    return pos;
}

LiteratureNtDeletion::~LiteratureNtDeletion() {}

string LiteratureNtDeletion::condensedInfo() {
    string toReturn = "NucDel:";
    toReturn += wt_nuc;
    toReturn += to_string(pos);
    toReturn += '-';
    toReturn += addContext(pos - 2, pos, wt_nuc);
    return toReturn;
}
void LiteratureNtDeletion::makeModel(string line) {
    wt_nuc = line.at(0);
    pos = stoi(line.substr(1, line.size() - 2));
}
string LiteratureNtDeletion::infoType() {
    return "Model";
}
