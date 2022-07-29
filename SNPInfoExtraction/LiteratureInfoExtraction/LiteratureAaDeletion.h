#pragma once
#include "../AaDeletion.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureAaDeletion : public AaDeletion, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureAaDeletion(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureAaDeletion();
    LiteratureAaDeletion(const LiteratureAaDeletion& other);
    InfoPipe* Clone();
    string condensedInfo();
    int getPos();
    string infoType();
};

LiteratureAaDeletion::LiteratureAaDeletion(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, false) {
    makeModel(line);
}
LiteratureAaDeletion::LiteratureAaDeletion(const LiteratureAaDeletion& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
}
InfoPipe* LiteratureAaDeletion::Clone() {
    return new LiteratureAaDeletion(*this);
}

int LiteratureAaDeletion::getPos() {
    return pos;
}

LiteratureAaDeletion::~LiteratureAaDeletion() {}

string LiteratureAaDeletion::condensedInfo() {
    string toReturn = "Del:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    toReturn += '-';
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
void LiteratureAaDeletion::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1, line.size() - 2));
}
string LiteratureAaDeletion::infoType() {
    return "Model";
}