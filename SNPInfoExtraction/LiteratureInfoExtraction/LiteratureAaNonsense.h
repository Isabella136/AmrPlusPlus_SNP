#pragma once
#include "../AaNonsense.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureAaNonsense : public AaNonsense, public LiteratureModel {
    private:
        void makeModel(string line);
    public:
        LiteratureAaNonsense(string line, string id, shared_ptr<MEGARes_database> dbSeq);
        ~LiteratureAaNonsense();
        LiteratureAaNonsense(const LiteratureAaNonsense& other);
        InfoPipe* Clone();
        string condensedInfo();
        string infoType();
        int getPos();
};

LiteratureAaNonsense::LiteratureAaNonsense(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, false) {
    makeModel(line);
}
LiteratureAaNonsense::~LiteratureAaNonsense() {}
LiteratureAaNonsense::LiteratureAaNonsense(const LiteratureAaNonsense& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
    this->nucleic = other.nucleic;
}
InfoPipe* LiteratureAaNonsense::Clone() {
    return new LiteratureAaNonsense(*this);
}
int LiteratureAaNonsense::getPos() {
    return pos;
}
void LiteratureAaNonsense::makeModel(string line) {
    wt_aa = line[0];
    pos = stoi(line.substr(1));
}
string LiteratureAaNonsense::condensedInfo() {
    string toReturn = "Non:";
    toReturn += wt_aa;
    toReturn += pos;
    toReturn += "*";
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
string LiteratureAaNonsense::infoType() {
    return "Model";
}