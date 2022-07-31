#pragma once
#include <string>
#include "../AaIntrinsic.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureAaIntrinsic : public AaIntrinsic, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureAaIntrinsic(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureAaIntrinsic();
    LiteratureAaIntrinsic(const LiteratureAaIntrinsic& other);
    InfoPipe* Clone();
    string condensedInfo();
    int getPos();
    string infoType();
};

LiteratureAaIntrinsic::LiteratureAaIntrinsic(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, false) {
    makeModel(line);
}
LiteratureAaIntrinsic::LiteratureAaIntrinsic(const LiteratureAaIntrinsic& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
    this->nucleic = other.nucleic;
}
InfoPipe* LiteratureAaIntrinsic::Clone() {
    return new LiteratureAaIntrinsic(*this);
}

int LiteratureAaIntrinsic::getPos() {
    return pos;
}

LiteratureAaIntrinsic::~LiteratureAaIntrinsic() {}

string LiteratureAaIntrinsic::condensedInfo() {
    string toReturn = "Amino:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    toReturn += addContext(pos - 2, pos, wt_aa);
    return toReturn;
}
void LiteratureAaIntrinsic::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1));
}
string LiteratureAaIntrinsic::infoType() {
    return "Model";
}