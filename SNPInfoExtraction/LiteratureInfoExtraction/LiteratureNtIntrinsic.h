#pragma once
#include <string>
#include "../NtIntrinsic.h"
#include "LiteratureModel.h"

using namespace std;

class LiteratureNtIntrinsic : public NtIntrinsic, public LiteratureModel {
private:
    void makeModel(string line);
public:
    LiteratureNtIntrinsic(string line, string id, shared_ptr<MEGARes_database> dbSeq);
    ~LiteratureNtIntrinsic();
    LiteratureNtIntrinsic(const LiteratureNtIntrinsic& other);
    InfoPipe* Clone();
    string condensedInfo();
    int getPos();
    string infoType();
};

LiteratureNtIntrinsic::LiteratureNtIntrinsic(string line, string id, shared_ptr<MEGARes_database> dbSeq) :LiteratureModel(id, dbSeq, true) {
    makeModel(line);
}
LiteratureNtIntrinsic::LiteratureNtIntrinsic(const LiteratureNtIntrinsic& other) {
    this->wt_nuc = other.wt_nuc;
    this->pos = other.pos;
    this->megID = other.megID;
    this->databaseSequences = other.databaseSequences;
    this->nucleic = other.nucleic;
    this->source = other.source;
}
InfoPipe* LiteratureNtIntrinsic::Clone() {
    return new LiteratureNtIntrinsic(*this);
}

int LiteratureNtIntrinsic::getPos() {
    return pos;
}

LiteratureNtIntrinsic::~LiteratureNtIntrinsic() {}

string LiteratureNtIntrinsic::condensedInfo() {
    string toReturn = "Nuc:";
    toReturn += wt_nuc;
    toReturn += to_string(pos);
    toReturn += addContext(pos - 2, pos, wt_nuc);
    return toReturn;
}
void LiteratureNtIntrinsic::makeModel(string line) {
    wt_nuc = line.at(0);
    pos = stoi(line.substr(1));
    source = "Literature";
}
string LiteratureNtIntrinsic::infoType() {
    return "Model";
}