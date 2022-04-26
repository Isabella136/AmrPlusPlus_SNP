#pragma once
#include <string>
#include <list>
#include "../ModelDeletion.h"
#include "KargvaModel.h"

using namespace std;

class KargvaModel : public ModelDeletion, public KargvaModel {
private:
    void makeModel(string line);
    string addContext();
public:
    KargvaModel(string line);
    void addToModel(string line);
    bool includes(string line);
    string condensedSNPinfo();
};

KargvaModel::KargvaModel(string line) {
    makeModel(line);
}
void KargvaModel::addToModel(string line) {
    throw std::exception("should not have been called: model type is deletion");
}
bool KargvaModel::includes(string line) {
    return false;
}
string KargvaModel::condensedSNPinfo() {
    return make_pair(make_pair(wt_aa, pos), mutant_aa);
}
void KargvaModel::makeModel(string line) {
    wt_aa = line.at(0);
    pos = stoi(line.substr(1, line.size() - 2));
    mutant_aa.push_back(line.at(line.size() - 1));
}