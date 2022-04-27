#pragma once
#include "../ModelInsertion.h"
#include "KargvaModel.h"

using namespace std;

class KargvaModelInsertion : public ModelInsertion, public KargvaModel {
private:
    void makeModel(string line);
public:
    KargvaModelInsertion(string line, string id);
    void addToModel(string line);
    bool includes(string line);
    string condensedSNPinfo();
};

KargvaModelInsertion::KargvaModelInsertion(string line, string id) :KargvaModel(id) {
    makeModel(line);
}
void KargvaModelInsertion::addToModel(string line) {
    mutant_aa.push_back(line.at(line.size() - 1));
}
bool KargvaModelInsertion::includes(string line) {
    string temp = "-";
    temp += pos;
    return temp == line.substr(0, line.size() - 1);
}
string KargvaModelInsertion::condensedSNPinfo() {
    string toReturn = "Ins:";
    toReturn += '-';
    toReturn += pos;
    for (int i = 0; i < mutant_aa.size; i++)
        toReturn += mutant_aa[i];
    toReturn += addContext(pos - 2, pos - 1);
    return toReturn;
}
void KargvaModelInsertion::makeModel(string line) {
    pos = stoi(line.substr(1, line.size() - 2));
    mutant_aa.push_back(line.at(line.size() - 1));
}