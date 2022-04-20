#pragma once
#include <string>
#include <list>
#include "../Model.h"

using namespace std;

class KargvaModel : public Model {
    private:
        void makeModel(string line);
    public:
        KargvaModel(string line);
        void addToModel(string line);
        bool includes(string line);
        pair<pair<char, int>, list<char>> condensedSNPinfo();
};

KargvaModel::KargvaModel(string line){
    makeModel(line);
}
void KargvaModel::addToModel(string line){
    mutant_aa.push_back(line.at(line.size()-1));
}
bool KargvaModel::includes(string line){
    string temp = "";
    temp += wt_aa;
    temp += pos;
    return temp == line.substr(0, line.size()-1);
}
pair<pair<char, int>, list<char>> KargvaModel::condensedSNPinfo(){
    return make_pair(make_pair(wt_aa, pos), mutant_aa);
}
void KargvaModel::makeModel(string line){
    wt_aa = line.at(0);
    pos = stoi(line.substr(1, line.size()-2));
    mutant_aa.push_back(line.at(line.size()-1));
}