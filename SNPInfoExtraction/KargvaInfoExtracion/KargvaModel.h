#pragma once
#include <string>
#include <list>
#include "../Model.h"

using namespace std;

class KargvaModel : public Model {
    protected:
        string addContext(int pos);
    public:
        KargvaModel();
        virtual void addToModel(string line) = 0;
        virtual bool includes(string line) = 0;
        virtual string condensedSNPinfo() = 0;
};

KargvaModel::KargvaModel(){}

list<char> KargvaModel::addContext(int pos){
    return
}