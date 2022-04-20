#pragma once
#include <string>
#include <list>

using namespace std;

class Model {
    protected:
        int pos;
        char wt_aa;
        list<char> mutant_aa;
        virtual void makeModel(string line) = 0;
    public:
        Model();
        virtual pair<pair<char, int>, list<char>> condensedSNPinfo() = 0;
};

Model::Model()
{
    pos = 0;
    wt_aa = 0;
}