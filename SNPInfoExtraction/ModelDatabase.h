#pragma once
#include <fstream>
#include <string>
#include <unordered_map>
#include <iostream>
#include "Model.h"

using namespace std;

class ModelDatabase{
    protected:
        unordered_map<string, list<Model*>> snpInfoDatabase;
    public:
        ModelDatabase();
        virtual void SNPInfo() = 0;
        unordered_map<string, list<Model*>> getDatabase();
};

ModelDatabase::ModelDatabase()
{

}
unordered_map<string, list<Model*>> ModelDatabase::getDatabase()
{
    return snpInfoDatabase;
}