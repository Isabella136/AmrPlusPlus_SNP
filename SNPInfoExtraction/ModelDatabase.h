#pragma once
#include <fstream>
#include <string>
#include <list>
#include <unordered_map>
#include <iostream>
#include "Model.h"

using namespace std;

class ModelDatabase{
    protected:
        unordered_map<string, list<InfoPipe*>> snpInfoDatabase;
    public:
        ModelDatabase();
        virtual void SNPInfo() = 0;
        unordered_map<string, list<InfoPipe*>> getDatabase();
};

ModelDatabase::ModelDatabase() {}
unordered_map<string, list<InfoPipe*>> ModelDatabase::getDatabase()
{
    return snpInfoDatabase;
}