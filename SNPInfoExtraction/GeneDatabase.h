#pragma once
#include <fstream>
#include <string>
#include <unordered_map>
#include <iostream>
#include "Gene.h"
#include "MetamarcInfoExtraction/MmarcDatabase.h"

using namespace std;

class GeneDatabase {
    private:
        unordered_map<string, Gene*>  genes;
        unordered_map<string, list<Model*>> snpInfoDatabase;
    public:
        GeneDatabase(MmarcDatabase& models);
        void combineDatabases(ModelDatabase& models);
        void print(string fileName);
        void addGene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
    
};
GeneDatabase::GeneDatabase(MmarcDatabase& models)
{
    snpInfoDatabase = models.getDatabase();
}
void GeneDatabase::combineDatabases(ModelDatabase& models)
{
    unordered_map<string, list<Model*>> temp = models.getDatabase();
    for (auto iter = temp.begin(); iter != temp.end(); ++iter) {
        if (snpInfoDatabase.find(iter->first) == snpInfoDatabase.end())
            snpInfoDatabase.emplace(iter->first, iter->second);
    }
}
void GeneDatabase::addGene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence)
{
    Gene* toAdd = new Gene(_geneName, _geneType, _geneClass, _geneMechanism, _geneGroup, _geneSequence);
    try
    {
        list<Model*> allSNPs = snpInfoDatabase.at(_geneName);
        for (auto iter = allSNPs.begin(); iter != allSNPs.end(); ++iter)
            toAdd->addSNP((*iter)->condensedSNPinfo().first, (*iter)->condensedSNPinfo().second);
    }
    catch (const out_of_range & oor) {}
    genes.emplace(_geneName, toAdd);
}
void GeneDatabase::print(string fileName)
{
    ofstream output;
    output.open(fileName);
    for (auto iter = genes.begin(); iter != genes.end(); ++iter)
        output << iter->second->getFASTA();
    output.close();
}