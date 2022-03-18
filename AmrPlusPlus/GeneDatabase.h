#pragma once
#include <fstream>
#include <string>
#include <unordered_map>
#include <iostream>
#include "Gene.h"

using namespace std;

class GeneDatabase{
    private:
        unordered_map<string, Gene*>  genes;
        unordered_map<string, pair<pair<char, int>, list<char>>> snpInfoDatabase;
    public:
        GeneDatabase(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
        void addGene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
        void SNPInfo(fstream file);
        void print(string fileName);
};

GeneDatabase::GeneDatabase(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence)
{
    
    addGene(_geneName, _geneType, _geneClass, _geneMechanism, _geneGroup, _geneSequence);
}
void GeneDatabase::addGene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence)
{
    Gene* toAdd = new Gene(_geneName, _geneType, _geneClass, _geneMechanism, _geneGroup, _geneSequence);
    genes.emplace(_geneName, toAdd);

}
void GeneDatabase::SNPInfo(fstream file)
{
    
}
void GeneDatabase::print(string fileName)
{
    ofstream output;
    output.open(fileName);
    for (auto iter = genes.begin(); iter != genes.end(); ++iter)
    {
        output << iter->second->getFASTA();
    }
    output.close();
}