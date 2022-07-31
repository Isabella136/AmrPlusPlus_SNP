#pragma once
#include <fstream>
#include <string>
#include <unordered_map>
#include <map>
#include <iostream>
#include "Gene.h"
#include "MetamarcInfoExtraction/MmarcDatabase.h"

using namespace std;

class GeneDatabase {
    private:
        map<int, Gene*>  genes;
        unordered_map<string, list<InfoPipe*>> snpInfoDatabase;
    public:
        GeneDatabase(ModelDatabase& models);
        void combineDatabases(ModelDatabase& models);
        void print(string csv, string fasta);
        void printIfNoSNP(string fileName);
        void addGene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
    
};
GeneDatabase::GeneDatabase(ModelDatabase& models)
{
    snpInfoDatabase = models.getDatabase();
}
void GeneDatabase::combineDatabases(ModelDatabase& models)
{
    unordered_map<string, list<InfoPipe*>> temp = models.getDatabase();
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
        list<InfoPipe*> allSNPs = snpInfoDatabase.at(_geneName);
        for (auto iter = allSNPs.begin(); iter != allSNPs.end(); ++iter)
            toAdd->addSNP((*iter)->condensedInfo());
    }
    catch (const out_of_range & oor) {}
    int megNumber = stoi(_geneName.substr(_geneName.find('_') + 1));
    genes.emplace(megNumber, toAdd);
}
void GeneDatabase::print(string csv, string fasta)
{
    ofstream outputCsv;
    outputCsv.open(csv);
    ofstream outputFasta;
    outputFasta.open(fasta);
    for (auto iter = genes.begin(); iter != genes.end(); ++iter) {
        outputCsv << iter->second->getHeader();
        outputFasta << iter->second->getFASTA();
    }
    outputCsv.close();
}
void GeneDatabase::printIfNoSNP(string fileName)
{
    ofstream output;
    output.open(fileName);
    for (auto iter = genes.begin(); iter != genes.end(); ++iter)
        output << iter->second->getFASTA2();
    output.close();
}