#pragma once
#include <fstream>
#include <string>
#include <unordered_map>
#include <iostream>
#include "GeneType.h"
#include "GeneClass.h"
#include "GeneMechanism.h"
#include "GeneGroup.h"
#include "Gene.h"

using namespace std;

class GeneDatabase{
    private:
        unordered_map<string, GeneType*>  geneTypes;
        unordered_map<string, GeneClass*>  geneClasses;
        unordered_map<string, GeneMechanism*>  geneMechanisms;
        unordered_map<string, GeneGroup*>  geneGroups;
        unordered_map<string, Gene*>  genes;
        unordered_map<string, pair<pair<char, int>, list<char>>> snpInfoDatabase;
    public:
        GeneDatabase(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
        void addGene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
        void makeAllOfInterest();
        void isItOfInterest(string category, fstream file, string ofInterest = "");
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
    
    //Updates geneTypes
    try {geneTypes.at(_geneType)->addGene(toAdd);}
	catch (const out_of_range& oor) {geneTypes.emplace(_geneType, toAdd);}

    //Updates geneClasses
    try {geneClasses.at(_geneClass)->addGene(toAdd);}
	catch (const out_of_range& oor) {geneClasses.emplace(_geneClass, toAdd);}

    //Updates geneMechanisms
    try {geneMechanisms.at(_geneMechanism)->addGene(toAdd);}
	catch (const out_of_range& oor) {geneMechanisms.emplace(_geneMechanism, toAdd);}

    //Updates geneGroups
    try {geneGroups.at(_geneGroup)->addGene(toAdd);}
	catch (const out_of_range& oor) {geneGroups.emplace(_geneGroup, toAdd);}

}
void GeneDatabase::isItOfInterest(string category, fstream file, string ofInterest = "")
{
    string header;
    string empty = "";
    string _geneName;
    string _geneType;
    string _geneClass;
    string _geneMechanism;
    string _geneGroup;
    string SNP;
    string restOfLine;
    std::getline(file, header);
    while (file.is_open())
    {
        std::getline(file, empty);
        std::getline(file, _geneName, ',');
        std::getline(file, _geneType, ',');
        std::getline(file, _geneClass, ',');
        std::getline(file, _geneMechanism, ',');
        std::getline(file, _geneGroup, ',');
        std::getline(file, SNP, ',');
        if (SNP != "RequiresSNPConfirmation")
            file.close();
        else {
            if (category == "type" && _geneType == ofInterest)
                geneTypes[ofInterest]->getGene(_geneName)->makeOfInterest();
            else if (category == "class" && _geneClass == ofInterest)
                geneClasses[ofInterest]->getGene(_geneName)->makeOfInterest();
            else if (category == "mechanism" && _geneMechanism == ofInterest)
                geneMechanisms[ofInterest]->getGene(_geneName)->makeOfInterest();
            else if (category == "group" && _geneGroup == ofInterest)
                geneGroups[ofInterest]->getGene(_geneName)->makeOfInterest();
            else if (category == "all")
                genes[_geneName]->makeOfInterest();
            else if (_geneName == ofInterest)
                genes[_geneName]->makeOfInterest();
        }
        std::getline(file, restOfLine);
    }
}
void GeneDatabase::makeAllOfInterest()
{
    for (auto iter = genes.begin(); iter != genes.end(); ++iter)
    {
        iter->second->makeOfInterest();
    }
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