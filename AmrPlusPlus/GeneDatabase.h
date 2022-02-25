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
    public:
        GeneDatabase(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
        void isItOfInterest(string category, fstream file);
        void SNPInfo(fstream file);
        void print();
};

GeneDatabase::GeneDatabase(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence)
{

}
void GeneDatabase::isItOfInterest(string category, fstream file)
{

}
void GeneDatabase::SNPInfo(fstream file)
{

}
void GeneDatabase::print()
{
    
}