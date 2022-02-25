#pragma once
#include <string>
#include <unordered_map>
#include "Gene.h"

using namespace std;

class GeneClass
{
	private:
		string geneClass;
		unordered_map<string,Gene*> genes;
	public:
		GeneClass(string _geneClass, Gene& gene);
		Gene* getGene(string geneName);
		void addGene(Gene& gene);
};

GeneClass::GeneClass(string _geneClass, Gene& gene)
{
	geneClass = _geneClass;
	addGene(gene);
}

Gene* GeneClass::getGene(string geneName)
{
	try 
	{
		return genes.at(geneName);
	}
	catch (const out_of_range& oor)
	{
		return nullptr;
	}
}

void GeneClass::addGene(Gene& gene)
{
	genes.emplace(gene.getName(), gene);
}