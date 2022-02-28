#pragma once
#include <string>
#include <unordered_map>
#include "Gene.h"

using namespace std;

class GeneType
{
	private:
		string geneType;
		unordered_map<string,Gene*> genes;
	public:
		GeneType(string _geneType, Gene* gene);
		Gene* getGene(string geneName);
		void addGene(Gene* gene);
};

GeneType::GeneType(string _geneType, Gene* gene)
{
	geneType = _geneType;
	addGene(gene);
}

Gene* GeneType::getGene(string geneName)
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

void GeneType::addGene(Gene* gene)
{
	genes.emplace(gene->getName(), gene);
}