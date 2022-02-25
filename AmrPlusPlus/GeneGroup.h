#pragma once
#include <string>
#include <unordered_map>
#include <iostream>
#include "Gene.h"

using namespace std;

class GeneGroup
{
	private:
		string geneGroup;
		unordered_map<string,Gene*> genes;
	public:
		GeneGroup(string _geneGroup, Gene& gene);
		Gene* getGene(string geneName);
		void addGene(Gene& gene);
};

GeneGroup::GeneGroup(string _geneGroup, Gene& gene)
{
	geneGroup = _geneGroup;
	addGene(gene);
}

Gene* GeneGroup::getGene(string geneName)
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

void GeneGroup::addGene(Gene& gene)
{
	genes.emplace(gene.getName(), gene);
}