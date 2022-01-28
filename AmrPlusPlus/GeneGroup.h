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
		unordered_map<string, Gene*> genes;
	public:
		GeneGroup(string geneName, string geneSequence, int delimiter4, int delimiter5);
		Gene* getGene(string geneName);
		void addGene(string geneName, string geneSequence);
};

GeneGroup::GeneGroup(string geneName, string geneSequence, int delimiter4, int delimiter5)
{
	geneGroup = geneName.substr(delimiter4 + 1, delimiter5 - delimiter4 - 1);
	addGene(geneName, geneSequence);
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
void GeneGroup::addGene(string geneName, string geneSequence)
{
	if (!getGene(geneName))
	{
		genes.emplace(geneGroup, new Gene(geneName, geneSequence));
	}
	else
	{
		cout << "Gene name \"" << geneName << "\" already exists";
	}
}