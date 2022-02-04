#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include "GeneType.h"

using namespace std;

unordered_map<string, GeneType*>  geneTypes;

Gene getGene(string name);

int main(){
	fstream database;
	database.open("../MEGARes_database.fasta");
	string name;
	string sequence;
	while (getline(database, name))
	{
		getline(database, sequence);
		int lastDelimiter = name.find_last_of('|');
		if (name.substr(lastDelimiter + 1) == "RequiresSNPConfirmation")
		{
			name = name.substr(1); //remove '>' character
			int delimiter1 = name.find('|');
			int delimiter2 = name.substr(delimiter1 + 1).find('|') + delimiter1 + 1;
			string type = name.substr(delimiter1 + 1, delimiter2 - delimiter1 - 1);
			try 
			{
				geneTypes.at(type)->addGene(name, sequence, delimiter2);
			}
			catch (const out_of_range& err)
			{
				geneTypes.emplace(type, new GeneType(name, sequence));
			}
		}
	}
	return 0;
}

Gene getGene(string name)
{
	int delimiter1 = name.find('|');
	int delimiter2 = name.substr(delimiter1 + 1).find('|') + delimiter1 + 1;
	string type = name.substr(delimiter1 + 1, delimiter2 - delimiter1 - 1);
	return *geneTypes.at(type)->getGene(name);
}