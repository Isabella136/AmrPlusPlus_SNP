#include <fstream>
#include <iostream>
#include <string>
#include "GeneType.h"

using namespace std;

Gene getGene(string name);

int main(){
	fstream database;
	database.open("MEGARes_databse.fasta");
	string name;
	string sequence;
	while (getline(database, name))
	{
		getline(database, sequence);
		int lastDelimiter = name.find_last_of('|');
		if (name.substr(lastDelimiter + 1) == "RequiresSNPConfirmation")
		{
			//Make and store object
		}
	}
}

Gene getGene(string name)
{

}