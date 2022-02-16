#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include "GeneType.h"

using namespace std;

unordered_map<string, GeneType*>  geneTypes;

Gene getGene(string name);
void makeDatabase();
void filterSAMfile(string fileName);

int main(){
	//makeDatabase();
	string fileName;
	cout << "Specify file name or type x to exit: \n";
	cin >> fileName;
	while (fileName != "x")
	{
		filterSAMfile(fileName);
		cout << "Specify file name or type x to exit: \n";
		cin >> fileName;
	}
	return 0;
}

void filterSAMfile(string fileName)
{
	ifstream samFileByCategory;
	ifstream samFileByLine;
	samFileByCategory.open(fileName);
	samFileByLine.open(fileName);
	ofstream outputSamFile;
	outputSamFile.open("Filtered_out_" + fileName, ofstream::trunc);
	string QNAME;
	string FLAG;
	string RNAME;
	string toWrite;
	string restOfLine;
	while (std::getline(samFileByCategory, QNAME, '\t'))
	{
		std::getline(samFileByLine, toWrite);
		if (QNAME.at(0) == '@')
		{
			outputSamFile << toWrite << '\n';
		}
		else
		{
			std::getline(samFileByCategory, FLAG, '\t');
			std::getline(samFileByCategory, RNAME, '\t');
			if (RNAME.find("RequiresSNPConfirmation") != string::npos)
			{
				outputSamFile << toWrite << '\n';
			}
		}
		std::getline(samFileByCategory, restOfLine);
	}
	samFileByCategory.close();
	samFileByLine.close();
	outputSamFile.close();
}

void makeDatabase()
{
	fstream database;
	database.open("../MEGARes_database.fasta");
	string name;
	string sequence;
	while (std::getline(database, name))
	{
		std::getline(database, sequence);
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
	database.close();
}

Gene getGene(string name)
{
	int delimiter1 = name.find('|');
	int delimiter2 = name.substr(delimiter1 + 1).find('|') + delimiter1 + 1;
	string type = name.substr(delimiter1 + 1, delimiter2 - delimiter1 - 1);
	return *geneTypes.at(type)->getGene(name);
}