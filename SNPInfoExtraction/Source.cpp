#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include "MetamarcInfoExtraction/MmarcDatabase.h"
#include "KargvaInfoExtracion/KargvaDatabase.h"
#include "GeneDatabase.h"

using namespace std;

int main()
{
	MmarcDatabase mmarc;
	KargvaDatabase kargva;
	GeneDatabase megaresSNP(mmarc);
	megaresSNP.combineDatabases(kargva);
	fstream database;
	database.open("../MEGARes_database.fasta");
	string header;
	string sequence;
	while (std::getline(database, header))
	{
		std::getline(database, sequence);
		int lastDelimiter = header.find_last_of('|');
		if (header.substr(lastDelimiter + 1) == "RequiresSNPConfirmation")
		{
			string _geneName = header.substr(1, header.find('|')-1);
			header = header.substr(header.find('|') + 1);
			string _geneType = header.substr(0, header.find('|'));
			header = header.substr(header.find('|') + 1);
			string _geneClass = header.substr(0, header.find('|'));
			header = header.substr(header.find('|') + 1);
			string _geneMechanism = header.substr(0, header.find('|'));
			header = header.substr(header.find('|') + 1);
			string _geneGroup = header.substr(0, header.find('|'));
			megaresSNP.addGene(_geneName, _geneType, _geneClass, _geneMechanism, _geneGroup, sequence);
		}
	}
	database.close();
	megaresSNP.print("../extracted_SNP_files/metamarcSNPinfo.fasta");
	return 0;
}