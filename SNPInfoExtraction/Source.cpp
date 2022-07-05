/*
	AMRPlusPlus_SNP_Verification
	Copyright(C) 2022  Nathalie Bonin

	This program is free software : you can redistribute it and /or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.If not, see https ://www.gnu.org/licenses/.
*/

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include "MetamarcInfoExtraction/MmarcDatabase.h"
#include "KargvaInfoExtracion/KargvaDatabase.h"
#include "LiteratureInfoExtraction/LiteratureDatabase.h"
#include "GeneDatabase.h"

using namespace std;

int main()
{
	MmarcDatabase mmarc;
	KargvaDatabase kargva;
	LiteratureDatabase literature;
	GeneDatabase megaresSNP(mmarc);
	megaresSNP.combineDatabases(kargva);
	megaresSNP.combineDatabases(literature);
	//GeneDatabase megaresSNP(kargva);
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
	megaresSNP.print("../extracted_SNP_files/SNPinfo.csv", "../extracted_SNP_files/SNPinfo.fasta");
	return 0;
}