#pragma once
#include <fstream>
#include <string>
#include <unordered_map>

using namespace std;

class CARD_database {
	private:
		unordered_map<string, string> databaseSequences;
	public:
		CARD_database();
		string getSequence(string id);
};
	
CARD_database::CARD_database()
{
	unordered_map<string, string> correctSequences;
	ifstream correctSearch;
	correctSearch.open("KargvaInfoExtracion/correctSequences.fasta");
	string line = "";
	string sequence = "";
	while (std::getline(correctSearch, line)) {
		std::getline(correctSearch, sequence);
		correctSequences.emplace(line.substr(1, line.find("|") - 1), sequence);
	}
	ifstream cardSearch;
	cardSearch.open("../CARD/card-data/protein_fasta_protein_variant_model.fasta");
	while (std::getline(cardSearch, line)) {
		std::getline(cardSearch, sequence);
		vector<string> header;
		header.push_back(line.substr(0, line.find("|")));
		while (line.find("|") != -1) {
			line = line.substr(line.find("|") + 1);
			header.push_back(line.substr(0, line.find("|")));
		}
		if (correctSequences.find(header[2]) != correctSequences.end())
			databaseSequences.emplace(header[2], correctSequences[header[2]]);
		else
			databaseSequences.emplace(header[2], sequence);
	}
	cardSearch.close();
}

string CARD_database::getSequence(string id) 
{
	if (databaseSequences.find(id) != databaseSequences.end())
		return databaseSequences[id];
	return "none";
}
