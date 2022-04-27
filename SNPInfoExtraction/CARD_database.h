#pragma once
#include <fstream>
#include <string>
#include <unordered_map>

using namespace std;

struct hash_pair {
	template <class T1, class T2>
	size_t operator()(const pair<T1, T2>& p) const
	{
		auto hash1 = hash<T1>{}(p.first);
		auto hash2 = hash<T2>{}(p.second);
		return hash1 ^ hash2;
	}
};

using namespace std;
class CARD_database {
	private:
		unordered_map<string, string, hash_pair> databaseSequences;
	public:
		CARD_database();
		string getSequence(string id);
};
	
CARD_database::CARD_database()
{
	ifstream cardSearch;
	cardSearch.open("../CARD/card-data/protein_fasta_protein_variant_model.fasta");
	string line = "";
	string sequence = "";
	while (std::getline(cardSearch, line)) {
		std::getline(cardSearch, sequence);
		vector<string> header;
		header.push_back(line.substr(0, line.find("|")));
		while (line.find("|") != -1) {
			line = line.substr(line.find("|") + 1);
			header.push_back(line.substr(0, line.find("|")));
		}
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
