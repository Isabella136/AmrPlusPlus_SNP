#pragma once
#include <fstream>
#include <string>
#include <unordered_map>
#include <algorithm>  

using namespace std;

class MEGARes_database {
	private:
		unordered_map<string, string> databaseSequences;
		string translate(string sequence);
	public:
		MEGARes_database();
		MEGARes_database(const MEGARes_database& other);
		string getSequence(string id);
		string getAASequence(string id);
};

MEGARes_database::MEGARes_database(const MEGARes_database& other) {
	this->databaseSequences = other.databaseSequences;
}
MEGARes_database::MEGARes_database()
{
	ifstream megaresSearch;
	megaresSearch.open("../MEGARes_database.fasta");
	string line = "";
	string sequence = "";
	while (std::getline(megaresSearch, line)) {
		std::getline(megaresSearch, sequence);
		databaseSequences.emplace(line.substr(1, line.find("|") - 1), sequence);
	}
	megaresSearch.close();
}

string MEGARes_database::getSequence(string id)
{
	if (databaseSequences.find(id) != databaseSequences.end())
		return databaseSequences[id];
	return "none";
}

string MEGARes_database::getAASequence(string id)
{
	if (databaseSequences.find(id) != databaseSequences.end())
		transform(databaseSequences[id].begin(), databaseSequences[id].end(), databaseSequences[id].begin(), ::toupper);
		return translate(databaseSequences[id]);
	return "none";
}

string MEGARes_database::translate(string sequence)
{
	string translated = "";
	while (sequence.length() > 0) {
		string codon = sequence.substr(0, 3);
		sequence = sequence.substr(3);
		if (codon[0] == 'T') {
			if (codon[1] == 'T') {
				if (codon[2] == 'T' || codon[2] == 'C') {
					translated += "F";
				}
				else {
					translated += "L";
				}
			}
			else if (codon[1] == 'C') {
				translated += "S";
			}
			else if (codon[1] == 'A') {
				if (codon[2] == 'T' || codon[2] == 'C') {
					translated += "Y";
				}
				else {
					translated += "*";
				}
			}
			else {
				if (codon[2] == 'T' || codon[2] == 'C') {
					translated += "C";
				}
				else if (codon[2] == 'A') {
					translated += "*";
				}
				else {
					translated += "W";
				}
			}
		}
		else if (codon[0] == 'C') {
			if (codon[1] == 'T') {
				translated += "L";
			}
			else if (codon[1] == 'C') {
				translated += "P";
			}
			else if (codon[1] == 'A') {
				if (codon[2] == 'T' || codon[2] == 'C') {
					translated += "H";
				}
				else {
					translated += "Q";
				}
			}
			else {
				translated += "R";
			}
		}
		else if (codon[0] == 'A') {
			if (codon[1] == 'T') {
				if (codon[2] == 'G') {
					translated += "M";
				}
				else {
					translated += "I";
				}
			}
			else if (codon[1] == 'C') {
				translated += "T";
			}
			else if (codon[1] == 'A') {
				if (codon[2] == 'T' || codon[2] == 'C') {
					translated += "N";
				}
				else {
					translated += "K";
				}
			}
			else {
				if (codon[2] == 'T' || codon[2] == 'C') {
					translated += "S";
				}
				else {
					translated += "R";
				}
			}
		}
		else if (codon[0] == 'G') {
			if (codon[1] == 'T') {
				translated += "V";
			}
			else if (codon[1] == 'C') {
				translated += "A";
			}
			else if (codon[1] == 'A') {
				if (codon[2] == 'T' || codon[2] == 'C') {
					translated += "D";
				}
				else {
					translated += "E";
				}
			}
			else {
				translated += "G";
			}
		}
	}
	return translated;
}