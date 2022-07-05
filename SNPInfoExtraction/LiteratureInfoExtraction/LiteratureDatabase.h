#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "LiteratureModel.h"
#include "LiteratureModelDeletion.h"
#include "LiteratureModelMissense.h"
#include "LiteratureMultipleModels.h"
#include "LiteratureMultipleNucModels.h"
#include "LiteratureMustList.h"
#include "LiteratureMustModel.h"
#include "LiteratureNucModelDeletion.h"
#include "LiteratureNucModelMissense.h"
#include "LiteratureNucMustModel.h"
#include "../ModelDatabase.h"

using namespace std;

class LiteratureDatabase : public ModelDatabase {
	private:
		shared_ptr<MEGARes_database> databaseSequences;
	public:
		LiteratureDatabase();
		void SNPInfo();
};

LiteratureDatabase::LiteratureDatabase() {
	databaseSequences = make_shared<MEGARes_database>();
	SNPInfo();
}

void LiteratureDatabase::SNPInfo() {
	ifstream snpList;
	snpList.open("LiteratureInfoExtraction/SNPinfo_literature.csv");
	string line;
	while (std::getline(snpList, line)) {
		vector<string> headerAndSNP;
		headerAndSNP.push_back(line.substr(0, line.find(',')));
		LiteratureMustList* must1 = nullptr;
		LiteratureMustList* must2 = nullptr;
		while (line.find(',') != -1) {
			line = line.substr(line.find(',') + 1);
			headerAndSNP.push_back(line.substr(0, line.find(',')));
		}
		for (int i = 5; i < headerAndSNP.size(); i++) {
			if (headerAndSNP[i] == "RequiresSNPConfirmation") break;
			if (headerAndSNP[i] == "") continue;
			Model* model = nullptr;
			if (headerAndSNP[i].find("NucMult:") != -1)
				model = new LiteratureNucMultipleModels(headerAndSNP[i].substr(8), headerAndSNP[0], databaseSequences);
			else if (headerAndSNP[i].find("Mult:") != -1)
				model = new LiteratureMultipleModels(headerAndSNP[i].substr(5), headerAndSNP[0], databaseSequences);
			else if (headerAndSNP[i].find("MustNuc:") != -1)
				if (must1 == nullptr) {
					must1 = new LiteratureMustList(headerAndSNP[i], headerAndSNP[0], databaseSequences);
				}
				else {
					must1->addToList(headerAndSNP[i], headerAndSNP[0], databaseSequences);
				}
			else if (headerAndSNP[i].find("Must") != -1) {
				if (headerAndSNP[i].find("Must**:") != -1) {
					if (must2 == nullptr) {
						must2 = new LiteratureMustList(headerAndSNP[i].substr(7), headerAndSNP[0], databaseSequences);
					}
					else
						must2->addToList(headerAndSNP[i].substr(7), headerAndSNP[0], databaseSequences);
				}
				else {
					if (must1 == nullptr) {
						must1 = new LiteratureMustList(headerAndSNP[i].substr(5), headerAndSNP[0], databaseSequences);
					}
					else
						must1->addToList(headerAndSNP[i].substr(5), headerAndSNP[0], databaseSequences);
				}
			}
			else if (headerAndSNP[i].find("NucDel:") != -1)
				model = new LiteratureNucModelDeletion(headerAndSNP[i].substr(7), headerAndSNP[0], databaseSequences);
			else if (headerAndSNP[i].find("Del:") != -1)
				model = new LiteratureModelDeletion(headerAndSNP[i].substr(4), headerAndSNP[0], databaseSequences);
			else if (headerAndSNP[i].find("Nuc:") != -1)
				model = new LiteratureNucModelMissense(headerAndSNP[i].substr(4), headerAndSNP[0], databaseSequences);
			else
				model = new LiteratureModelMissense(headerAndSNP[i].substr(4), headerAndSNP[0], databaseSequences);
			if (model != nullptr) {
				if (snpInfoDatabase.find(headerAndSNP[0]) == snpInfoDatabase.end()) 
					snpInfoDatabase.emplace(headerAndSNP[0], list<Model*>());
				snpInfoDatabase.at(headerAndSNP[0]).push_back(model);
				
			}
		}
		if (must1 != nullptr) {
			if (snpInfoDatabase.find(headerAndSNP[0]) == snpInfoDatabase.end()) 
				snpInfoDatabase.emplace(headerAndSNP[0], list<Model*>());
			snpInfoDatabase.at(headerAndSNP[0]).push_back(must1);
			
		}
		if (must2 != nullptr) {
			if (snpInfoDatabase.find(headerAndSNP[0]) == snpInfoDatabase.end()) 
				snpInfoDatabase.emplace(headerAndSNP[0], list<Model*>());
			snpInfoDatabase.at(headerAndSNP[0]).push_back(must2);
			
		}
	}
	snpList.close();
}