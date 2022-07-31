#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "LiteratureModel.h"
#include "LiteratureAaDeletion.h"
#include "LiteratureAaMissense.h"
#include "LiteratureAaMultiple.h"
#include "LiteratureAaNonsense.h"
#include "LiteratureNtMultiple.h"
#include "LiteratureIntrinsicList.h"
#include "LiteratureNtDeletion.h"
#include "LiteratureNtMissense.h"
#include "../FrameshiftInfo.h"
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
		LiteratureIntrinsicList* must = nullptr;
		while (line.find(',') != -1) {
			line = line.substr(line.find(',') + 1);
			headerAndSNP.push_back(line.substr(0, line.find(',')));
		}
		for (int i = 5; i < headerAndSNP.size(); i++) {
			if (headerAndSNP[i] == "RequiresSNPConfirmation") break;
			if (headerAndSNP[i] == "") continue;
			InfoPipe* info = nullptr;
			if (headerAndSNP[i].find("NucMult:") != -1)
				info = new LiteratureNtMultiple(headerAndSNP[i].substr(8), headerAndSNP[0], databaseSequences);
			else if (headerAndSNP[i].find("Mult:") != -1)
				info = new LiteratureAaMultiple(headerAndSNP[i].substr(5), headerAndSNP[0], databaseSequences);
			else if (headerAndSNP[i].find("MustNuc:") != -1)
				if (must == nullptr) {
					must = new LiteratureIntrinsicList(headerAndSNP[i], headerAndSNP[0], databaseSequences);
				}
				else {
					must->addToList(headerAndSNP[i], headerAndSNP[0], databaseSequences);
				}
			else if (headerAndSNP[i].find("Must:") != -1) {
				if (must == nullptr) {
					must = new LiteratureIntrinsicList(headerAndSNP[i].substr(5), headerAndSNP[0], databaseSequences);
				}
				else
					must->addToList(headerAndSNP[i].substr(5), headerAndSNP[0], databaseSequences);
			}
			else if (headerAndSNP[i].find("FS-") != -1)
				info = new FrameshiftInfo(headerAndSNP[i]);
			else if (headerAndSNP[i].find("Non:") != -1)
				info = new LiteratureAaNonsense(headerAndSNP[i].substr(4), headerAndSNP[0], databaseSequences);
			else if (headerAndSNP[i].find("NucDel:") != -1)
				info = new LiteratureNtDeletion(headerAndSNP[i].substr(7), headerAndSNP[0], databaseSequences);
			else if (headerAndSNP[i].find("Del:") != -1)
				info = new LiteratureAaDeletion(headerAndSNP[i].substr(4), headerAndSNP[0], databaseSequences);
			else if (headerAndSNP[i].find("Nuc:") != -1)
				info = new LiteratureNtMissense(headerAndSNP[i].substr(4), headerAndSNP[0], databaseSequences);
			else
				info = new LiteratureAaMissense(headerAndSNP[i].substr(4), headerAndSNP[0], databaseSequences);
			if (info != nullptr) {
				if (snpInfoDatabase.find(headerAndSNP[0]) == snpInfoDatabase.end()) 
					snpInfoDatabase.emplace(headerAndSNP[0], list<InfoPipe*>());
				snpInfoDatabase.at(headerAndSNP[0]).push_back(info);
				
			}
		}
		if (must != nullptr) {
			if (snpInfoDatabase.find(headerAndSNP[0]) == snpInfoDatabase.end()) 
				snpInfoDatabase.emplace(headerAndSNP[0], list<InfoPipe*>());
			snpInfoDatabase.at(headerAndSNP[0]).push_back(must);
		}
	}
	snpList.close();
}