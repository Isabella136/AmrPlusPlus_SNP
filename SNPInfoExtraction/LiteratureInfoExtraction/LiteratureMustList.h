#pragma once
#include <string>
#include <map>
#include "LiteratureNucMustModel.h"
#include "LiteratureMustModel.h"

using namespace std;

class LiteratureMustList: public Model {
	private:
		map<int, LiteratureModel*> mustHaveList;
	public:
		LiteratureMustList();
		LiteratureMustList(string line, string id, shared_ptr<MEGARes_database> dbSeq);
		~LiteratureMustList();
		LiteratureMustList(const LiteratureMustList& other);
		void addToList(string line, string id, shared_ptr<MEGARes_database> dbSeq);
		LiteratureMustList* Clone();
		string condensedSNPinfo();
		
};

LiteratureMustList::LiteratureMustList() {}
LiteratureMustList::LiteratureMustList(string line, string id, shared_ptr<MEGARes_database> dbSeq) {
	addToList(line, id, dbSeq);
}
LiteratureMustList::~LiteratureMustList() {}
LiteratureMustList::LiteratureMustList(const LiteratureMustList& other) {
	this->mustHaveList = other.mustHaveList;
}
void LiteratureMustList::addToList(string line, string id, shared_ptr<MEGARes_database> dbSeq) {
	LiteratureModel* must;
	if (line.find("MustNuc:") != -1)
		must = new LiteratureNucMustModel(line.substr(8), id, dbSeq);
	else
		must = new LiteratureMustModel(line, id, dbSeq);
	mustHaveList.emplace(must->getPos(), must);
}
LiteratureMustList* LiteratureMustList::Clone() {
	return new LiteratureMustList(*this);
}

string LiteratureMustList::condensedSNPinfo() {
	string toReturn = "Must:";
	for (auto iter = mustHaveList.begin(); iter != mustHaveList.end(); ++iter) {
		toReturn += iter->second->condensedSNPinfo();
		toReturn += ";";
	}
	return toReturn;
}