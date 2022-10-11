#pragma once
#include <string>
#include <map>
#include "LiteratureNtIntrinsic.h"
#include "LiteratureAaIntrinsic.h"

using namespace std;

class LiteratureIntrinsicList : public Model {
	private:
		map<int, LiteratureModel*> mustHaveList;
	public:
		LiteratureIntrinsicList();
		LiteratureIntrinsicList(string line, string id, shared_ptr<MEGARes_database> dbSeq);
		~LiteratureIntrinsicList();
		LiteratureIntrinsicList(const LiteratureIntrinsicList& other);
		void addToList(string line, string id, shared_ptr<MEGARes_database> dbSeq);
		InfoPipe* Clone();
		list<int> getFirstPos();
		string condensedInfo();
		string infoType();
};
list<int> LiteratureIntrinsicList::getFirstPos() {
	return {0, *(++(mustHaveList.begin()->second->getFirstPos()).begin())};	//intrinsic come before snp(1);
}
LiteratureIntrinsicList::LiteratureIntrinsicList() {}
LiteratureIntrinsicList::LiteratureIntrinsicList(string line, string id, shared_ptr<MEGARes_database> dbSeq) {
	addToList(line, id, dbSeq);
}
LiteratureIntrinsicList::~LiteratureIntrinsicList() {}
LiteratureIntrinsicList::LiteratureIntrinsicList(const LiteratureIntrinsicList& other) {
	for (auto iter = other.mustHaveList.begin(); iter != other.mustHaveList.end(); ++iter) {
		this->mustHaveList.emplace(iter->first, dynamic_cast<LiteratureModel*>(iter->second->Clone()));
	}
}
void LiteratureIntrinsicList::addToList(string line, string id, shared_ptr<MEGARes_database> dbSeq) {
	LiteratureModel* must;
	if (line.find("MustNuc:") != -1)
		must = new LiteratureNtIntrinsic(line.substr(8), id, dbSeq);
	else
		must = new LiteratureAaIntrinsic(line, id, dbSeq);
	mustHaveList.emplace(must->getPos(), must);
}
InfoPipe* LiteratureIntrinsicList::Clone() {
	return new LiteratureIntrinsicList(*this);
}

string LiteratureIntrinsicList::condensedInfo() {
	string toReturn = "Must:";
	for (auto iter = mustHaveList.begin(); iter != mustHaveList.end(); ++iter) {
		toReturn += iter->second->condensedInfo();
		toReturn += ";";
	}
	return toReturn.substr(0,toReturn.length()-1);
}
string LiteratureIntrinsicList::infoType() {
	return "List";
}
