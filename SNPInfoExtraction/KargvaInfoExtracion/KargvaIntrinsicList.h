#pragma once
#include <map>
#include "KargvaAaIntrinsic.h"

using namespace std;

class KargvaIntrinsicList : public Model {
	private:
		map<int, KargvaAaIntrinsic*> mustHaveList;
	public:
		KargvaIntrinsicList(string line, string id, shared_ptr<CARD_database> dbSeq);
		~KargvaIntrinsicList();
		KargvaIntrinsicList(const KargvaIntrinsicList& other);
		void addToList(string line, string id, shared_ptr<CARD_database> dbSeq);
		InfoPipe* Clone();
		string condensedInfo();
		string infoType();
};

KargvaIntrinsicList::KargvaIntrinsicList(string line, string id, shared_ptr<CARD_database> dbSeq) {
	addToList(line, id, dbSeq);
}
KargvaIntrinsicList::~KargvaIntrinsicList() {}
KargvaIntrinsicList::KargvaIntrinsicList(const KargvaIntrinsicList& other) {
	for (auto iter = other.mustHaveList.begin(); iter != other.mustHaveList.end(); ++iter) {
		this->mustHaveList.emplace(iter->first, iter->second->Clone());
	}
}
void KargvaIntrinsicList::addToList(string line, string id, shared_ptr<CARD_database> dbSeq) {
	KargvaAaIntrinsic* must = new KargvaAaIntrinsic(line.substr(5), id, dbSeq);
	mustHaveList.emplace(must->getPos(), must);
}
InfoPipe* KargvaIntrinsicList::Clone() {
	return new KargvaIntrinsicList(*this);
}
string KargvaIntrinsicList::condensedInfo() {
	string toReturn = "Must:";
	for (auto iter = mustHaveList.begin(); iter != mustHaveList.end(); ++iter) {
		toReturn += iter->second->condensedInfo();
		toReturn += ";";
	}
	return toReturn.substr(0, toReturn.length() - 1);
}
string KargvaIntrinsicList::infoType() {
	return "Model";
}
