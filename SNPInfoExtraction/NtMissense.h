#pragma once
#include <string>
#include <vector>
#include "ModelSNP.h"

using namespace std;

class NtMissense : public virtual ModelSNP {
	protected:
		char wt_nuc = 0;
		vector<char> mutant_nuc;
		virtual void makeModel(string line) = 0;
		vector<char> sortMutant(vector<char> mtVector);
	public:
		NtMissense();
		virtual ~NtMissense();
		virtual string condensedInfo() = 0;
};

NtMissense::NtMissense() {}
NtMissense::~NtMissense() {}
vector<char> NtMissense::sortMutant(vector<char> mtVector) {
	if (mtVector.size() > 1) {
		vector<char> mtVectorA;
		vector<char> mtVectorB;
		for (int i = 0; i < mtVector.size(); i++) {
			if (i < mtVector.size() / 2)
				mtVectorA.push_back(mtVector[i]);
			else
				mtVectorB.push_back(mtVector[i]);
		}
		sortMutant(mtVectorA);
		sortMutant(mtVectorB);
		mtVector.clear();
		while (!mtVectorA.empty() || !mtVectorB.empty()) {
			if (mtVectorA.empty()) {
				mtVector.push_back(mtVectorB[0]);
				mtVectorB.erase(mtVectorB.begin());
			}
			else if (mtVectorB.empty()) {
				mtVector.push_back(mtVectorA[0]);
				mtVectorA.erase(mtVectorA.begin());
			}
			else if (mtVectorB[0] < mtVectorA[0]) {
				mtVector.push_back(mtVectorB[0]);
				mtVectorB.erase(mtVectorB.begin());
			}
			else {
				mtVector.push_back(mtVectorA[0]);
				mtVectorA.erase(mtVectorA.begin());
			}
		}
	}
	return mtVector;
}