#pragma once
#include <string>
#include <vector>
#include "ModelSingle.h"

using namespace std;

class ModelNucMissense : public virtual ModelSingle {
	protected:
		char wt_nuc = 0;
		vector<char> mutant_nuc;
		virtual void makeModel(string line) = 0;
	public:
		ModelNucMissense();
		virtual ~ModelNucMissense();
		virtual string condensedSNPinfo() = 0;
};

ModelNucMissense::ModelNucMissense() {}
ModelNucMissense::~ModelNucMissense() {}