#pragma once
#include <string>
#include "ModelSingle.h"

using namespace std;

class ModelNucDeletion: public virtual ModelSingle {
	protected:
		char wt_nuc = 0;
		virtual void makeModel(string line) = 0;
	public:
		ModelNucDeletion();
		virtual ~ModelNucDeletion();
		virtual string condensedSNPinfo() = 0;
};

ModelNucDeletion::ModelNucDeletion() {}
ModelNucDeletion::~ModelNucDeletion() {}
