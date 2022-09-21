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
	public:
		NtMissense();
		virtual ~NtMissense();
		virtual string condensedInfo() = 0;
};

NtMissense::NtMissense() {}
NtMissense::~NtMissense() {}