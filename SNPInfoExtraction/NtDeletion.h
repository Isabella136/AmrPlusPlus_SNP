#pragma once
#include <string>
#include "ModelSingle.h"

using namespace std;

class NtDeletion : public virtual ModelSingle {
	protected:
		char wt_nuc = 0;
		virtual void makeModel(string line) = 0;
	public:
		NtDeletion();
		virtual ~NtDeletion();
		virtual string condensedInfo() = 0;
};

NtDeletion::NtDeletion() {}
NtDeletion::~NtDeletion() {}
