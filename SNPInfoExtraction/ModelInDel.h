#pragma once
#include <string>
#include <list>
#include "Model.h"

using namespace std;

class ModelInDel : public virtual Model {
protected:
	list<int> pos;
	virtual void makeModel(string line) = 0;
	virtual ~ModelInDel();
public:
	ModelInDel();
	list<int> getFirstPos();
	virtual string condensedInfo() = 0;
};

ModelInDel::ModelInDel() {}
ModelInDel::~ModelInDel() {}
list<int> ModelInDel::getFirstPos() {
	return {3, *(pos.begin())}; //deletion come after hypersusceptible(2) and before insertion(4)
}
