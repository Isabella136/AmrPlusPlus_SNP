#pragma once
#include <string>
#include <list>

using namespace std;

class Model {
	public:
		Model();
		virtual string condensedSNPinfo() = 0;
};

Model::Model() {}
