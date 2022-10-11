#pragma once
#include <string>

using namespace std;

class InfoPipe {
	protected:
		string source = "CARD";
	public:
		InfoPipe();
		virtual ~InfoPipe();
		virtual InfoPipe* Clone() = 0;
		virtual string infoType() = 0;
		virtual string condensedInfo() = 0;
		string getSource();
};

InfoPipe::InfoPipe() {}
InfoPipe::~InfoPipe() {}
string InfoPipe::getSource() {
	return source;
}
