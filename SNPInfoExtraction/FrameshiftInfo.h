#pragma once
#include <string>
#include "InfoPipe.h"

using namespace std;

class FrameshiftInfo : public InfoPipe {
	protected:
		string detail;
	public:
		FrameshiftInfo(string input);
		~FrameshiftInfo();
		FrameshiftInfo(const FrameshiftInfo& other);
		InfoPipe* Clone();
		string condensedInfo();
		string infoType();
};

FrameshiftInfo::FrameshiftInfo(string input) {
	detail = input;
}
FrameshiftInfo::FrameshiftInfo(const FrameshiftInfo& other) {
	this->detail = other.detail;
}
InfoPipe* FrameshiftInfo::Clone() {
	return new FrameshiftInfo(*this);
}
FrameshiftInfo::~FrameshiftInfo() {}

string FrameshiftInfo::infoType() {
	return "Frameshift";
}
string FrameshiftInfo::condensedInfo() {
	return detail;
}