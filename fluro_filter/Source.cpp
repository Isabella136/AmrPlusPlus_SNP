#include <fstream>
#include <iostream>
#include <string>

using namespace std;

void filterForFluro(string fileName);

int main() {
	string fileName;
	cout << "Specify file name or type x to exit: \n";
	cin >> fileName;
	while (fileName != "x")
	{
		filterForFluro(fileName);
		cout << "Specify file name or type x to exit: \n";
		cin >> fileName;
	}
	return 0;
}

void filterForFluro(string fileName)
{
	ifstream samFileByCategory;
	ifstream samFileByLine;
	samFileByCategory.open(fileName);
	samFileByLine.open(fileName);
	ofstream outputSamFile;
	outputSamFile.open("Fluro_" + fileName, ofstream::trunc);
	string QNAME;
	string FLAG;
	string RNAME;
	string toWrite;
	string restOfLine;
	while (std::getline(samFileByCategory, QNAME, '\t'))
	{
		std::getline(samFileByLine, toWrite);
		if (QNAME.at(0) != '@')
		{
			std::getline(samFileByCategory, FLAG, '\t');
			std::getline(samFileByCategory, RNAME, '\t');
			if (RNAME.find("Fluoroquinolones") != string::npos)
			{
				outputSamFile << toWrite << '\n';
			}
		}
		else
			outputSamFile << toWrite << '\n';
		std::getline(samFileByCategory, restOfLine);
	}
	samFileByCategory.close();
	samFileByLine.close();
	outputSamFile.close();
}