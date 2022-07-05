#pragma once
#include <string>
#include "../Model.h"
#include "../MEGARes_database.h"

using namespace std;

class LiteratureModel : public virtual Model {
protected:
    shared_ptr<MEGARes_database> databaseSequences;
    bool nucleic = false;
    string megID = "";
    string addContext(int lastLeftIndex, int firstRightIndex, char wt);
public:
    LiteratureModel();
    LiteratureModel(string id, shared_ptr<MEGARes_database> dbSeq, bool nuc);
    ~LiteratureModel();
    virtual int getPos() = 0;
    virtual string condensedSNPinfo() = 0;
};
LiteratureModel::LiteratureModel() {}

LiteratureModel::LiteratureModel(string id, shared_ptr<MEGARes_database> dbSeq, bool nuc) {
    megID = id;
    databaseSequences = move(dbSeq);
    nucleic = nuc;
}

LiteratureModel::~LiteratureModel() {}

string LiteratureModel::addContext(int lastLeftIndex, int firstRightIndex, char wt) {
    string sequence = "";
    if (!nucleic)
        sequence = databaseSequences->getAASequence(megID);
    else
        sequence = databaseSequences->getSequence(megID);
    if (sequence == "none")
        return "__";
    string toReturn = "_";
    int firstLeftIndex = lastLeftIndex - 4;
    if (lastLeftIndex < 4)
        firstLeftIndex = 0;
    int leftContextSize = lastLeftIndex - firstLeftIndex + 1;
    if (leftContextSize != 0)
        toReturn += sequence.substr(firstLeftIndex, leftContextSize);
    toReturn += "_";
    int lastRightIndex = firstRightIndex + 4;
    if (lastRightIndex > sequence.size() - 1)
        lastRightIndex = sequence.size() - 1;
    int rightContextSize = lastRightIndex - firstRightIndex + 1;
    if (rightContextSize != 0)
        toReturn += sequence.substr(firstRightIndex, rightContextSize);
    return toReturn;
}