#pragma once
#include <string>
#include "../Model.h"
#include "../CARD_database.h"

using namespace std;

class KargvaModel : public virtual Model {
    protected:
        static CARD_database databaseSequences;
        string cardID = "";
        string addContext(int lastLeftIndex, int firstRightIndex);
    public:
        KargvaModel(string id);
        virtual void addToModel(string line) = 0;
        virtual bool includes(string line) = 0;
        virtual string condensedSNPinfo() = 0;
};

KargvaModel::KargvaModel(string id){
    cardID = id;
}

string KargvaModel::addContext(int lastLeftIndex, int firstRightIndex) {
    string sequence = databaseSequences.getSequence(cardID);
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