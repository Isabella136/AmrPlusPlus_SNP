#pragma once
#include "InfoPipe.h"
#include "ModelSNP.h"
#include "ModelInDel.h"
#include "AaInsertion.h"
#include "LiteratureInfoExtraction/LiteratureNtMultiple.h"
#include "LiteratureInfoExtraction/LiteratureAaMultiple.h"
#include "KargvaInfoExtracion/KargvaAaMultiple.h"
#include <list>

using namespace std;


class Node {
private:
	Node* parent;
	Node* left;
	Node* right;
	InfoPipe* node;
	void findPos();
	list<int> pos;
public:
	Node(InfoPipe* node);
	list<int> getPos();
	Node* getChild(bool left);
	InfoPipe* getContent();
	void addChild(bool left, Node* child);
	void addParent(Node* parent);
};

class NodeTree {
private:
	Node* top;
	list<InfoPipe*> traverse(Node* current);
public:
	NodeTree(InfoPipe* top);
	void addChild(InfoPipe* child);
	list<InfoPipe*> returnOrderedNodes();
};

Node::Node(InfoPipe* node) {
	this->node = node;
	this->parent = nullptr;
	this->left = nullptr;
	this->right = nullptr;
	findPos();
}

void Node::addChild(bool left, Node* child) {
	if (left)
		this->left = child;
	else
		this->right = child;
}
void Node::addParent(Node* parent) {
	this->parent = parent;
}

void Node::findPos() {
	Model* temp = dynamic_cast<Model*>(node);
	pos = temp->getFirstPos();
}

list<int> Node::getPos() {
	return pos;
}

InfoPipe* Node::getContent() {
	return node;
}

Node* Node::getChild(bool left) {
	if (left)
		return this->left;
	return this->right;
}

NodeTree::NodeTree(InfoPipe* top) {
	this->top = new Node(top);
}

void NodeTree::addChild(InfoPipe* child) {
	Node* toAdd = new Node(child);
	bool parentFound = false;
	bool left = false;
	Node* currentNode = top;
	while (!parentFound) {
		list<int> toAddPosList = toAdd->getPos();
		list<int> currentNodePosList = currentNode->getPos();
		auto toAddPos = toAddPosList.begin();
		auto currentNodePos = currentNodePosList.begin();
		if (*toAddPos < *currentNodePos)
			left = true;
		else if (*toAddPos == *currentNodePos) {
			++toAddPos;
			++currentNodePos;
			if (*toAddPos < *currentNodePos)
				left = true;
			else
				left = false;
		}
		else
			left = false;
		if (currentNode->getChild(left) == nullptr)
			parentFound = true;
		else
			currentNode = currentNode->getChild(left);
	}
	currentNode->addChild(left, toAdd);
	toAdd->addParent(currentNode);
}

list<InfoPipe*> NodeTree::returnOrderedNodes() {
	list<InfoPipe*> toReturn = traverse(top);
	return toReturn;
}

list<InfoPipe*> NodeTree::traverse(Node* current) {
	list<InfoPipe*> toReturn;
	if (current->getChild(true) != NULL)
		toReturn = traverse(current->getChild(true));
	toReturn.push_back(current->getContent());
	if (current->getChild(false) != NULL) {
		list<InfoPipe*> temp = traverse(current->getChild(false));
		temp.splice(temp.begin(), toReturn);
		toReturn.splice(toReturn.begin(), temp);
	}
	return toReturn;
}