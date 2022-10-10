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
	string pos = "";
public:
	Node(InfoPipe* node);
	string getPos();
	Node* getChild(bool left);
	void addChild(bool left, Node* child);
	void addParent(Node* parent);
};

class NodeTree {
private:
	Node* top;
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

string Node::getPos() {
	return pos;
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
		if (strcmp(toAdd->getPos().c_str(), currentNode->getPos().c_str()) < 0)
			left = true;
		else
			left = false;
		if (currentNode->getChild(left) == nullptr)
			parentFound = true;
		else
			currentNode = currentNode->getChild(left);
	}
	currentNode->addChild(left, toAdd);
}

list<InfoPipe*> NodeTree::returnOrderedNodes() {
	
}