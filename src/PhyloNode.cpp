/*
 * PhyloNode.cpp
 *
 *  Created on: 24 Mar 2020
 *      Author: nw14
 */

#include "PhyloNode.h"
PhyloNode::PhyloNode() {
	//parent=NULL;
	//list<PhyloNode *> children;
	ndiv=0;
//	ndriver=0;
	tbirth=0;
	id=0;
}


PhyloNode::~PhyloNode() {
	// TODO Auto-generated destructor stub

}

shared_ptr<PhyloNode> PhyloNode::getSibling(){

	shared_ptr<PhyloNode> parent2=parent.lock();
	if (parent2) {
		for(const shared_ptr<PhyloNode> & child: parent2->children){
			if(child.get()!=this){
				return child;
			}
		}
	}
	return NULL;
}


int PhyloNode::addChild(shared_ptr<PhyloNode> child){
	children.push_back(child);
	return 0;
}
int PhyloNode::removeChild(shared_ptr<PhyloNode> child){
	children.remove(child);
  return 0;
}


int PhyloNode::addEvent(shared_ptr<Event> event){
	events.push_back(event);
	return 0;
}

void PhyloNode::clear(){
//	ndriver=0;
		  //out->nmut=0;
	ndiv=0;
	children.clear();
	events.clear();
	//parent=NULL;
	parent.reset();
	//out->nchild=0;
	tbirth=0.0;
	id=0;
}
