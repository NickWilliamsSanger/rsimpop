/*
 * PhyloNode.h
 *
 *  Created on: 24 Mar 2020
 *      Author: nw14
 */

#include <list>
#include <memory>
#include <vector>
#include "Event.h"

using namespace std;

#ifndef PHYLONODE_H_
#define PHYLONODE_H_

class PhyloNode {

public:
	PhyloNode();

	virtual ~PhyloNode();
	weak_ptr<PhyloNode> parent;
	list<shared_ptr<PhyloNode>> children;
	int id;//Set as the initial id
	int ndiv;
	//int ndriver;
	double tbirth;
	vector<shared_ptr<Event>> events;
	//int addChild(PhyloNode * child);
	int addChild(shared_ptr<PhyloNode> child);
	int removeChild(shared_ptr<PhyloNode> child);
	int addEvent(shared_ptr<Event> event);
	shared_ptr<PhyloNode> getSibling();
	void clear();

};

#endif /* PHYLONODE_H_ */
