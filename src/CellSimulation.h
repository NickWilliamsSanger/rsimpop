/*
 * CellSimulation.h
 *
 *  Created on: 24 Mar 2020
 *      Author: nw14
 */
#include <stack>
#include <vector>
#include <memory>
#include <set>
#include "PhyloNode.h"
#include "CellCompartment.h"
#include "RandomNumberGenerator.h"
#include "Event.h"

using namespace std;

#ifndef CELLSIMULATION_H_
#define CELLSIMULATION_H_

class CellSimulation {
public:
	CellSimulation(int * edges,
	            int * ndivs,
	            double * tBirth,
				int nedge,
				int ntips,
				vector<Event> events,
				vector<shared_ptr<CellCompartment>> cellCompartments,
				double startTime,
				double driverRate);

	virtual ~CellSimulation();

	stack<shared_ptr<PhyloNode>> nodeStack;
	shared_ptr<PhyloNode> root;
	vector<shared_ptr<CellCompartment>> compartments;


	shared_ptr<PhyloNode> getNewPhyloNode();
	double currentTime;
	void populate_edge_info(int * edgesOut,int * ndivsOut,int * stateOut,int * driverIDOut,
			double * tBirthOut,int max_nedge,int * nedgeOut,int *nInternalNodes,
			int * numTip,vector<Event> & eventsOut);
	int run(double stopTime,bool bStopAtEquilibrium,bool bStopIfEmpty);
	void recycle(shared_ptr<PhyloNode> node);
	vector<tuple<double,int,int>> getPopulationTrace();
	void snap();
	void die(shared_ptr<PhyloNode>);
	pair<shared_ptr<PhyloNode>,shared_ptr<PhyloNode>> divide(shared_ptr<PhyloNode> node);
	void deleteTips(std::set<int> tipsToDelete);
	void setCurrentDriverID(int id);
	/**
	 * Retrieves current driver ID and increments.
	 */
	int incrementCurrentDriverID();

private:

	/**
	 * The following recursively updates the edge matrix. The idea is that the root has number #tips+1 and each tip in the range 1:tips and then
	 * internal nodes are > #tips+1
	 */
	void update_edge_matrix(int * edges,
			int * ndivsOut,
			int * stateOut,
			int * driverIDOut,
			double *tBirthOut,
			int parent_counter,
			int * internal_counter,
			int * tip_counter,
			int * row_counter,
			shared_ptr<PhyloNode> thisNode,
			int max_size,vector<Event> & eventsOut,int parentState,int driverIDState);
	void assignCellCompartments();
	void setCompartmentInfoRecursively(shared_ptr<PhyloNode> thisNode,int compartment,int drivers,int * tip_idx,int ntip);
	void setCompartmentInfo();
	void resetNtips();
	void resetNtips(bool bCheckMatch);
	RandomNumberGenerator * rndGen;
	int ntips;
	int ndrivers;
	double driverRate;
	int driverID=1;
	//double * timeStamp;
	//int * popTrace;
	vector<tuple<double,int,int>> populationTrace;;///TODO implement at compartment/driver level

};

#endif /* CELLSIMULATION_H_ */
