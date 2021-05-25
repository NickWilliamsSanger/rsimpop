/*
 * CellCompartment.h
 *
 *  Created on: 20 Mar 2020
 *      Author: nw14
 */

#include <vector>
#include <iostream>
#include <map>
#include "PhyloNode.h"
#include "RandomNumberGenerator.h"

#ifndef CELLCOMPARTMENT_H_
#define CELLCOMPARTMENT_H_
using namespace std;
class CellSimulation;
class CellCompartment {

private:
	vector<vector<shared_ptr<PhyloNode>>> subCompartments;
	vector<bool> bSubActive;
	RandomNumberGenerator * rndGen;
	double * prob;
	vector<int> nonEmptyCompartmentIndices;
	vector<int> emptyCompartmentIndices;
	int numNonEmptyCompartments;
	int numEmptyCompartments;
	std::map<int, int> idxByID;




public:
	CellCompartment(int id,int targetPopSize,double divisionRate,std::vector<std::pair<double,int>> fitnessID);
	virtual ~CellCompartment();
	int id;
	int mTargetPopSize;
	double mDivisionRate;
	double mTotalDivRate=0.0;
	double mTotalDeathRate=0.0;
	double alpha=0.2;
	std::vector<double> mFitness;
	int nsub;
	int totalpop;
	bool active=true;
	bool atEquilibrium=false;
	void addNode(shared_ptr<PhyloNode> node,int subid);
	double getTotalRate();
	double getTotalDivisionRate();
	void setRates();
	void setNumNonEmptyIndices();
	int getSub(int ID);





	void doEvent(CellSimulation & sim);
	void checkPop();
	int die(shared_ptr<PhyloNode> node,CellSimulation & sim );
	vector<shared_ptr<PhyloNode>> getNodes();
	void clear();
	vector<pair<bool,int>> getSubCounts();
	void printInfo();
    //friend ostream & operator << (ostream &out, const CellCompartment &c);
	void addDriver(CellSimulation & sim,double ts,double fitness);
};
//ostream & operator << (ostream &out, const CellCompartment &c);
#endif /* CELLCOMPARTMENT_H_ */
