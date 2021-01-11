/*
 * CellCompartment.cpp
 *
 *  Created on: 20 Mar 2020
 *      Author: nw14
 */

#include "CellCompartment.h"
//#include <R.h>
#include <Rmath.h>
#include <string>
using namespace std;
CellCompartment::CellCompartment(int id,
		int targetPopSize,double divisionRate,std::vector<double> fitness):
		id(id),mTargetPopSize(targetPopSize),mDivisionRate(divisionRate),mFitness(fitness){
		nsub=fitness.size();
		for(int i=0;i<nsub;i++){
			vector<shared_ptr<PhyloNode>> tmp;
			subCompartments.push_back(tmp);
			bSubActive.push_back(false);
		}
		totalpop=0;
		rndGen=RandomNumberGenerator::getInstance();
		//Get random number generator
		prob=new double[nsub];
		if(mDivisionRate < -0.001){
			active=false;
		//	printf("Compartment %d is not active\n",id);
		}else{
			active=true;
		//	printf("Compartment %d is active\n",id);
		}
		//setRates();
		//printInfo();
	}

CellCompartment::~CellCompartment() {

	//printf("destroying %d",id);

}
//double rrexp(double lambda){
//  return exp_rand()/lambda;
//}

void CellCompartment::addNode(shared_ptr<PhyloNode> node,int subid){
	subCompartments[subid].push_back(node);
	totalpop++;
	bSubActive[subid]=true;//not ideal dirty way of setting bSubActive flag if we start with non-empty subcompartment
	if(totalpop>=mTargetPopSize){
		atEquilibrium=true;
	}
}

void CellCompartment::setNumNonEmptyIndices(){
	//printf("In setNumNonEmptyIndices\n");
	for(int i=0;i<subCompartments.size();i++){
		if(subCompartments[i].size()>0){
			nonEmptyCompartmentIndices.push_back(i);
		}
	}
	numNonEmptyCompartments=nonEmptyCompartmentIndices.size();

}

void CellCompartment::checkPop(){
	int chk=0;
	for(const vector<shared_ptr<PhyloNode>> & subCompartment : subCompartments){
			 chk+=subCompartment.size();
	}
	if(chk!=totalpop){
		printf("pop mismatch:%d != %d\n",chk,totalpop);
		throw "CellCompartment::checkPop: inconsistency in population tracking";
	}
}

void CellCompartment::printInfo(){

	printf("id=%d\nrate=%5.4f\ntargetPop=%d\nmTotalDeathRate=%7.6f\nmTotalDivRate=%7.6f\n%s\n",id,mDivisionRate,mTargetPopSize,mTotalDeathRate,mTotalDivRate,active ? "ACTIVE" : "INACTIVE");
	for(int k=0;k<nsub;k++){
		printf("fitness[%d]=%5.4f,count[%d]=%d\n",k,mFitness[k],k,subCompartments[k].size());
	}
}


void CellCompartment::clear(){
	for(auto & sub:subCompartments){
		sub.clear();
	}
	totalpop=0;
}

vector<shared_ptr<PhyloNode>> CellCompartment::getNodes(){
	vector<shared_ptr<PhyloNode>> nodes;
	for(const vector<shared_ptr<PhyloNode>> & subCompartment : subCompartments){
				 nodes.insert(nodes.end(),subCompartment.begin(),subCompartment.end());
	}
	return nodes;
}

vector<pair<bool,int>> CellCompartment::getSubCounts(){
	vector<pair<bool,int>> counts;
	int k=0;
	for(const vector<shared_ptr<PhyloNode>> & subCompartment : subCompartments){
				counts.push_back(pair<bool,int>(bSubActive[k],subCompartment.size()));
				k++;
	}
	return counts;

}


double CellCompartment::getTotalDivisionRate(){
	if(mDivisionRate<0){
		return 0.0;
	}
	double tot=0.0;
	int i;
	int ii;
	for(i=0;i<numNonEmptyCompartments;i++){
		ii=nonEmptyCompartmentIndices[i];
		tot+=mDivisionRate*(1+mFitness[ii])*subCompartments[ii].size();
	}
	/*
	for(i=0;i<nsub;i++){
		tot+=mDivisionRate*(1+mFitness[i])*subCompartments[i].size();
	}
    */


	return tot;
}


double CellCompartment::getTotalRate(){
	double tot=mTotalDivRate+mTotalDeathRate;
	if(tot<0.0){
		printInfo();
		throw "getTotalRate: UnexpectedNegative";

	}
	return mTotalDivRate+mTotalDeathRate;
	/**
	if(!active){
		return 0.0;
	}
	double tot=getTotalDivisionRate();
	if(atEquilibrium){
		return 2*tot;//Deaths are equal....  Come back to this...
	}else{
		return tot;
	}*/
}

void CellCompartment::setRates(){
	//We want 20% (alpha) of the deviation to go per day - so rate_delta needs to be rate_delta=0.05*pop_delta
	//printf("In set rates\n");
	if(!active){
		mTotalDivRate=0.0;
		mTotalDeathRate=0.0;
	}else{

		mTotalDivRate=getTotalDivisionRate();
		double tmp=(mTotalDivRate-alpha*(mTargetPopSize-totalpop));
		mTotalDeathRate=tmp>0?tmp:0.0;
	}
}

/* Moved to CellSimulation.cpp
void CellCompartment::doEvent(CellSimulation & sim);
*/
