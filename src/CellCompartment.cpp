/*
 * CellCompartment.cpp
 *
 *  Created on: 20 Mar 2020
 *      Author: nw14
 */

#include "CellCompartment.h"
#include "CellSimulation.h"
//#include <R.h>
#include <Rmath.h>
#include <string>
#include "Event.h"
using namespace std;
CellCompartment::CellCompartment(int id,
                                 int targetPopSize,double divisionRate,std::vector<std::pair<double,int>> fitnessID):
  id(id),mTargetPopSize(targetPopSize),mDivisionRate(divisionRate){
  nsub=fitnessID.size();
  for(int i=0;i<nsub;i++){
    idxByID[fitnessID[i].second]=i;
    mFitness.push_back(fitnessID[i].first);
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
  }else{
    active=true;
  }
}

CellCompartment::~CellCompartment() {

	//printf("destroying %d",id);

}

//double rrexp(double lambda){
//  return exp_rand()/lambda;
//}


int CellCompartment::getSub(int ID){
	return idxByID[ID];
}


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
	nonEmptyCompartmentIndices.clear();
	emptyCompartmentIndices.clear();
	for(int i=0;i<subCompartments.size();i++){
		if(subCompartments[i].size()>0){
			nonEmptyCompartmentIndices.push_back(i);
		}else{
			emptyCompartmentIndices.push_back(i);
		}
	}
	numNonEmptyCompartments=nonEmptyCompartmentIndices.size();
	numEmptyCompartments=emptyCompartmentIndices.size();


}

void CellCompartment::checkPop(){
	int chk=0;
	int k=0;
	for(const vector<shared_ptr<PhyloNode>> & subCompartment : subCompartments){
		// printf("%d=%d,",k++,subCompartment.size());
			 chk+=subCompartment.size();
	}
    //printf("::totals=%d %d\n",chk,totalpop);
	if(chk!=totalpop){
		printf("pop mismatch:%d != %d\n",chk,totalpop);
		throw "CellCompartment::checkPop: inconsistency in population tracking";
	}
}

void CellCompartment::printInfo(){
	printf("id=%d\nrate=%5.4f\ntargetPop=%d\nmTotalDeathRate=%7.6f\nmTotalDivRate=%7.6f\n%s\n",id,mDivisionRate,mTargetPopSize,mTotalDeathRate,mTotalDivRate,active ? "ACTIVE" : "INACTIVE");
	for(int k=0;k<nsub;k++){
		printf("fitness[%d]=%5.4f,count[%d]=%lu\n",k,mFitness[k],k,subCompartments[k].size());
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
	return tot;
}


double CellCompartment::getTotalRate(){
	double tot=mTotalDivRate+mTotalDeathRate;
	if(tot<0.0){
		printInfo();
		throw "getTotalRate: UnexpectedNegative";

	}
	if(atEquilibrium){
		return mTotalDivRate+mTotalDeathRate;
	}else{
	  //Following means we have pure exponential growth until equilibrium
		return mTotalDivRate;
	}
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


void CellCompartment::addDriver(CellSimulation & sim,double ts,double fitness){
	checkPop();
	setNumNonEmptyIndices();
	double prob[nsub];
	for(int i=0;i<numNonEmptyCompartments;i++){
		prob[i]=subCompartments[nonEmptyCompartmentIndices[i]].size();
	}
	int ii=rndGen->sample(numNonEmptyCompartments,prob,false);
	int i=nonEmptyCompartmentIndices[ii];
	int sz=subCompartments[i].size();
	int k=rndGen->sample(sz);
	shared_ptr<PhyloNode> node=subCompartments[i][k];
	//Node ID of event only matters when provided externally.
	int driverid=sim.incrementCurrentDriverID();
	Event event(-1,ts,id,driverid);
	shared_ptr<Event> thisEvent=std::make_shared<Event>(event);
	node->addEvent(thisEvent);
	//Create a new sub-compartment for this node and remove from old compartment..

	fitness=fitness+mFitness[i];//Additive model is HARD coded here.
	if(numEmptyCompartments==0){
		subCompartments.push_back(vector<shared_ptr<PhyloNode>>());
		subCompartments[nsub].push_back(node);
		mFitness.push_back(fitness);
		idxByID[driverid]=sz;
		bSubActive.push_back(true);
		nsub++;
		printf("adding new compartment: %d fit=%3.2f\n",nsub,fitness);
	}else{
		int l=emptyCompartmentIndices[0];
		subCompartments[l].push_back(node);
		mFitness[l]=fitness;
		bSubActive[l]=true;
		idxByID[driverid]=l;
		printf("recycling compartment: %d fit=%3.2f\n",nsub,fitness);
	}
    //Remove from vector in old compartment
	if(k<sz-1){
		subCompartments[i][k]=subCompartments[i][sz-1];
	}
	subCompartments[i].pop_back();
	//checkPop();
	//Finally reset nonEmptyCompartments etc
	setNumNonEmptyIndices();
	checkPop();

}




