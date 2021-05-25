/*
 * rsimpop.cpp
 *
 *  Created on: 25 Mar 2020
 *      Author: nw14
 */

#include <math.h>
#include <map>
#include <list>
#include <set>
#include <utility>
#include "CellSimulation.h"
#include "CellCompartment.h"
#include "PhyloNode.h"
#include "Event.h"
extern "C" {
/**
 * Functions called by R.
 */

void setSimData(vector<Event> & events,vector<shared_ptr<CellCompartment>> & cellCompartments,
		int * eventnode,double * eventts,int * eventval,int * eventdriverid,
		int * nevents,int * compinfoval,double * compinfofitness,
		int * ncomp,int* compartmentsize,double * compartmentrate,
		int * compartmentval,int * ncompartment,int * driverid){
	for(int i=0;i<*nevents;i++){
		events.push_back(Event(eventnode[i],eventts[i],eventval[i],eventdriverid[i]));
	}
	std::map<int, vector<std::pair<double,int>>> fitnessByCompartment;
	for(int i=0;i<*ncomp;i++){
		auto it=fitnessByCompartment.find(compinfoval[i]);
		if(it==fitnessByCompartment.end()){
			vector<std::pair<double,int>> tmp;
			fitnessByCompartment[compinfoval[i]]=tmp;
		}
		fitnessByCompartment[compinfoval[i]].push_back(std::pair<double,int>(compinfofitness[i],driverid[i]));

	}
	for(int i=0;i<*ncompartment;i++){
		if(compartmentval[i]!=i){
			printf("Should be contiguous %d %d\n",compartmentval[i],i);
			throw "compartment should have contiguous ordered value 0..n";
		}
		cellCompartments.push_back(std::make_shared<CellCompartment>(CellCompartment(i,
				compartmentsize[i],
				compartmentrate[i],
				fitnessByCompartment[compartmentval[i]]
				)));
	}
}

void populateEvents(const vector<Event> & eventsOut,int * neventOut,int * eventvalOut,int * eventdriveridOut,double * eventtsOut,int * eventnodeOut){
	int k=0;
	*neventOut=eventsOut.size();
	for(const Event & evento : eventsOut){
				eventvalOut[k]=evento.value;
				eventdriveridOut[k]=evento.driverid;
				eventtsOut[k]=evento.timeStamp;
				eventnodeOut[k++]=evento.node;
	}
}

void populateCompartmentInfo(const CellSimulation & sim,int ncomp,int * nCompPop){
	int k=0;
	for(const auto & compartment : sim.compartments){
		vector<pair<bool,int>> counts=compartment->getSubCounts();
		for(const auto & count: counts ){
			if(k< ncomp){
				nCompPop[k++]=count.second;
			}else{
				throw "populateCompartmentInfo:Too many compartments!";
			}
		}
				//compartment->printInfo();
	}
}
void initRSimPop(int * seed){
	RandomNumberGenerator::createInstance(*seed);
}

void sim_pop2(
		int * edges,
		int * ndivs,
		double * tBirth,
		int * ntips,
		int * nedge,
		int * eventnode,
		int * eventval,
		int * eventdriverid,
		double * eventts,
		int * nevents,
		int * compartmentval,
		double *  compartmentrate,
		int * compartmentsize,
		int * ncompartment,
		int * compinfoval,
		double * compinfofitness,
		int * driverid,
		int * ndriver,
		int *  ncomp,
		double * params,
		int * nparams,
		int * max_size,
		int * max_events,
		int * edgesOut,
		int * nDivsOut,
		int * statusOut,
		int * driverIDOut,
		double * tBirthOut,
		int * nedgeOut,
		int * nInternalNodeOut,
		int * eventvalOut,
		int * eventdriveridOut,
		double * eventtsOut,
		int * eventnodeOut,
		int * neventOut,
		double * timestampOut,
		int * nPopSizeOut,
		int * nDriverOut,
		int * nEventsCount,
		int * ntipsOut,
		int * nCompPop,
		int * status)
{

	try{
		//TODO:
		//RandomNumberGenerator::createInstance(-1);
		int n_days=round(params[0]);
		double startTime=params[3];
		int b_stop_at_pop_size=round(params[2]);
		int b_stop_if_empty=round(params[1]);
		double driverAcquisitionRate=params[4];
		if(driverAcquisitionRate>1e-5){
			throw "driverAcquisitionRate too high!";
		}
		vector<Event> events;
		vector<shared_ptr<CellCompartment>> cellCompartments;
		setSimData(events,cellCompartments,eventnode,eventts,eventval,eventdriverid,nevents,compinfoval,compinfofitness,
				ncomp,compartmentsize,compartmentrate,compartmentval,ncompartment,driverid);
		//printf("allocating cellSim\n");

		CellSimulation sim(edges,
				ndivs,
				tBirth,
				*nedge,
				*ntips,
				events,
				cellCompartments,
				startTime,
				driverAcquisitionRate);
		//Add in Events.
		*status=sim.run((double) n_days,(bool) b_stop_at_pop_size, (bool) b_stop_if_empty);
		sim.snap();
		vector<Event> eventsOut;
		sim.populate_edge_info(edgesOut,nDivsOut,statusOut,driverIDOut,tBirthOut,*max_size,nedgeOut,nInternalNodeOut,ntipsOut,eventsOut);
		populateEvents(eventsOut,neventOut,eventvalOut,eventdriveridOut,eventtsOut,eventnodeOut);
		populateCompartmentInfo(sim,*ncomp,nCompPop);
		auto popTrace=sim.getPopulationTrace();

		int k=0;
		for(tuple<double,int,int> timepoint : popTrace){
			timestampOut[k]=std::get<0>(timepoint);//.first;
			nPopSizeOut[k]=std::get<1>(timepoint);//timepoint.second;
			nDriverOut[k++]=std::get<2>(timepoint);//
		}
		*nEventsCount=popTrace.size();
		//printf("status=%d\n",*status);
		//*status=0;

	} catch (const char* msg) {
		//cerr << msg << endl;
		printf("%s\n",msg);
		*status=-1;
	} catch (const std::exception& ex) {
		printf("%s\n",ex.what());
		*status=-1;
	}

}

/**
 * Restricts the tree (here specified by edges,muts and drivers and times) to nodes connected to the specified subtips.
 */
void sub_sample(int * edges,
		int * ndivs,
		double * tBirth,
		int * ntips,
		int * nedge,
		int * eventnode,
		int * eventval,
		int * eventdriverid,
		double * eventts,
		int * nevents,
		int * compartmentval,
		double *  compartmentrate,
		int * compartmentsize,
		int * ncompartment,
		int * compinfoval,
		double * compinfofitness,
		int * driverid,
		int * ndriver,
		int *  ncomp,
		int * subtips,
		int * nsubtips,
		int * max_size,
		int * edgesOut,
		int * nDivsOut,
		int * stateOut,
		int * driverIDOut,
		double * tBirthOut,
		int * nedgeOut,
		int * nInternalNodeOut,
		double * eventtsOut,
		int * eventvalOut,
		int * eventdriveridOut,
		int * eventnodeOut,
		int * neventOut,
		int * ntipsOut,
		int * nCompPop,
		int * status)
{

	try{
		///RandomNumberGenerator::createInstance(1234567);
		vector<Event> events;
		vector<shared_ptr<CellCompartment>> cellCompartments;
		setSimData(events,cellCompartments,eventnode,eventts,eventval,eventdriverid,nevents,compinfoval,compinfofitness,
				ncomp,compartmentsize,compartmentrate,compartmentval,ncompartment,driverid);
		printf("Specifying %d tips to keep!\n",*nsubtips);
		std::set<int> tipsToDelete;
		int i;
		for(i=0;i<*ntips;i++){
			tipsToDelete.insert(i+1);
		}
		int nkeep=*nsubtips;
		//printf("keeping:");
		for(i=0;i<nkeep;i++){
			tipsToDelete.erase(subtips[i]);
		}

		CellSimulation sim(edges,
				ndivs,
				tBirth,
				*nedge,
				*ntips,
				events,
				cellCompartments,
				0.0,
				0.0);
		sim.deleteTips(tipsToDelete);
		vector<Event> eventsOut;
		sim.populate_edge_info(edgesOut,nDivsOut,stateOut,driverIDOut,tBirthOut,*max_size,nedgeOut,nInternalNodeOut,ntipsOut,eventsOut);
		populateEvents(eventsOut,neventOut,eventvalOut,eventdriveridOut,eventtsOut,eventnodeOut);
		populateCompartmentInfo(sim,*ncomp,nCompPop);
		*status=0;
	} catch (const char* msg) {
		//cerr << msg << endl;
		printf("%s\n",msg);
		*status=1;
	} catch (const std::exception& ex) {
		printf("%s\n",ex.what());
		*status=1;
	}

}





}
