/*
 * RandomNumberGenerator.cpp
 *
 *  Created on: 16 Apr 2020
 *      Author: nw14
 */

#include "RandomNumberGenerator.h"
#include <ctime>

RandomNumberGenerator::RandomNumberGenerator(int seed): gen(rd()),d(1.0),u(0.0,1.0){
	//printf("In RandomNumberGenerator Constructor!\n");
	if(seed<0){
		seed=time( NULL );
	}
	printf("RandomNumberGenerator:seed=%d\n",seed);
	this->seed=seed;
	gen.seed(seed);
}

void RandomNumberGenerator::createInstance(int seed){
	//printf("createInstance!\n");
	inst_=new RandomNumberGenerator(seed);
}
// Define the static Singleton pointer
RandomNumberGenerator* RandomNumberGenerator::inst_ = NULL;

RandomNumberGenerator* RandomNumberGenerator::getInstance() {
	//printf("getInstance!\n");
  if (inst_ == NULL) {
  //    inst_ = new RandomNumberGenerator();
	  printf("Should never get here  -- call createInstance(seed)!\n");
	  throw "RandomNumberGenerator::getInstance:Assertion error";
  }
  return(inst_);
}

double RandomNumberGenerator::getExponential(double lambda){
	return d(gen)/lambda;
}
double RandomNumberGenerator::getUniform(){
	return u(gen);
}

int RandomNumberGenerator::sample(int n){
	double rnd=u(gen);
	int val=(int) std::floor(n*rnd);
	if(val==n){  //should never happen!
		printf("RandomNumberGenerator::sample failed bounds check! %d %d %15.14f\n",n,val,rnd);
		throw "RandomNumberGenerator::sample failed bounds check!";
		//exit(EXIT_FAILURE);
		val=n-1;
	}
	return val;
}

int RandomNumberGenerator::sample(int n,double * prob,bool isNormalised){
	double rnd=u(gen);
	if(rnd>maxu){
		maxu=rnd;
	}
	///printf("rnd %7.6f\n",rnd);
	double ptot=0.0;
	if(!isNormalised){
		for(int i=0;i<n;i++){
			ptot+=prob[i];
		}
	}else{
		ptot=1.0;
	}
	double tot=0.0;
	for(int i=0;i<n;i++){
		tot+=(prob[i]/ptot);
		if(rnd<tot){
			return i;
		}
	}
	throw "RandomNumberGenerator::sample(,,) Unexpected path!";
	return -1;
}
