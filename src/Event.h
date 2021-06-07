/*
 * Event.h
 *
 *  Created on: 1 Apr 2020
 *      Author: nw14
 */

#include <string>

#ifndef EVENT_H_
#define EVENT_H_

class Event {

public:
	Event(int node,double timeStamp,int value,int driverid,int uid):
			node(node),timeStamp(timeStamp),value(value),driverid(driverid),uid(uid){};
	int node;
	double timeStamp;

	int value;//Compartment
	int driverid;//Driver id
	int uid; //unique id
	//string toString(){ return sprintf("%d:%d",node,value)};
};

#endif /* EVENT_H_ */
