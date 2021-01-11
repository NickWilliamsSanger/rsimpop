/*
 * RandomNumberGenerator.h
 *
 *  TODO: Fix up to better conform to singleton pattern.
 *  Created on: 16 Apr 2020
 *      Author: nw14
 */
#include <random>

#ifndef RANDOMNUMBERGENERATOR_H_
#define RANDOMNUMBERGENERATOR_H_

//using namespace std;
//TODO: incorporate R based random variables..

class RandomNumberGenerator {
public:
   // This is how clients can access the single instance
   static RandomNumberGenerator* getInstance();
   static void createInstance(int seed);
   double getExponential(double lambda);
   int sample(int n);
   int sample(int n,double * prob,bool isNormalised);
   double getUniform();
   double maxu=0.0;

protected:
   int seed;

private:
   static RandomNumberGenerator* inst_;   // The single instance
   RandomNumberGenerator(int seed);// private constructor
   RandomNumberGenerator(const RandomNumberGenerator&);
   RandomNumberGenerator& operator=(const RandomNumberGenerator&);
   std::random_device rd;
   std::mt19937 gen;
   std::exponential_distribution<double> d;
   std::uniform_real_distribution<double> u;
};


#endif /* RANDOMNUMBERGENERATOR_H_ */
