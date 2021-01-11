/*
 * utils.h
 *
 *  Created on: 23 Mar 2020
 *      Author: nw14
 */

#include <string>
#include <iostream>
#include <random>
#ifndef UTILS_H_
#define UTILS_H_

/*using namespace std;
template <typename T>
string join(const T& v, const string& delim);
*/
template <typename T>


string join(const T& v, const string& delim) {
    ostringstream s;
    for (const auto& i : v) {
        if (&i != &v[0]) {
            s << delim;
        }
        s << i;
    }
    return s.str();
}


//double rrexp(double lambda){
//return exp_rand()/lambda;


//}




#endif /* UTILS_H_ */
