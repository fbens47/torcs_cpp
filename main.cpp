#include <iostream>
// #include <omp.h>
#include "carrier.h"
#include <armadillo>
#include <chrono>
#include <ctime>
#include "math.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>



#define CONSTANT_STRING "Hello World!"
#define HELLO() std::cout << CONSTANT_STRING << std::endl;

// https://bisqwit.iki.fi/story/howto/openmp/#PrefaceImportanceOfMultithreading

int main() {
    // HELLO()
    // ~ 22.1s runtime
    double freqMin(2000);
    double freqMax(64000);
    unsigned int nFreq(640);
    double duration(0.5);
    double samplingRate(192000.0);
    unsigned int intSize;
    intSize = sizeTimeMatrix(duration, samplingRate);
    arma::mat spectralModulation;
    spectralModulation << -1.5 << -0.9 << -0.3 << 0 << 0.6 << 1.2 << 1.8 << arma::endr;
    // std::cout << spectralModulation << std::endl;
    arma::mat temporalModulation = 4 * arma::linspace(1, 6, 6);
    arma::mat freq = freqSeries(nFreq, freqMin, freqMax);
    arma::mat oct = octSeries(freq, freqMin);
    arma::mat time = timeSeries(duration, intSize);
    arma::mat phases = makePhases(nFreq);
    arma::mat phasesE = makePhases(6);
    arma::mat carrier = makeCarrierEmpty(intSize, nFreq);
    arma::mat envelops = carrier;
    makeTorcs(carrier, envelops, &phases, &phasesE, time, &freq, &oct, &temporalModulation, 1.2);
    return 0;
}

//std::time_t begin = std::time(nullptr);
//    // std::chrono::time_point<std::chrono::system_clock> start, end;
//    // start = std::chrono::system_clock::now();
//    //
//    // auto start = std::chrono::steady_clock::now();
//    // auto end = std::chrono::steady_clock::now();
//    // auto diff = end - start;
//    // auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(diff).count();
//    // end = std::chrono::system_clock::now();
//    // long long elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
//    // std::cout << "Elapsed Time: " << elapsed_time << std::endl;
//    // std::time_t end = std::time(nullptr);
//    // std::cout << end - begin << std::endl;


/// vector< vector<int> > vvi;
//
//Then you need to use two iterators to traverse it, the first the iterator of the "rows", the second the iterators of the "columns" in that "row":
//
////assuming you have a "2D" vector vvi (vector of vector of int's)
//vector< vector<int> >::iterator row;
//vector<int>::iterator col;
//for (row = vvi.begin(); row != vvi.end(); row++) {
//    for (col = row->begin(); col != row->end(); col++) {
//        // do stuff ...
//    }
//}