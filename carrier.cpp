//
// Created by flavien feral on 2018-12-15.
//

#include "carrier.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <armadillo>
// http://www.cplusplus.com/reference/cmath/
// http://arma.sourceforge.net/armadillo_nicta_2010.pdf
// https://solarianprogrammer.com/2017/03/24/getting-started-armadillo-cpp-linear-algebra-windows-mac-linux/
// randu(rows, cols) * 2 * M_PI
// linspace(start, end, n) eg: (start O end duration n dt)
// Carrier % Envelops (% Schur product or Hadamard product)
// arma::mat carrier
// arma::mat ripples...
// arma::vec phases
// Advanced Penetration Testing: Hacking the World's Most Secure Networks
// https://www.amazon.fr/Hacker-Playbook-Practical-Penetration-Testing/dp/1980901759/ref=pd_bxgy_14_img_2?_encoding=UTF8&pd_rd_i=1980901759&pd_rd_r=f2bbbac6-009f-11e9-b20d-6137d26eb85e&pd_rd_w=96Si7&pd_rd_wg=WffMU&pf_rd_p=0fca4f02-2421-4308-92a8-78657ba3b2e5&pf_rd_r=BXDXWAQQP3NVMHDNXBV0&psc=1&refRID=BXDXWAQQP3NVMHDNXBV0

// Column major order for Armadillo!!!

unsigned int sizeTimeMatrix(double duration, double samplingRate){
    auto size = (unsigned int) (duration * samplingRate);
    return size;
}

void makeTorcs(arma::mat &carrier, arma::mat &envelops, arma::mat *phasesCarrier, arma::mat *phasesEnvelops,
               arma::mat &time, arma::mat *freq, arma::mat *oct, arma::mat *tempMod, double specMod){
    makeCarrier(carrier, freq, phasesCarrier, time);
    makeEnvelops(envelops, oct, phasesEnvelops, time, tempMod, specMod);
    arma::mat torcs = carrier % envelops;
    torcs = sum(torcs, 1);  // Good here
    std::cout << torcs.n_cols << std::endl;
    std::cout << torcs.n_rows << std::endl;
}


void makeCarrier(arma::mat &empty, arma::mat *freq, arma::mat *phases, arma::mat &time){

    arma::mat::iterator itFreq;
    arma::mat::iterator itPhi;

    unsigned int i(0);
    for (itFreq = freq->begin(), itPhi = phases->begin(); itFreq != freq->end(); ++itFreq, ++itPhi, ++i){
        empty.col(i) = makeOneRowCarrier(time, *itFreq, *itPhi);
    }
}

void makeEnvelops(arma::mat &empty, arma::mat *oct, arma::mat *phases, arma::mat &time, arma::mat *tempMod,
                  double specMod){
    double ajk(30./6);
    double a0(30);
    arma::mat::iterator itOct;
    unsigned int n(0);
    for (itOct = oct->begin(); itOct != oct->end(); ++itOct, ++n){
        arma::mat tmp = a0 * arma::ones(time.n_elem);
        arma::mat::iterator itPhi;
        arma::mat::iterator itTMod;
        for (itPhi = phases->begin(), itTMod = tempMod->begin(); itPhi != phases->end(); ++itPhi, ++itTMod){
            tmp += (ajk/2.) * arma::cos(2 * M_PI * (*itTMod * time + specMod * *itOct) + *itPhi);

        }
        empty.col(n) = 20 * arma::log(tmp);
    }
}

arma::mat makeCarrierEmpty(unsigned int Size, unsigned int nFreq){
    arma::mat empty;
    empty.resize(Size, nFreq);
    empty.fill(0);
    return empty;
}

arma::mat makePhases(unsigned int nFreq){
    arma::mat phases = 2 * M_PI * arma::randu(nFreq);
    return phases;
}

arma::mat makeOneRowCarrier(arma::mat &time, double freq, double phi){
    arma::mat row = arma::sin(2 * M_PI * time * freq + phi);
    return row;
}


arma::mat timeSeries(double duration, unsigned int Size){
    arma::mat timeMat = arma::linspace(0, duration, Size);
    return timeMat;
}

arma::mat freqSeries(unsigned int nFreq, double freqMin, double freqMax){
    arma::mat freqMat = arma::linspace(freqMin, freqMax, nFreq);
    return freqMat;
}

arma::mat octSeries(arma::mat &freq, double freqMin){
    arma::mat oct = arma::log2(freq / freqMin);
    return oct;
}







