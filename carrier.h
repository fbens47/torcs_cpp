//
// Created by flavien feral on 2018-12-15.
//

#ifndef TORCS_CARRIER_H
#define TORCS_CARRIER_H

#include <vector>
#include "carrier.h"
#include "math.h"
#include <armadillo>

/**
 * Calculate the time matrix (from 0 to total_duration)
 * @param duration: a double
 * @param Size: unsigned int, number of values, from the sampling rate
 * @return: a matrix (nRows=1, nCols=Size), from 0 to duration.
 */
arma::mat timeSeries(double duration, unsigned int Size);

/**
 * Calculate a range of frequencies that will be used
 * @param nFreq: number of frequencies
 * @param freqMin: the minimum frequency
 * @param freqMax: the maximum frequency
 * @return: a matrix (nRows=1, nCols=nFreq), from freqMin to freqMax.
 */
arma::mat freqSeries(unsigned int nFreq, double freqMin, double freqMax);

/**
 * Octave is log2(freq_i / freq_0)
 * @param freq: an armadillo matrix (from freqSeries)
 * @param freqMin: the minimum frequency
 * @return a matrix (same shape as freq Matrix) from 0 to log2(freqMax / freqMin)
 */
arma::mat octSeries(arma::mat &freq, double freqMin);

/**
 * Typo: OneCol /!\, calculatate the carrier values on one column (embedded in a for loop in another
 * function -> makeCarrier), formula is sin(2*pi*f*t)
 * @param time: matrix (nRows=1, nCols=nFreq) of time values
 * @param freq: scalar (arma::as_scalar() or iterator).
 * @param phi: scalar (arma::as_scalar() or iterator).
 * @return: matrix (nRows=1, nCols=nFreq) for one freq.
 */
 // TODO: I should select a random value of phi instead of creating a matrix with nFreq * phi values (memory efficient)
 //  -> OBSOLETE?
arma::mat makeOneRowCarrier(arma::mat &time, double freq, double phi);

/**
 * Calculate the carrier: for each frequency -> loop
 * @param empty: matrix of zeros (nRows=length of time matrix, nCols=nFreq)
 * @param freq: pointer to freq matrix
 * @param phases: pointer to phases matrix
 * @param time: reference to time matrix
 */
void makeCarrier(arma::mat &empty, arma::mat *freq, arma::mat *phases, arma::mat &time);

/**
 * make a matrix filled with zeros
 * @param Size: length of time matrix
 * @param nFreq: number of frequencies
 * @return: array filled with 0
 */
arma::mat makeCarrierEmpty(unsigned int Size, unsigned int nFreq);

/**
 * random values (uniform distribution) between 0 and 2*pi
 * @param nFreq
 * @return matrix (nRows=1, nCols=nFreq) with phases
 */
arma::mat makePhases(unsigned int nFreq);

/**
 * calculate the envelops: a0 + Sum((ajk / 2) * cos(2pi(tm * t + specmod) + phi)) sum on 6 envelops (for loop)
 * @param empty: empty array
 * @param oct: pointer to octave matrix
 * @param phases: pointer to phases matrix
 * @param time: reference to time matrix
 * @param tempMod: pointer to temporal modulation matrix
 * @param specMod: double
 */
void makeEnvelops(arma::mat &empty, arma::mat *oct, arma::mat *phases, arma::mat &time, arma::mat *tempMod,
                  double specMod);

/**
 * Use makeEnvelops and makeCarrier, Schur product between two matrix at the end
 * @param carrier
 * @param envelops
 * @param phasesCarrier
 * @param phasesEnvelops
 * @param time
 * @param freq
 * @param oct
 * @param tempMod
 * @param specMod
 */
void makeTorcs(arma::mat &carrier, arma::mat &envelops, arma::mat *phasesCarrier, arma::mat *phasesEnvelops,
               arma::mat &time, arma::mat *freq, arma::mat *oct, arma::mat *tempMod, double specMod);


/**
 * calculate the size of time matrix
 * @param duration: double
 * @param samplingRate: double
 * @return unisgned int specifying the size of the time matrix.
 */
unsigned int sizeTimeMatrix(double duration, double samplingRate);
#endif //TORCS_CARRIER_H
