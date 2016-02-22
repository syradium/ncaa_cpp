//
// Created by user on 6/28/15.
//

#ifndef DIPLOMA_MEP_H
#define DIPLOMA_MEP_H
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/thread/thread.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/math/tools/roots.hpp>
#include <fstream>
#include <armadillo>
#include <boost/math/constants/constants.hpp>
#include "constants.h"
#include "gradient.h"

void recordPath( const arma::mat& path, const char* fileName);

void recordPoint( const arma::mat& path, const std::vector<std::vector<double>> & pathN, const std::vector<std::vector<double>> & pathM, const std::vector<double> & energies, const char* fileName, const int& pointIndex );

std::vector<double> pathTangent(const int& imageNum, const arma::mat& path, const std::vector<double> & energies);
std::vector<double> calculateStepQuick(const std::vector<double> & perpV, int n, arma::mat& velocities);
std::vector<double> calculateStepQuick2(const std::vector<std::vector<double> >& path, const std::vector<double> & perpV, int n, arma::mat& velocities);
void buildMEP(
			    std::vector<double> thetaAngles, std::vector<double> phiAngles,
				std::vector<double> magneticMoments, std::vector<double> electronsNumber,
              std::vector<double> E0, std::vector<double> U0, std::vector<std::vector<double>> hopingIntegrals, int numberOfImages,
              std::vector<double> initialState, std::vector<double> finalState, std::string fileName);

std::vector<std::vector<double> > findMinima(const std::vector<double> & magneticMoments, const std::vector<double> & electronsNumber,
		const std::vector<double> & E0, const std::vector<double> & U0, const std::vector<std::vector<double>> & hopingIntegrals);

std::vector<std::vector<double> > findMinimaVelocities(const dvector& magneticMoments, const dvector& electronsNumber,
		const dvector& E0, const dvector & U0, const std::vector<dvector> & hopingIntegrals);

std::vector<std::vector<double> > findMinimaVelocitiesThreaded(const dvector& magneticMoments, const dvector& electronsNumber,
		const dvector& E0, const dvector& U0, const std::vector<dvector> & hopingIntegrals);

#endif //DIPLOMA_MEP_H
