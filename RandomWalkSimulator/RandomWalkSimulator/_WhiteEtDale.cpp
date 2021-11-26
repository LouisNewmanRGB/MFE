#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "Scripts.h"

#include <array>
#include <vector>
#include <npy.hpp>

#include "SimulationFiniteGradient.h"
#include "Sphere.h"
#include "Util.h"

void WhiteEtDale() {
	std::array<double, 9> diffustionTimes = { 12, 18, 24, 30, 36, 42, 48, 54, 60 };
	int nPart = 1e4;
	std::array<double, 8> nucleusFractions = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
	double TE = 140; //ms
	int nStep = 195; //8
	double bValue = 4; //ms / um2
	std::string saveFileName = "../results/white_et_dale_cpp2.npy";

	//environment parameters
	double permeabilityCE = 0;
	double permeabilityNC = 0.5; //0.1 #0.5 #um / ms
	double diffusivityN = 1.31; //um2 / ms
	double diffusivityC = 0.48; //um2 / ms
	double diffusivityExtra = 1.82; //um2 / ms
	double T2N = 63.29; //ms
	double T2C = 23.87; //ms
	double T2Extra = 150; //ms
	double cellRadius = 5; //um

	double envSize = 5 * cellRadius;
	Environment	env = Environment(T2Extra, diffusivityExtra, envSize, envSize, envSize);
	Sphere cytoplasm = Sphere(0, 0, 0, T2C, diffusivityC, permeabilityCE, cellRadius);

	std::vector<double> resultsVector;
	//TODO test with reserve
	//resultsVector.resize(10 * nPart*size(diffustionTimes)*size(nucleusFractions));

	for (double Delta : diffustionTimes) {
		double timeStep = TE / nStep;
		double delta = Delta;
		double delay = TE / 2 - Delta;
		std::vector<double> sequenceTimes = { delay , TE / 2, TE - delay, TE};
		std::vector<Eigen::Vector3d> sequenceValues = { Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{-1, 0, 0}, Eigen::Vector3d{1, 0, 0} , Eigen::Vector3d{0, 0, 0} };
		for (double fraction : nucleusFractions) {
			double nucleusRadius = std::cbrt(fraction) * cellRadius;
			Sphere nucleus = Sphere(0, 0, 0, T2N, diffusivityN, permeabilityNC, nucleusRadius);
			std::vector<Eigen::Vector3d> startingPos(nPart);
			for (int i = 0; i < nPart; i++) {
				startingPos[i] = Util::getRandomDirection()*Util::getRandomQuadratic(cellRadius);
			}
			std::vector<std::shared_ptr<AbstractCompartment>> comp = { std::make_shared<Sphere>(nucleus), std::make_shared<Sphere>(cytoplasm) };

			std::shared_ptr<SimulationFiniteGradient> sim = std::make_shared<SimulationFiniteGradient>(
				SimulationFiniteGradient{ nStep, timeStep, startingPos, env, comp, sequenceTimes, sequenceValues });
			//std::cout << "DONE BUILDING!" << std::endl;
			sim->run();
			std::cout << "DONE RUNNING!" << std::endl;
			std::vector<double> results = sim->getResults();
			resultsVector.insert(resultsVector.end(), results.begin(), results.end());
		}
	}

	//saving results
	std::cout << size(resultsVector) << std::endl;
	const long unsigned shapeArray[] = { size(diffustionTimes), size(nucleusFractions), nPart, 11 };
	npy::SaveArrayAsNumpy(saveFileName, false, 4, shapeArray, resultsVector);
}
