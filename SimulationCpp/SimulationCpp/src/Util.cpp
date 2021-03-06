#include <cstdlib> //for random
#include <Eigen/Dense>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846 //TODO IMPROVE
#endif

#include "Util.h"

Eigen::Vector3d Util::getRandomDirection() {
	double phi = Util::getRandomUniform(0, 2 * M_PI);
	double cosTheta = Util::getRandomUniform(-1, 1);
	double theta = acos(cosTheta);
	//double theta = Util::getRandomUniform(0, M_PI);
	return Eigen::Vector3d{ sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };
	//double x = Util::getRandomUniform(-1, 1);
	//double y = Util::getRandomUniform(-1, 1);
	//double z = Util::getRandomUniform(-1, 1);
	//double n = sqrt(x*x + y * y + z * z);
	//return Eigen::Vector3d{ x / n, y / n, z / n };
}

Eigen::Vector3d Util::getRandomUniformPosition(Eigen::Vector3d envSize) {
	return Eigen::Vector3d{ Util::getRandomUniform(-envSize[0] / 2, envSize[0] / 2),
				Util::getRandomUniform(-envSize[1] / 2, envSize[1] / 2) , Util::getRandomUniform(-envSize[2] / 2, envSize[2] / 2) };
}

double Util::getRandomUniform(double a, double b) {
	double random01 = double(rand()) / (double(RAND_MAX) + 1.0);
	return a + random01 * (b - a);
}

double Util::getRandomQuadratic(double radius) {
	//returns a random number drawn from the interval[0, radius] following a quadradic pdf
	return radius * std::cbrt(Util::getRandomUniform(0, 1));
}
	
