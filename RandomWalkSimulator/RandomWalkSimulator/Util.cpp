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
	return Eigen::Vector3d{ sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };
}

double Util::getRandomUniform(double a, double b) {
	double random01 = double(rand()) / (double(RAND_MAX) + 1.0);
	return a + random01 * (b - a);
}

double Util::getRandomQuadratic(double radius) {
	//returns a random number drawn from the interval[0, radius] following a quadradic pdf
	return radius * std::cbrt(Util::getRandomUniform(0, 1));
}
	
