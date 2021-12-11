#ifndef UTIL_H
#define UTIL_H

#include <Eigen/Dense>

class Util
{
public:
	static Eigen::Vector3d getRandomDirection();
	static Eigen::Vector3d getRandomUniformPosition(Eigen::Vector3d envSize);
	static double getRandomUniform(double a, double b);
	static double getRandomQuadratic(double radius);
};

#endif

