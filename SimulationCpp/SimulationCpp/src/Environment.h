#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "AbstractCompartment.h"

class Simulation;

class Environment : public AbstractCompartment
{
private:
	Eigen::Vector3d m_size;
	std::array<Eigen::Vector3d, 2> m_aabb;
	double m_tolerance;

public:
	Environment(double T2, double diffusivity, double sizeX, double sizeY, double sizeZ);

	bool contains(Eigen::Vector3d pos);
	double findIntersection(std::array<Eigen::Vector3d, 2> ray, double maxDistance);
	bool collide(std::shared_ptr<AbstractParticle3D> particle, Eigen::Vector3d oldPos, Eigen::Vector3d intersection, double timeStep);

	Eigen::Vector3d getSize() { return m_size; }
};

#endif