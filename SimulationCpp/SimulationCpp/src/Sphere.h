#ifndef SPHERE_H
#define SPHERE_H

#include "AbstractCompartment.h"

class Simulation;

class Sphere : public AbstractCompartment
{
private:
	double m_permeability;
	double m_radius;

public:
	Sphere(double x, double y, double z, double T2, double diffusivity, double permeability, double radius);

	bool contains(Eigen::Vector3d pos);
	double findIntersection(std::array<Eigen::Vector3d, 2> ray, double maxDistance);
	bool collide(std::shared_ptr<AbstractParticle3D> particle, Eigen::Vector3d oldPos, Eigen::Vector3d intersection, double timeStep);
};

#endif

