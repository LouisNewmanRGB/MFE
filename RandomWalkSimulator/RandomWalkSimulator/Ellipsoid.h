#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "AbstractCompartment.h"

class Ellipsoid : public AbstractCompartment
{
private:
	double m_permeability;
	double m_maxRadius;
	Eigen::Matrix3d m_invR;
	Eigen::Matrix3d m_A;

public:
	Ellipsoid(double x, double y, double z, double T2, double diffusivity, double permeability, Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c);

	bool contains(Eigen::Vector3d pos);
	double findIntersection(std::array<Eigen::Vector3d, 2> ray, double maxDistance);
	void collide(std::shared_ptr<AbstractParticle3D> particle, Eigen::Vector3d oldPos, Eigen::Vector3d intersection, AbstractSimulation &sim);
};

#endif