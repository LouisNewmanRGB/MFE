#ifndef PLANES_H
#define PLANES_H

#include "AbstractCompartment.h"
class Planes : public AbstractCompartment
{
private:
	double m_halfSpacing;

public:
	Planes(double T2, double diffusivity, double spacing);

	bool contains(Eigen::Vector3d pos);
	double findIntersection(std::array<Eigen::Vector3d, 2> ray, double maxDistance);
	void collide(std::shared_ptr<AbstractParticle3D> particle, Eigen::Vector3d oldPos, Eigen::Vector3d intersection, AbstractSimulation &sim);
};

#endif

