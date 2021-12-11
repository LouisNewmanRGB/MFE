#ifndef ABSTRACTCOMPARTMENT_H
#define ABSTRACTCOMPARTMENT_H

#include <array>
#include <Eigen/Dense>

#include "Abstract3DObject.h"

class AbstractParticle3D;
class Simulation;

class AbstractCompartment : public Abstract3DObject
{
protected:
	double m_T2;
	double m_diffusivity;
	std::shared_ptr<AbstractCompartment> m_parentComp;

public:
	AbstractCompartment(double x, double y, double z, double T2, double diffusivity);

	virtual bool contains(Eigen::Vector3d pos) = 0;
	virtual double findIntersection(std::array<Eigen::Vector3d, 2> ray, double maxDistance) = 0;
	virtual bool collide(std::shared_ptr<AbstractParticle3D> particle, Eigen::Vector3d oldPos, Eigen::Vector3d intersection, double timeStep) = 0;

	void setParentComp(std::shared_ptr<AbstractCompartment> parentComp);

	double getT2() { return m_T2; }
	double getDiffusivity() { return m_diffusivity; }

};

#endif