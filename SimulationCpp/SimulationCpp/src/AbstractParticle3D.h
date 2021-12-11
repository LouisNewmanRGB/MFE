#ifndef PARTICLE3D_H
#define PARTICLE3D_H

#include <Eigen/Dense>
#include <vector>

#include "Abstract3DObject.h"

class AbstractCompartment;

class AbstractParticle3D : public Abstract3DObject
{
protected:
	Eigen::Vector3d m_velocity;
	Eigen::Vector3d m_truePos;
	Eigen::Vector3d m_startingPos;
	std::shared_ptr<AbstractCompartment> m_compartment;
	double m_T2Signal;

public:
	AbstractParticle3D(double x, double y, double z);

	void setVelocity(Eigen::Vector3d velocity);
	void changeCompartment(std::shared_ptr<AbstractCompartment> compartment, double timeStep);
	void changeSpeed(double newSpeed);

	virtual void move(double time) = 0;
	virtual std::vector<double> getResults() = 0;

	std::shared_ptr<AbstractCompartment> getCompartment() { return m_compartment; }
	Eigen::Vector3d getVelocity() { return m_velocity; }
	Eigen::Vector3d getTruePos() { return m_truePos; }
	double getSpeed() { return m_velocity.norm(); }
	double getT2Signal() { return m_T2Signal; }
};

#endif
