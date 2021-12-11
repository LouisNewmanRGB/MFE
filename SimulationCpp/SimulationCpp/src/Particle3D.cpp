#include <math.h>
#include <iostream>

#include "AbstractParticle3D.h"
#include "AbstractCompartment.h"

AbstractParticle3D::AbstractParticle3D(double x, double y, double z)
	: Abstract3DObject{ x, y, z },
	m_truePos{m_pos},
	m_startingPos{m_pos},
	m_T2Signal{1.}
{
}

void AbstractParticle3D::setVelocity(Eigen::Vector3d velocity) {
	m_velocity = velocity;
}

void AbstractParticle3D::changeCompartment(std::shared_ptr<AbstractCompartment> compartment, double timeStep) {
	m_compartment = compartment;
	changeSpeed(sqrt(6 * m_compartment->getDiffusivity() / timeStep));
}

void AbstractParticle3D::changeSpeed(double newSpeed) {
	m_velocity = newSpeed * m_velocity / getSpeed();
}