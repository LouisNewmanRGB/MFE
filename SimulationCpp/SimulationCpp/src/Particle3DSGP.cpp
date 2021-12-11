#include "Particle3DSGP.h"
#include "AbstractCompartment.h"
#include <iostream>

Particle3DSGP::Particle3DSGP(double x, double y, double z)
:AbstractParticle3D(x, y, z)
{
}

void Particle3DSGP::move(double time) {
	m_pos += time * m_velocity;
	m_truePos += time * m_velocity;
	m_T2Signal *= exp(-time / m_compartment->getT2()); //T2 signal attenuation
}

std::vector<double> Particle3DSGP::getResults() {
	return std::vector<double> { m_startingPos[0], m_startingPos[1], m_startingPos[2],
		m_pos[0], m_pos[1], m_pos[2], m_truePos[0], m_truePos[1], m_truePos[2], getT2Signal()};
}
