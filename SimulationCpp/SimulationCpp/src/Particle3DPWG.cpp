#include "Particle3DPWG.h"
#include "Simulation.h"

#include <iostream>

Particle3DPWG::Particle3DPWG(double x, double y, double z, std::shared_ptr<std::vector<double>> sequenceTimes, std::shared_ptr<std::vector<Eigen::Vector3d>> sequenceValues)
:AbstractParticle3D{x, y, z},
m_phase{0},
m_totalElapsedTime{0},
m_sequenceIndex{0},
m_maxIndex{(int) size(*sequenceTimes.get()) - 1 },
m_sequenceTime{ (*sequenceTimes.get())[0] },
m_sequenceValue{ (*sequenceValues.get())[0] },
m_sequenceTimes{sequenceTimes},
m_sequenceValues{sequenceValues}
{
}

double Particle3DPWG::integrateGradient(double time) {
	return time * m_sequenceValue.dot(m_truePos + m_velocity * time / 2);
}

void Particle3DPWG::move(double time) {
	double t = m_totalElapsedTime;
	double tFinal = t + time;
	while (m_sequenceIndex < m_maxIndex && m_sequenceTime < tFinal) {
		double dt = m_sequenceTime - t;
		if (!m_sequenceValue.isZero()) {
			m_phase += integrateGradient(dt);
		}
		t = m_sequenceTime;
		m_sequenceIndex++;
		m_sequenceTime = (*m_sequenceTimes.get())[m_sequenceIndex];
		m_sequenceValue = (*m_sequenceValues.get())[m_sequenceIndex];
		m_truePos += dt * m_velocity;
	}
	double dtFinal = tFinal - t;
	if (!m_sequenceValue.isZero()) {
		m_phase += integrateGradient(dtFinal);
	}
	m_pos += time * m_velocity;
	m_truePos += dtFinal * m_velocity;
	m_T2Signal *= exp(-time / m_compartment->getT2()); //T2 signal attenuation
	m_totalElapsedTime += time;
}

std::vector<double> Particle3DPWG::getResults() {
	return std::vector<double> { m_startingPos[0], m_startingPos[1] , m_startingPos[2],
		m_pos[0], m_pos[1], m_pos[2], m_truePos[0], m_truePos[1],  m_truePos[2], getT2Signal(), getPhase()};
}
