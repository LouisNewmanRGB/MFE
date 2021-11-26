#include <initializer_list>
#include <iostream>

#include "Environment.h"
#include "AbstractSimulation.h"

Environment::Environment(double T2, double diffusivity, double sizeX, double sizeY, double sizeZ)
	: AbstractCompartment{ 0., 0., 0., T2, diffusivity },
	m_size{ Eigen::Vector3d{sizeX, sizeY, sizeZ} },
	m_aabb{ std::array<Eigen::Vector3d, 2> {-m_size/2, m_size/2} },
	m_tolerance{ std::max({sizeX, sizeY, sizeZ})*AbstractSimulation::m_TOL}
{
}

bool Environment::contains(Eigen::Vector3d pos) {
	return true;
}

double Environment::findIntersection(std::array<Eigen::Vector3d, 2> ray, double maxDistance) {
	Eigen::Vector3d direction = ray[1];
	Eigen::Vector3d dirFraction = {};
	for (int i = 0; i < 3; i++) {
		if (direction[i] == 0.0) {
			dirFraction[i] = std::numeric_limits<double>::infinity();
		}
		else {
			dirFraction[i] = 1 / direction[i];
		}
	}

	double t1 = (m_aabb[0][0] - ray[0][0]) * dirFraction[0];
	double t2 = (m_aabb[1][0] - ray[0][0]) * dirFraction[0];
	double t3 = (m_aabb[0][1] - ray[0][1]) * dirFraction[1];
	double t4 = (m_aabb[1][1] - ray[0][1]) * dirFraction[1];
	double t5 = (m_aabb[0][2] - ray[0][2]) * dirFraction[2];
	double t6 = (m_aabb[1][2] - ray[0][2]) * dirFraction[2];

	double tmin = std::max({ std::min(t1, t2), std::min(t3, t4), std::min(t5, t6) });
	double tmax = std::min({ std::max(t1, t2), std::max(t3, t4), std::max(t5, t6) });

	if (tmax < 0) {
		return std::numeric_limits<double>::infinity();
	}
	if (tmin > tmax) {
		return std::numeric_limits<double>::infinity();
	}
	if (tmin > 0) {
		return tmin;
	}
	else {
		return tmax;
	}

}


void Environment::collide(std::shared_ptr<AbstractParticle3D> particle, Eigen::Vector3d oldPos, Eigen::Vector3d intersection, AbstractSimulation &sim) {
	Eigen::Vector3d line = intersection + particle->getVelocity() / particle->getSpeed()*m_tolerance;
	Eigen::Vector3d newPos = intersection; //making a copy
	for (int i = 0; i < 3; i++) {
		if (line[i] > m_size[i] / 2) {
			newPos[i] = -m_size[i] / 2;
		}
		else if (line[i] < -m_size[i] / 2) {
			newPos[i] = m_size[i] / 2;
			}
	}
	particle->setPos(newPos);
	std::shared_ptr<AbstractCompartment> newComp = sim.findCompartment(particle->getPos());
	particle->changeCompartment(newComp, sim.getTimeStep());
}
