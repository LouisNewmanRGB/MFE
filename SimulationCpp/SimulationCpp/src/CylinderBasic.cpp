#include "CylinderBasic.h"

#include "AbstractParticle3D.h"
#include "Util.h"
#include "Simulation.h"
#include <iostream>
#include <cmath>

CylinderBasic::CylinderBasic(double x, double y, double T2, double diffusivity, double permeability, double radius)
	:AbstractCompartment{ x, y, 0, T2, diffusivity },
	m_permeability{ permeability },
	m_radius{ radius }
{
}

bool CylinderBasic::contains(Eigen::Vector3d pos) {
	return pow((pos[0] - m_pos[0]), 2) + pow((pos[1] - m_pos[1]), 2) < m_radius*m_radius;
}

double CylinderBasic::findIntersection(std::array<Eigen::Vector3d, 2> ray, double maxDistance) {
	double ox = ray[0][0] - m_pos[0];
	double oy = ray[0][1] - m_pos[1];
	if (std::sqrt(pow(ox, 2) + pow(oy, 2) ) - m_radius <= maxDistance) {
	//if ((ray[0] - m_pos).norm() - m_maxRadius <= maxDistance) {
		double a = pow(ray[1][0], 2) + pow(ray[1][1], 2);
		double b = 2*(ray[1][0] * ox + ray[1][1] * oy);
		double c = pow(ox, 2) + pow(oy, 2) - m_radius*m_radius;
		double delta = b * b - 4 * a*c;
		if (delta > 0) {
			double t = (-b - sqrt(delta)) / (2 * a);
			if (t > 0) {
				return t;
			}
			else {
				double t = (-b + sqrt(delta)) / (2 * a);
				if (t > 0) {
					return t;
				}
			}
		}
		else if (delta == 0) {
			double t = -b / (2 * a);
			if (t > 0) {
				return t;
			}
		}
	}
	//no intersection
	return std::numeric_limits<double>::infinity();
}

bool CylinderBasic::collide(std::shared_ptr<AbstractParticle3D> particle, Eigen::Vector3d oldPos, Eigen::Vector3d intersection, double timeStep) {
	double transmissionProba = 4 * m_permeability / particle->getSpeed();
	if (transmissionProba > 1) {
		std::cout << "Error: greater than one transmission probability!" << std::endl;
	}
	assert(transmissionProba <= 1);

	if (Util::getRandomUniform(0, 1) >= transmissionProba) {
		//deflection
		Eigen::Vector3d normal = (intersection - Eigen::Vector3d{ m_pos[0], m_pos[1], intersection[2] }).normalized();
		particle->setVelocity(particle->getVelocity() - 2 * normal*normal.dot(particle->getVelocity()));
	}
	else {
		//transmission
		if (contains(oldPos)) {
			//we are leaving this compartment and going to parent
			particle->changeCompartment(m_parentComp, timeStep);
		}
		else {
			//we are entering this compartment
			particle->changeCompartment(std::make_shared<CylinderBasic>(*this), timeStep);
		}
	}
	return false;
}