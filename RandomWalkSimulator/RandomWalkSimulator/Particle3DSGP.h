#ifndef PARTICLE3DSGP_H
#define PARTICLE3DSGP_H

#include "AbstractParticle3D.h"

class Particle3DSGP : public AbstractParticle3D
{
public:
	Particle3DSGP(double x, double y, double z);

	void move(double time);
	std::vector<double> getResults();
};

#endif

