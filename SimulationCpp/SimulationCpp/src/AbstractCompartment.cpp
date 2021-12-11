#include "AbstractCompartment.h"

AbstractCompartment::AbstractCompartment(double x, double y, double z, double T2, double diffusivity)
	: Abstract3DObject{ x, y, z },
	m_T2{ T2 },
	m_diffusivity{ diffusivity }
{
}

void AbstractCompartment::setParentComp(std::shared_ptr<AbstractCompartment> parentComp) {
	m_parentComp = parentComp;
}