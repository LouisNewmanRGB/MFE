# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _SimulationCpp
else:
    import _SimulationCpp

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


class Simulation(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, nStep, timeStep, envT2, envDiffusivity, envSizeX, envSizeY, envSizeZ):
        _SimulationCpp.Simulation_swiginit(self, _SimulationCpp.new_Simulation(nStep, timeStep, envT2, envDiffusivity, envSizeX, envSizeY, envSizeZ))

    def addSphere(self, x, y, z, T2, diffusivity, permeability, radius):
        return _SimulationCpp.Simulation_addSphere(self, x, y, z, T2, diffusivity, permeability, radius)

    def addEllipsoid(self, *args):
        return _SimulationCpp.Simulation_addEllipsoid(self, *args)

    def addPlanes(self, T2, diffusivity, spacing):
        return _SimulationCpp.Simulation_addPlanes(self, T2, diffusivity, spacing)

    def createStartingPositions(self, nPart, insideCompartments):
        return _SimulationCpp.Simulation_createStartingPositions(self, nPart, insideCompartments)

    def createStartingPositionsAtPos(self, nPart, x, y, z):
        return _SimulationCpp.Simulation_createStartingPositionsAtPos(self, nPart, x, y, z)

    def createSequenceSGP(self):
        return _SimulationCpp.Simulation_createSequenceSGP(self)

    def createSequencePWG(self, sequenceTimes, sequenceVectors):
        return _SimulationCpp.Simulation_createSequencePWG(self, sequenceTimes, sequenceVectors)

    def run(self, seed=-1, partPrintNumber=-1):
        return _SimulationCpp.Simulation_run(self, seed, partPrintNumber)

    def getResults(self):
        return _SimulationCpp.Simulation_getResults(self)
    __swig_destroy__ = _SimulationCpp.delete_Simulation

# Register Simulation in _SimulationCpp:
_SimulationCpp.Simulation_swigregister(Simulation)
cvar = _SimulationCpp.cvar
Simulation.m_TOL = _SimulationCpp.cvar.Simulation_m_TOL



