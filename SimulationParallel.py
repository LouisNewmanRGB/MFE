import multiprocessing
import random
import numpy as np

from Simulation import Simulation
from SimulationResults import SimulationResults

class SimulationParallel():
    def __init__(self, nStep, timeStep, particles, environment, compartments=[], nProcess = 8):
        self.nStep = nStep
        self.timeStep = timeStep
        self.particles = particles
        self.nPart = len(particles)
        self.startingPositions = np.array([particle.getTruePos() for particle in particles])
        self.environment = environment
        self.compartments = compartments
        self.nProcess = nProcess
        self.processes = []

        #for results
        self.simulations = None
        self.results = None

    def getResults(self):
        return self.results

    def run(self, seed=None, calcData = False, partPrintNumber = None):
        if seed != None:
            random.seed()

        manager = multiprocessing.Manager()
        result_dict = manager.dict()
        for number in range(self.nProcess):
            partPerProc = self.nPart // self.nProcess
            remainder = self.nPart % self.nProcess
            if number < remainder:
                first, last = number*(partPerProc + 1), (number+1)*(partPerProc + 1)
            else:
                first, last = remainder + number*partPerProc, remainder + (number+1)*partPerProc
            proc = multiprocessing.Process(target=SimulationParallel.runProcess, args=(result_dict, number, partPrintNumber, \
                    self.nStep, self.timeStep, self.particles[first:last], self.environment, self.compartments,))
            #proc = multiprocessing.Process(target=runProcess, args=(number, self.nStep, self.timeStep, last - first, self.environment, self.compartments,))
            proc.start()
            self.processes.append(proc)

        for proc in self.processes:
            proc.join()
        self.simulations = np.array(result_dict.values())

        particles = []
        for sim in self.simulations:
            particles.extend(sim.getParticles())
        self.results = SimulationResults(self.startingPositions, particles)

    #static method
    def runProcess(result_dict, number, partPrintNumber, nStep, timeStep, particles, environment, compartments):
        print("Starting process", number)
        #sim = Simulation(nStep, timeStep, [Particle3D(0, 0, 0) for i in range(nPart)], environment, compartments)
        sim = Simulation(nStep, timeStep, particles, environment, compartments)
        sim.run(seed=None, calcData=False, partPrintNumber=partPrintNumber)
        result_dict[number] = sim
        print("Finished process", number)
