from math import log
from random import random, seed


class SSAException(Exception):
    """A simple exception."""
    pass


class SSA:
    """Container for Stochastic Simulation Algorithms"""

    def __init__(self, a=None, version=None):
        """Initializer for all three SSA methods."""
        if a is not None or version is not None:
            seed(a=a, version=version)

    def direct(self, model):
        """Indefinite generator of direct-method trajectories"""

        yield model.trajectory
        model.reset()

    def first_reaction(self, model):
        """Indefinite generator of 1st-reaction trajectories"""
        pass

    def first_family(self, model):
        """Indefinite generator of 1st-family trajectories"""
        pass


class Model(dict):
    def __init__(self, exitlambda=None, **species):
        self.exitlambda = exitlambda
        
        
    def reset(self):
        for key in self.keys():
            del self[key][1:]