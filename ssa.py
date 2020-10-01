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
            seed(a=a, version=(version or 2))

    def direct(self, model, count=1):
        """Indefinite generator of direct-method trajectories"""
        while count > 0:
            while not model.exit():
                propen = dict()
                for k,v in model.propen.items():
                    propen[k] = v(model)
                partition = sum(propen.values())
                for k in propensities.keys():
                    propen[k] = propen[k] / partition
                model["time"].append(
                    log(1 / random()) / partition
                )
                # WIP propensity MC
            yield model.trajectory
            model.reset()
            count -= 1

    def first_reaction(self, model):
        """Indefinite generator of 1st-reaction trajectories"""
        pass

    def first_family(self, model):
        """Indefinite generator of 1st-family trajectories"""
        pass


class Model(dict):
    """
    - self is a dict of N species and a temporal param
    - self.stoich is a dict of M stoichiometry lists
    - self.propen is dict of M propensity lambdas
    """
    def __init__(self, initia, propen, stoich):
        self.propen = propen.items().sort(key=lambda t: t[0])
        self.stoich = stoich.items().sort(key=lambda t: t[0])
        super().__init__(**initia)

    def exit(self):
        pass
        
    def reset(self):
        for key in self.keys():
            del self[key][1:]