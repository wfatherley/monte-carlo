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
        """Generator of direct-method trajectories"""
        while count > 0:
            while not model.exit():

                # init step: reaction probabilties and partition func
                probs = list((k, v(model)) for k,v in model.propen)
                partition = sum(t[1] for t in probs)
                probs = list(
                    (k, v / partition) for k,v in probs
                ).reverse()
                
                # monte carlo step: next reaction time
                model["time"].append(
                    log(1.0 / random()) / partition
                )

                # monte carlo step: next reaction
                next_reaction = partition * random()
                curr_reaction = 0.0
                while curr_reaction < next_reaction:
                    curr_reaction += probs.pop()
                reaction_stoich = model.stoich[probs.pop()[0]]
                for k,v in reaction_stoich:
                    model[k] += v
                    
            yield model.trajectory
            model.reset()
            count -= 1

    def first_reaction(self, model):
        """Generator of 1st-reaction trajectories"""
        pass

    def first_family(self, model):
        """Generator of 1st-family trajectories"""
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