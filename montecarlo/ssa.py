"""Stochastic Simulation Algorithm (SSA)"""
from math import log
from random import random

from . import MCException, Mersenne


class SSA:
    """Container for SSAs"""

    def __init__(self, model):
        """Initialize container with model"""
        self.random = Mersenne()
        self.model = model
        
    def direct(self):
        """Indefinite generator of direct-method trajectories"""
        while True:
            while not self.model.exit():
                weights = list(
                    propensity(self.model) for propensity
                    in self.model.propensities.values()
                )
                partition = sum(w for w in weights)
                sojourn = log(
                    1.0 / self.random.floating()
                ) / partition
                self.model["time"].append(
                    self.model["time"][-1] + sojourn
                )
                partition = partition * self.random.floating()
                j = 0
                while partition >= 0.0:
                    partition -= weights.pop(0)
                    j += 1
                reaction_stoich = self.model.stoichiometry[j]
                for species, delta in reaction_stoich.items():
                    self[species].append(
                        self[species][-1] + delta
                    )
            yield self.model.items()
            self.model.reset()
            
    def first_reaction(self):
        """Indefinite generator of 1st-reaction trajectories"""
        while True:
            while not self.exit():
                times = list(
                    (k,  log(1.0 / random()) / v(model))
                    for k,v in self.propen
                ).sort(key=lambda t: t[1])
                model["time"].append(times[0][1])
                reaction_stoich = self.stoich[times[0][0]]
                for species, delta in reaction_stoich:
                    self[species] += delta
            yield self.items()
            self.reset()

    def first_family(self, model, families=2):
        """Indefinite generator of 1st-family trajectories"""
        if len(model.propen) < families:
            warn("Too many families, using direct-method")
            return self.direct(model)
        elif families == 1:
            return self.direct(model)
        elif len(model.propen) == families:
            return self.first_reaction(model)
        else:
            family_size = len(model.propen) // families
            orphans = len(model.propen) % families
            family_indices = list()
            head = 0
            tail = family_size - 1
            for i in range(families):
                family_indicies.append((i, (head, tail)))
                head += family_size - 1
                tail += family_size - 1
            if orphans > 0:
                family_indicies[-1][1] += oprhans
            while True:
                while not model.exit():

                    # monte carlo step: generate random floats
                    random_floats = list(
                        random() for _ in range(familes+1)
                    )

                yield model.items()
                model.reset()


class SSAModel(dict):
    """Container for SSA model"""

    def __init__(
        self, initial_conditions, propensities, stoichiometry
    ):
        """Initialize model"""
        for k,v in propensities.items():
            stoichiometry[k] = (
                stoichiometry[k], v
            )
        self.valid_reactions = sorted(stoichiometry.items())
        self.excluded_reactions = dict()
        super().__init__(**initial_conditions)

    def exit(self):
        """Return True to break out of trajectory"""
        if len(self.valid_reactions) == 0:
            return True
        else:
            return False

    @property
    def propensities(self):
        """Return map of valid propensities"""
        pass

    def reset(self):
        """Clear the trajectory"""
        for key in self:
            del self[key][1:]

    @property
    def stoichiometry(self):
        """Return stoichiometry map"""
        pass


__all__ = ["SSA", "SSAModel"]