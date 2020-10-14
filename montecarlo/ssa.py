from math import log
from random import random
from warnings import warn


class Base(dict):
    """Container for Stochastic Simulation Algorithms"""

    def __init__(
        self,
        initial_conditions,
        propensities,
        stoichiomemtry
    ):
        """Initialize SSA"""
        self.propen = list(propensities.items())
        self.stoich = list(stoichiometry.items())
        super().__init__(**initial_conditions)
        
    def exit(self, *args, **kwargs):
        """Return True if conditions met, else False"""
        raise NotImplementedError
        
    def reset(self, *args, **kwargs):
        """Clean up trajectory on or after exit"""
        raise NotImplementedError
        
    def direct(self):
        """Indefinite generator of direct-method trajectories"""
        while True:
            while not self.exit():
                
                # init step: evaluate propensities and partition
                weights = list((k, v(self)) for k,v in self.propen)
                partition = sum(tup[1] for tup in weights)
                
                # monte carlo step 1: next reaction time
                sojourn = log(1.0 / random()) / partition
                self["time"].append(
                    self["time"][-1] + sojourn
                )
                
                # monte carlo step 2: next reaction
                partition = partition * random()
                j = len(weights) - 1
                while partition >= 0.0:
                    partition -= weights.pop()[1]
                    j -= 1
                reaction_stoich = self.stoich[j][1]
                
                # final step: update reaction species
                for species, delta in reaction_stoich.items():
                    self[species].append(
                        self[species][-1] + delta
                    )
                
            yield self.items()
            self.reset()
            
    def first_reaction(self):
        """Indefinite generator of 1st-reaction trajectories"""
        while True:
            while not self.exit():

                # monte carlo step: generate reaction times
                times = list(
                    (k,  log(1.0 / random()) / v(model))
                    for k,v in self.propen
                ).sort(key=lambda t: t[1])

                # update next reaction time
                model["time"].append(times[0][1])

                # update reaction species
                reaction_stoich = self.stoich[times[0][0]]
                for species, delta in reaction_stoich:
                    self[species] += delta

            yield self.items()
            self.reset()

    def first_family(self, model, families=2):
        """Indefinite generator of 1st-family trajectories"""
        
        # switch ladder for choosing method
        if len(model.propen) < families:
            warn("Too many families, using direct-method")
            return self.direct(model)
        elif families == 1:
            return self.direct(model)
        elif len(model.propen) == families:
            return self.first_reaction(model)
        else:

            # compute family partition info
            family_size = len(model.propen) // families
            orphans = len(model.propen) % families

            # gather family indicies of propensities
            family_indices = list()
            head = 0
            tail = family_size - 1
            for i in range(families):
                family_indicies.append((i, (head, tail)))
                head += family_size - 1
                tail += family_size - 1
            if orphans > 0:
                family_indicies[-1][1] += oprhans

            # main loops
            while True:
                while not model.exit():

                    # monte carlo step: generate random floats
                    random_floats = list(
                        random() for _ in range(familes+1)
                    )

                yield model.items()
                model.reset()


__all__ = ["Base"]