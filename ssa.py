from math import log
from random import random, seed


class SSA:
    """Container for Stochastic Simulation Algorithms"""

    def __init__(self, a=None, version=None):
        """Initializer for SSA methods."""
        if a is not None or version is not None:
            seed(a=a, version=(version or 2))

    def direct(self, model, count=1):
        """Generator of direct-method trajectories"""
        while count > 0:
            while not model.exit():

                # init step: reaction probabilties and partition func
                weights = list((k, v(model)) for k,v in model.propen)
                partition = sum(t[1] for t in weights)
                weights = list((k, v / partition) for k,v in weights)
                
                # monte carlo step: next reaction time
                model["time"].append(
                    log(1.0 / random()) / partition
                )

                # monte carlo step: next reaction
                next_reaction = partition * random()
                curr_reaction = 0.0
                while curr_reaction < next_reaction:
                    curr_reaction += weights.pop()
                reaction_stoich = model.stoich[weights.pop()[0]]

                # update reaction species
                for species, delta in reaction_stoich:
                    model[species] += delta

            yield model.trajectory
            model.reset()
            count -= 1

    def first_reaction(self, model, count=1):
        """Generator of 1st-reaction trajectories"""
        while count > 0:
            while not model.exit():

                # monte carlo step: generate reaction times
                times = list(
                    (k,  log(1.0 / random()) / v(model))
                    for k,v in model.propen
                ).sort(key=lambda t: t[1])

                # update next reaction time
                model["time"].append(times[0][1])

                # update reaction species
                reaction_stoich = model.stoich[times[0][0]]
                for species, delta in reaction_stoich:
                    model[species] += delta

            yield model.trajectory
            model.reset()
            count -= 1

    def first_family(self, model, count=1, families=1):
        """Generator of 1st-family trajectories"""
        while count > 0:
            
            # build list of families
            families_list = list()
            family_size = len(model.propen) // families
            i = 0
            for j in range(families):
                if j < families - 1:
                    families_list.append(
                        (j, model.propen[i:i+family_size])
                    )
                else:
                    families_list.append(
                        (j, model.propen[i:i+family_size])
                    )

            while not model.exit():
                
            yield model.trajectory
            model.reset()
            count -= 1



class Model(dict):
    """
    WIP
    - self is a dict of N species and a temporal param
    - self.stoich is a dict of M stoichiometry lists
    - self.propen is dict of M propensity lambdas
    """
    def __init__(self, initia, propen, stoich):
        self.propen = propen.items().sort(
            key=lambda t: t[0], reverse=True
        )
        self.stoich = stoich.items()
        super().__init__(**initia)

    def exit(self):
        pass
        
    def reset(self):
        for key in self.keys():
            del self[key][1:]