from math import log
from random import random, seed


class SSA:
    """Container for Stochastic Simulation Algorithms"""

    def __init__(self, a=None, version=None):
        """Initializer for SSA methods."""
        if a is not None or version is not None:
            seed(a=a, version=(version or 2))

    def direct(self, model):
        """Indefinite generator of direct-method trajectories"""
        while True:
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

    def first_reaction(self, model):
        """Indefinite generator of 1st-reaction trajectories"""
        while True:
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

    def first_family(self, model, families=2):
        """Indefinite generator of 1st-family trajectories"""

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

        while True:
            while not model.exit():

            yield model.trajectory
            model.reset()
            count -= 1



class Epidemic(dict):
    """Simple model to meet design of class SSA"""
    def __init__(
        self,
        initial_conditions: dict,
        propensities: dict,
        stoichiometry: dict
    ):
        self.propen = sorted(
            propensities.items(), reverse=True
        )
        self.stoich = list(stoichiometry.items())
        super().__init__(**initial_conditions)
        
    def exit(self):
        pass
        
    def reset(self):
        for key in self.keys():
            del self[key][1:]


def first_family(model, families=2):
    """Indefinite generator of 1st-family trajectories"""
    # build list of families
    families_list = list()
    family_size = len(model.propen) // families
    i = 0
    for j in range(families):
        if j == families - 1:
            family_size += len(model.propen) % families
        families_list.append(
            (j, model.propen[i:i+family_size])
        )
        i += family_size
        return families_list
        
ic = dict(s=99,i=1,r=0)
p = {
    0: lambda d: 0.5 * d["s"] * d["i"] / sum(d.values()),
    1: lambda d: 1.0 * d["i"],
}
s = {
    0: {"s": -1, "i": 1, "r": 0},
    1: {"s": 0, "i": -1, "r": 1}
}
