from math import log
from random import random, seed
from warnings import warn


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
                
                # monte carlo step 1: next reaction time
                model["time"].append(
                    log(1.0 / random()) / partition
                )

                # monte carlo step 2: next reaction
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

            while True:
                while not model.exit():

                    # monte carlo step: generate random floats
                    random_floats = list(
                        random() for _ in range(familes+1)
                    )

                yield model.trajectory
                model.reset()


class Model(dict):
    """Simple model to meet design of class SSA"""
    def __init__(
        self,
        initial_conditions: dict,
        propensities: dict,
        stoichiometries: dict
    ):
        self.propen = sorted(
            propensities.items(), reverse=True
        )
        self.stoich = list(stoichiometry.items())
        super().__init__(**initial_conditions)

    def exit(self):
        pass
        
    def reset(self):
        pass


__all__ = ["Model", "SSA"]