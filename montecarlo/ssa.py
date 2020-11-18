"""Stochastic Simulation Algorithm (SSA)"""
from math import log

from .random import Mersenne


class SSA:
    """Container for SSAs"""

    def __init__(self, model, seed=1234):
        """Initialize container with model"""
        self.model = model
        self.random = Mersenne(seed=seed)

    def direct(self):
        """Indefinite generator of direct-method trajectories"""
        while True:
            while not self.model.exit():
                
                # evaluate weights and partition
                weights = [
                    (rxn, sto, pro(self.model))
                    for (rxn, sto, pro)
                    in self.model.reactions
                ]
                partition = sum(w[-1] for w in weights)

                # evaluate sojourn time
                sojourn = log(
                    1.0 / self.random.floating()
                ) / partition
                self.model["time"].append(
                    self.model["time"][-1] + sojourn
                )

                # evaluate the reaction
                partition = partition * self.random.floating()
                while partition >= 0.0:
                    rxn, sto, pro = weights.pop(0)
                    partition -= pro
                for species, delta in sto.items():
                    self.model[species].append(
                        self.model[species][-1] + delta
                    )

                self.model.curate()
            yield self.model
            self.model.reset()
            
    def first_reaction(self):
        """Indefinite generator of 1st-reaction trajectories"""
        while True:
            while not self.exit():
                times = [
                    (
                        k,
                        log(
                            1.0 / self.random.floating()
                        ) / pro(self.model)
                    )
                    for (rxn, sto, pro) in self.reactions
                ]
                times.sort(key=lambda t: t[1])
                self.model["time"].append(times[0][1])
                reaction_stoich = self.stoich[times[0][0]]
                for species, delta in reaction_stoich:
                    self[species] += delta
            yield self.model
            self.reset()

    # def first_family(self, model, families=2):
    #     """Indefinite generator of 1st-family trajectories"""
    #     if len(model.propen) < families:
    #         warn("Too many families, using direct-method")
    #         return self.direct(model)
    #     elif families == 1:
    #         return self.direct(model)
    #     elif len(model.propen) == families:
    #         return self.first_reaction(model)
    #     else:
    #         family_size = len(model.propen) // families
    #         orphans = len(model.propen) % families
    #         family_indices = list()
    #         head = 0
    #         tail = family_size - 1
    #         for i in range(families):
    #             family_indicies.append((i, (head, tail)))
    #             head += family_size - 1
    #             tail += family_size - 1
    #         if orphans > 0:
    #             family_indicies[-1][1] += oprhans
    #         while True:
    #             while not model.exit():

    #                 # monte carlo step: generate random floats
    #                 random_floats = list(
    #                     random() for _ in range(familes+1)
    #                 )

    #             yield model.items()
    #             model.reset()


class SSAModel(dict):
    """Container for SSA model"""

    def __init__(
        self, initial_conditions, propensities, stoichiometry
    ):
        """Initialize model"""
        super().__init__(**initial_conditions)
        self.reactions = list()
        self.excluded_reactions = list()
        for reaction,propensity in propensities.items():
            if propensity(self) == 0.0:
                self.excluded_reactions.append(
                    (
                        reaction,
                        stoichiometry[reaction],
                        propensity
                    )
                )
            else:
                self.reactions.append(
                    (
                        reaction,
                        stoichiometry[reaction],
                        propensity
                    )
                )

    def exit(self):
        """Return True to break out of trajectory"""

        # return True if no more reactions
        if len(self.reactions) == 0: return True

        # return False if there are more reactions
        else: return False

    def curate(self):
        """Validate and invalidate model reactions"""
        
        # evaulate possible reactions
        reactions = []
        while len(self.reactions) > 0:
            reaction = self.reactions.pop()
            if reaction[2](self) == 0:
                self.excluded_reactions.append(reaction)
            else:
                reactions.append(reaction)
        reactions.sort()
        self.reactions = reactions

        # evaluate impossible reactions
        excluded_reactions = []
        while len(self.excluded_reactions) > 0:
            reaction = self.excluded_reactions.pop()
            if reaction[2](self) > 0:
                self.reactions.append(reaction)
            else:
                excluded_reactions.append(reaction)
        excluded_reactions.sort()
        self.excluded_reactions = excluded_reactions

    def reset(self):
        """Clear the trajectory"""

        # reset species to initial conditions
        for key in self: del self[key][1:]

        # reset reactions per initial conditions
        self.curate()


__all__ = ["SSA", "SSAModel"]