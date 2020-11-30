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
            while not self.model.exit():

                # evaluate next reaction times
                times = [
                    (
                        log(
                            1.0 / self.random.floating()
                        ) / pro(self.model),
                        sto
                    )
                    for (rxn, sto, pro) in self.model.reactions
                ]
                times.sort()

                # evaluate reaction time
                self.model["time"].append(
                    self.model["time"][-1] + times[0][0]
                )

                # evaluate reaction
                for species, delta in times[0][1].items():
                    self.model[species].append(
                        self.model[species][-1] + delta
                    )

                self.model.curate()
            yield self.model
            self.model.reset()


class SSAModel(dict):
    """Container for SSA model"""

    def __init__(
        self,
        initial_conditions,
        propensities,
        stoichiometry,
        max_duration=None
    ):
        """Initialize model"""
        super().__init__(**initial_conditions)
        self.max_duration = max_duration
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
        else: return False

        # return True if there is no time left
        if self.max_duration is not None:
            if self["time"][-1] >= self.max_duration:
                return True

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