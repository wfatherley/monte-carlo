"""Stochastic Simulation Algorithm (SSA)"""
from math import log

from .random import Mersenne


class SSA:
    """Container for the SSA"""

    def __init__(self, model, seed=1234):
        """Initialize container with model"""
        self.model = model
        self.random = Mersenne(seed=seed)

    def direct(self):
        """Indefinite generator of direct-method trajectories"""
        while True:
            while not self.model.exit():
                weights = [
                    (rxn, sto, pro(self.model))
                    for (rxn, sto, pro)
                    in self.model.events
                ]
                partition = sum(w[-1] for w in weights)
                sojourn = log(
                    1.0 / self.random.floating()
                ) / partition
                self.model["time"].append(
                    self.model["time"][-1] + sojourn
                )
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
                times = [
                    (
                        log(
                            1.0 / self.random.floating()
                        ) / pro(self.model),
                        sto
                    )
                    for (rxn, sto, pro) in self.model.events
                ]
                times.sort()
                self.model["time"].append(
                    self.model["time"][-1] + times[0][0]
                )
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
        self.events = list()
        self.excluded_events = list()
        for event, propensity in propensities.items():
            if propensity(self) == 0.0:
                self.excluded_events.append(
                    (
                        event,
                        stoichiometry[event],
                        propensity
                    )
                )
            else:
                self.events.append(
                    (
                        event,
                        stoichiometry[event],
                        propensity
                    )
                )
        self.events.sort()
        self.excluded_events.sort()

    def exit(self):
        """Return True to break out of trajectory"""
        if len(self.events) == 0:
            return True
        else:
            return False
        if self.max_duration is not None:
            if self["time"][-1] >= self.max_duration:
                return True

    def curate(self):
        """Validate and invalidate elementary events"""
        events = []
        while len(self.events) > 0:
            event = self.events.pop()
            if event[2](self) == 0:
                self.excluded_events.append(event)
            else:
                events.append(event)
        self.events = events
        self.events.sort()
        excluded_events = []
        while len(self.excluded_events) > 0:
            event = self.excluded_events.pop()
            if event[2](self) > 0:
                self.events.append(event)
            else:
                excluded_events.append(event)
        self.excluded_events = excluded_events
        self.excluded_events.sort()

    def reset(self):
        """Clear the trajectory"""
        for key in self:
            del self[key][1:]
        self.curate()


__all__ = ["SSA", "SSAModel"]