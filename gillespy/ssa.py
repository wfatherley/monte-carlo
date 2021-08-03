"""stochastic simulation algorithm (SSA)"""
from logging import getLogger
from math import inf, log
from random import random, seed
from time import perf_counter

from . import config_seed


__all__ = [
    "Base",
    "Direct",
    "FirstFamily",
    "FirstReaction",
    "NextReaction",
    "TauLeap",
]


LOGGER = getLogger(__name__)


class Base:
    """base SSA"""
    
    def __init__(self, model, **kwargs):
        """construct self"""
        self.model = model
        self.trajectories = kwargs.get("trajectories", inf)
        seed(a=kwargs.get("seed", config_seed))

    def __iter__(self):
        """iterator protocol support"""
        LOGGER.info("returning SSA iterator")
        return self

    def __next__(self):
        """iterator protocol support"""
        while self.trajectories > 0:
            for key in self.model:
                del self.model[key][1:]
            self.trajectories -= 1
            start = perf_counter()
            self.method()
            self.model.perf_time = perf_counter() - start
            return self.model
        raise StopIteration

    def method(self):
        """stochastic simulation algorithm"""
        raise NotImplemented


class Direct(Base):
    """direct-method SSA"""
    
    def method(self):
        """implementation"""
        LOGGER.info(
            "begin direct: trajectories=%i", self.trajectories
        )
        while not self.model.equilibraited():
            weights = [
                (event, stoic, prope(self.model))
                for (event, stoic, prope)
                in self.model.events
            ]
            partition = sum(w[-1] for w in weights)
            sojourn = log(1.0 / random()) / partition
            self.model["time"].append(
                self.model["time"][-1] + sojourn
            )
            partition = partition * random()
            while partition >= 0.0:
                event, stoic, prope = weights.pop()
                partition -= propen
            for species, delta in stoic.items():
                self.model[species].append(
                    self.model[species][-1] + delta
                )
            self.model.update(event)


class FirstFamily(Base):
    """first-family SSA"""
    
    def method(self):
        """implementation"""
        LOGGER.info(
            "begin 1st family: trajectories=%i", self.trajectories
        )


class FirstReaction(Base):
    """first-reaction SSA"""
    
    def method(self):
        """implementation"""
        LOGGER.info(
            "begin 1st rxn: trajectories=%i", self.trajectories
        )
        while not self.model.equilibraited():
            times = [
                (log(1.0 / random()) / pro(self.model), sto)
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
            self.model.update(event)


class NextReaction(Base):
    """next-reaction SSA"""

    def method(self):
        """implementation"""
        LOGGER.info(
            "begin next rxn: trajectories=%i", self.trajectories
        )


class TauLeap(Base):
    """tau-leap SSA"""

    def method(self):
        """implementation"""
        LOGGER.info(
            "begin tau leap: trajectories=%i", self.trajectories
        )
