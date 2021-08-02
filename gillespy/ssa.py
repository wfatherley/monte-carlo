"""stochastic simulation algorithm (SSA)"""
from logging import getLogger
from math import inf, log
from random import random, seed
from time import perf_counter

from . import config_seed, GillespieException


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
        while not self.model.equilibraited():
            weights = [
                (eve, sto, pro(self.model))
                for (eve, sto, pro)
                in self.model.events
            ]
            partition = sum(w[-1] for w in weights)
            sojourn = log(1.0 / random()) / partition
            self.model["time"].append(
                self.model["time"][-1] + sojourn
            )
            partition = partition * random()
            while partition >= 0.0:
                rxn, sto, pro = weights.pop()
                partition -= pro
            for species, delta in sto.items():
                self.model[species].append(
                    self.model[species][-1] + delta
                )
        return self.model


class FirstFamily(Base):
    """first-family SSA"""
    
    def method(self):
        """implementation"""
        pass


class FirstReaction(Base):
    """first-reaction SSA"""
    
    def method(self):
        """implementation"""
        for key in self.model:
            del self.model[key][1:]
        start = perf_counter()
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
        self.model.perf_time = perf_counter() - start
        return self.model


class NextReaction(Base):
    """next-reaction SSA"""
    pass


class TauLeap(Base):
    """tau-leap SSA"""
    pass
