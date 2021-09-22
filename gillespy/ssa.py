"""stochastic simulation algorithm objects"""
from logging import getLogger
from math import inf, log
from random import random, seed
from time import perf_counter

from .util import default_seed


logger = getLogger(__name__)


class Base:
    """base stochastic simulation algorithm"""
    
    def __init__(self, model, **kwargs):
        """construct self"""
        self.model = model
        self.trajectories = kwargs.get("trajectories", inf)
        self.seed = kwargs.get("seed", default_seed)
        seed(a=self.seed)
        logger.info(
            "seed set: model_id=%s, seed_value=%s",
            model.id,
            str(seed)
        )

    def __iter__(self):
        """iterator protocol support"""
        logger.info(
            "returning simulation iterator: model_id=%s",
            self.model.id
        )
        return self

    def __next__(self):
        """iterator protocol support"""
        while self.trajectories > 0:
            for key in self.model:
                del self.model[key][1:]
            self.trajectories -= 1
            start = perf_counter()
            self.method()
            logger.info(
                "simulation complete: model_id=%s, perf_time=%i",
                model.id,
                perf_counter() - start
            )
            return self.model
        raise StopIteration

    def method(self):
        """stochastic simulation algorithm"""
        raise NotImplemented


class Direct(Base):
    """direct-method stochastic simulation algorithm"""
    
    def method(self):
        """implementation"""
        while not self.model.equilibriated():
            weights = [
                (event, stoic, prope(self.model))
                for (event, stoic, prope)
                in self.model["valid_events"]
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
    """first-family stochastic simulation algorithm"""

    def __init__(self, model, **kwargs):
        """construct self"""
        super().__init__(model, **kwargs)
        self.family_count = kwargs.get("family_count", 2)
    
    def method(self):
        """implementation"""
        while not self.model.equilibriated():
            families = list()
            family_size = len(
                self.model["valid_events"]
            ) // family_count
            for i in range(self.family_count):
                j = i * family_size
                family = [
                    (event, stoic, prope(self.model))
                    for (event, stoic, prope)
                    in self.model["valid_events"][j:j + family_size]
                ]
                if len(
                    self.model["valid_events"][j + family_size:]
                ) < family_size:
                    family += self.model[
                        "valid_events"
                    ][j + family_size:]
                families.append(
                    (
                        log(1.0 / random()) / sum(
                            f[-1] for f in family
                        ),
                        family
                    )
                )
                family = min(families)
                self.model["time"].append(
                    self.model["time"][-1] + family[0]
                )
                partition = sum(f[-1] for f in family[-1]) * random()
                while partition >= 0.0:
                    event, stoic, prope = weights.pop()
                    partition -= propen
                for species, delta in stoic.items():
                    self.model[species].append(
                        self.model[species][-1] + delta
                    )
            self.model.update(event)


class FirstReaction(Base):
    """first-reaction stochastic simulation algorithm"""
    
    def method(self):
        """implementation"""
        while not self.model.equilibriated():
            times = [
                (log(1.0 / random()) / pro(self.model), sto)
                for (eve, sto, pro) in self.model["valid_events"]
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
