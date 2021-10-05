"""method and associated objects"""
from datetime import datetime, timezone
from logging import getLogger, ERROR
from math import inf, log
from random import random, seed
from time import perf_counter


logger = getLogger(__name__)
logger.setLevel(ERROR)


class Base:
    """base stochastic simulation algorithm"""
    
    def __init__(self, model, **kwargs):
        """construct self"""
        self.model = model
        self.trajectories = kwargs.get("trajectories", inf)
        self.seed = (
            kwargs.get("seed")
            or int(datetime.now(timezone.utc).timestamp())
        )
        seed(a=self.seed)
        logger.info(
            "seed set: model_id=%s, seed_value=%s",
            model.id,
            str(self.seed)
        )

    def __iter__(self):
        """return self"""
        logger.info(
            "returning simulation iterator: model_id=%s",
            self.model.id
        )
        return self

    def __next__(self):
        """return stochastic simulation trajectory"""
        while self.trajectories > 0:
            start = perf_counter()
            self.method()
            stop = perf_counter()
            logger.info(
                "simulation complete: model_id=%s, perf_time=%f",
                self.model.id,
                stop - start
            )
            self.trajectories -= 1
            return self.model
        raise StopIteration

    def method(self):
        """stochastic simulation algorithm"""
        raise NotImplemented


class Direct(Base):
    """direct-method stochastic simulation algorithm"""
    
    def method(self):
        """implementation"""
        self.model.reset()
        while not self.model.equilibriated():
            print(self.model.dep_graph)
            weights = [
                (event, stoic, prope(self.model))
                for (event, stoic, prope)
                in self.model.valid_events.values()
            ]
            partition = sum(w[-1] for w in weights)
            sojourn = log(1.0 / random()) / partition
            self.model["time"].append(
                self.model["time"][-1] + sojourn
            )
            partition = partition * random()
            while partition >= 0.0:
                event, stoic, prope = weights.pop()
                partition -= prope
            for species, delta in stoic.items():
                self.model[species].append(
                    self.model[species][-1] + delta
                )
            self.model.update_events(event)


class FirstReaction(Base):
    """first-reaction stochastic simulation algorithm"""
    
    def method(self):
        """implementation"""
        self.model.reset()
        while not self.model.equilibriated():
            times = [
                (log(1.0 / random()) / pro(self.model), sto, eve)
                for (eve, sto, pro)
                in self.model.valid_events.values()
            ]
            times.sort()
            self.model["time"].append(
                self.model["time"][-1] + times[0][0]
            )
            for species, delta in times[0][1].items():
                self.model[species].append(
                    self.model[species][-1] + delta
                )
            self.model.update_events(times[0][-1])
