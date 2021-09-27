"""model and associated objects"""
from datetime import datetime, timezone
from logging import getLogger, ERROR
from math import inf
from secrets import choice
from string import ascii_letters, digits

from .util import (
    default_seed, species_re, zero_propensity, GillespyException
)


__all__ = ["Model",]


logger = getLogger(__name__)
logger.setLevel(ERROR)


class Model(dict):
    """model for stochastic simulation algorithm"""

    equilibrium_hooks = list()
    invalid_events = dict()
    valid_events = dict()
    
    def __getitem__(self, key):
        """return `key` state entity"""
        if key == "soujorn":
            return self["time"][0] + [
                self["time"][k] - self["time"][k-1]
                for k in range(1, len(self["time"]))
            ]
        try:
            return super().__getitem__(key)
        except KeyError:
            logger.error(
                "bad key passed: model_id=%s, key=%s", self.id, key
            )
            raise

    def __init__(self, **kwargs):
        """construct self"""
        self.id = "".join(
            choice(ascii_letters + digits) for _ in range(32)
        )
        self.duration = kwargs.pop("duration", inf)
        self.steps = 0
        self.max_steps = kwargs.pop("max_steps", inf)
        self.equilibrium_hooks.extend(
            kwargs.pop("equilibrium_hooks", [])
        )
        self.seed = (
            kwargs.pop("seed", None)
            or int(datetime.now(timezone.utc).timestamp())
        )
        propensity = kwargs.pop("propensity", {})
        stoichiometry = kwargs.pop("stoichiometry", {})
        kwargs.update(kwargs.pop("state", {}))
        super().__init__(**kwargs)
        self.build_events(
            propensity=propensity, stoichiometry=stoichiometry
        )

    def __setitem__(self, key, value):
        """vaildate and set `key` state entity"""
        raise Exception
        if not isinstance(value, list):
            logger.error(
                "state values must be lists: model_id=%s", self.id
            )
            raise GillespieException("invalid state object")
        if len(value) < 1:
            logger.error(
                "state value list must not be empty: model_id=%s",
                self.id
            )
            raise GillespieException("invalid state object")
        if not isinstance(value[-1], (int, complex, float)):
            logger.error(
                "non-numeric state value entry: model_id=%s", self.id
            )
            raise GillespieException("invalid state object")
        if not species_re.fullmatch(key):
            logger.error(
                "invalid state key: model_id=%s", self.id
            )
            raise GillespieException("invalid state object")
        super().__setitem__(key, value)

    def build_dependency_graph(self):
        """set dependency graph attribute"""
        dep_graph = dict()
        for event, objects in [
            (tup[0], tup) for tup in (
                list(self.invalid_events.values())
                + list(self.valid_events.values())
            )
        ]:
            dep_graph[event] = [
                species for species, delta
                in objects[1].items() if delta > 0
            ]
        for event, event_species in dep_graph.items():
            event_deps = list()
            for other_event in dep_graph.keys():
                if not set(event_species).isdisjoint(
                    set(dep_graph[other_event])
                ):
                    event_deps.append(other_event)
            dep_graph[event] = event_deps
        self.dep_graph = dep_graph

    def build_propensity_lambda(self, propensity):
        """return anonymous, evaluable propensity"""
        return eval(
            "lambda d: " + species_re.sub(
                lambda mo: "d['" + mo.group(0) + "'][-1]", propensity
            )
        )

    def build_events(self, propensity={}, stoichiometry={}):
        """add or modify model events"""
        if set(propensity.keys()) != set(stoichiometry.keys()):
            logger.error(
                "mismatched event names: model_id=%s", self.id
            )
            raise GillespyException("bad event objects")
        for eve, sto, pro in tuple(
            (eve, stoichiometry[eve], propensity[eve])
            for eve in stoichiometry.keys()
        ):
            pro = self.build_propensity_lambda(pro)
            if pro(self) > zero_propensity:
                self.valid_events[eve] = (eve,sto,pro)
            else:
                self.invalid_events[eve] = (eve,sto,pro)
        self.build_dependency_graph()

    def equilibriated(self):
        """return True if simulation over else False"""
        if self["time"][-1] >= self.duration:
            logger.info("exit on duration: model_id=%s", self.id)
            for key in self:
                self[key].pop()
            return True
        elif self.steps == self.max_steps:
            logger.info("exit on steps: model_id=%s", self.id)
            return True
        elif not any(self.valid_events):
            return True
        elif any(h() for h in self.equilibrium_hooks):
            return True
        return False

    def equilibrium_hook(self, func):
        """decorate equilibrium hook"""
        def wrapped_hook():
            """return True if `func` detects equilibrium"""
            if func(self) is True:
                logger.info(
                    func.__doc__.strip() + ": model_id=%s", self.id
                )
                return True
            return False
        self.equilibrium_hooks.append(wrapped_hook)

    def reset(self):
        """reset self"""
        self.steps = 0
        for k in self:
            del self[k][1:]

    def update(self, event):
        """update model given event"""
        self.steps += 1
        for dep_event in self.dep_graph[event]:
            try:
                pro = self.valid_events[dep_event][2](self)
                if pro <= zero_propensity:
                    eve = self.valid_events.pop(dep_event)
                    self.invalid_events[dep_event] = eve
            except KeyError:
                pro = self.invalid_events[dep_event][2](self)
                if pro > zero_propensity:
                    eve = self.invalid_events.pop(dep_event)
                    self.valid_events[dep_event] = eve
