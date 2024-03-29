"""model and associated objects"""
from logging import getLogger, ERROR
from math import inf
from secrets import choice
from string import ascii_letters, digits
from re import compile, ASCII


logger = getLogger(__name__)


species_re = compile(r"([a-zA-Z]{1}\w{,31})", flags=ASCII)


class GillespianModel(dict):
    """gillespian model"""

    dependency_map = dict()
    equilibrium_hooks = list()
    invalid_events = dict()
    valid_events = dict()
    
    def __getitem__(self, key):
        """return `key` state entity"""
        try:
            if key == "sojourn":
                return [self["time"][0]] + [
                    self["time"][k] - self["time"][k-1]
                    for k in range(1, len(self["time"]))
                ]
            return super().__getitem__(key)
        except KeyError:
            logger.error(
                "bad or missing key: model_id=%s, key=%s",
                self.id,
                key
            )
            raise
        except IndexError:
            return self["time"]

    def __init__(self, **kwargs):
        """construct self"""
        self.id = kwargs.pop("id", None) or "".join(
            choice(ascii_letters + digits) for _ in range(32)
        )
        self.duration = kwargs.pop("duration", inf)
        self.steps = 0
        self.max_steps = kwargs.pop("max_steps", inf)
        self.equilibrium_hooks.extend(
            kwargs.pop("equilibrium_hooks", [])
        )
        propensity = kwargs.pop("propensity", {})
        stoichiometry = kwargs.pop("stoichiometry", {})
        kwargs.update(kwargs.pop("state", {}))
        if "time" not in kwargs:
            raise Exception("state object requires `time` series")
        super().__init__(**kwargs)
        self.build_events(
            propensity=propensity, stoichiometry=stoichiometry
        )
        self.build_dependency_map()

    def build_dependency_map(self):
        """set dependency map attribute"""
        dep_map = dict()
        for event, objects in [
            (tup[0], tup) for tup in (
                list(self.invalid_events.values())
                + list(self.valid_events.values())
            )
        ]:
            dep_map[event] = [
                species for species, delta
                in objects[1].items() if abs(delta) > 0
            ]
        for event, event_species in dep_map.items():
            event_deps = list()
            for other_event in dep_map.keys():
                if not set(event_species).isdisjoint(
                    set(dep_map[other_event])
                ):
                    event_deps.append(other_event)
            self.dependency_map[event] = event_deps

    def build_events(self, propensity={}, stoichiometry={}):
        """add or modify model events"""
        if set(propensity.keys()) != set(stoichiometry.keys()):
            logger.error(
                "mismatched event names: model_id=%s", self.id
            )
            raise Exception("bad event objects")
        for eve, sto, pro in tuple(
            (eve, stoichiometry[eve], propensity[eve])
            for eve in stoichiometry.keys()
        ):
            pro = self.build_propensity(pro)
            if pro(self) > 10**-15:
                self.valid_events[eve] = (eve,sto,pro)
            else:
                self.invalid_events[eve] = (eve,sto,pro)

    def build_propensity(self, propensity):
        """return anonymous, evaluable propensity"""
        if isinstance(propensity, str):
            return eval(
                "lambda d: " + species_re.sub(
                    lambda mo: "d['" + mo.group(0) + "'][-1]",
                    propensity
                )
            )
        return propensity

    def equilibriated(self):
        """return True if simulation over else False"""
        if self["time"][-1] >= self.duration:
            logger.info("exit on duration: model_id=%s", self.id)
            return True
        if self.steps == self.max_steps:
            logger.info("exit on steps: model_id=%s", self.id)
            return True
        if len(self.valid_events) == 0:
            return True
        if any(h() for h in self.equilibrium_hooks):
            return True
        return False

    def equilibrium_hook(self, func):
        """decorate equilibrium hook"""
        def wrapped_hook():
            """return True if `func` detects equilibrium"""
            if func(self) is True:
                func_doc = func.__doc__.strip() or func.__name__
                logger.info(func_doc + ": model_id=%s", self.id)
                return True
            return False
        self.equilibrium_hooks.append(wrapped_hook)

    def initialize(self):
        """initialize self"""
        if self.steps > 0:
            for k in self:
                del self[k][1:]
            for eve, sto, pro in (
                list(self.valid_events.values())
                + list(self.invalid_events.values())
            ):
                if pro(self) > 10**-15:
                    self.valid_events[eve] = (eve,sto,pro)
                else:
                    self.invalid_events[eve] = (eve,sto,pro)
            self.steps = 0

    def update(self, event, stoichiometry, sojourn):
        """update self given event"""
        self.steps += 1
        self.update_time(sojourn)
        self.update_species(stoichiometry)
        self.update_events(event)

    def update_events(self, event):
        """update model events"""
        for dep_event in self.dependency_map[event]:
            try:
                pro = self.valid_events[dep_event][2](self)
                if pro <= 10**-15:
                    eve = self.valid_events.pop(dep_event)
                    self.invalid_events[dep_event] = eve
            except KeyError:
                pro = self.invalid_events[dep_event][2](self)
                if pro > 10**-15:
                    eve = self.invalid_events.pop(dep_event)
                    self.valid_events[dep_event] = eve

    def update_species(self, stoichiometry):
        """update model species"""
        for species, delta in stoichiometry.items():
            self[species].append(self[species][-1] + delta)

    def update_time(self, sojourn):
        """update model time"""
        self["time"].append(self["time"][-1] + sojourn)


__all__ = ["GillespianModel", "species_re",]
