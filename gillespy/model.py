"""model and associated objects"""
from io import IOBase
from json import load as json_ld, loads as json_lds
from logging import getLogger
from math import inf
from queue import PriorityQueue
from re import compile, ASCII
from secrets import choice
from string import ascii_letters, digits

from .util import zero_propensity, GillespyException


__all__ = ["BaseModel", "Model", "EfficientModel"]


logger = getLogger(__name__)


species_re = compile(r"(^[a-zA-Z]{1}\w{,31})", flags=ASCII)


def process_data(data):
    """return processed propensity object"""
    dep_graph = dict()
    for k,v in data.items():
        dep_graph[k] = [
            s.group(0) for s
            in species_re.findall()
        ]
        data[k] = eval(
            "lambda d: " + species_re.sub(
                lambda mo: "d['" + mo.group(0) + "']",
                v,
                flags=ASCII
            )
        )
    for event, event_species in dep_graph.items():
        event_deps = list()
        for other_event in dep_graph.keys():
            if set(v).isdisjoint(
                set(dep_graph[other_event])
            ):
                continue
            event_deps.append(other_event)
        dep_graph[event] = event_deps
    return data


class Propensity:
    """immutable propensity descriptor"""

    def __delete__(self, obj):
        """delete propensities"""
        pass

    def __get__(self, obj, obj_owner=None):
        """return propensities"""
        pass

    def __set__(self, obj, value):
        """set propensities"""
        pass

    def __set_name__(self, obj_owner, name):
        """remember attribute name"""
        pass


class Stoichiometry:
    """immutable stoichiometry descriptor"""

    def __delete__(self, obj):
        """delete stoichiometries"""
        pass

    def __get__(self, obj, obj_owner=None):
        """return stoichiometries"""
        pass

    def __set__(self, obj, value):
        """set stoichiometries"""
        pass

    def __set_name__(self, obj_owner, name):
        """remember attribute name"""
        pass


class Model(dict):
    """model container for Gillespian SSAs"""

    propensity = Propensity()
    stoichiometry = Stoichiometry()
    
    def __getitem__(self, key):
        """get series data"""
        if key in self:
            return super().__getitem__(key)
        elif key == "soujorn":
            return self["time"][0] + [
                self["time"][k] - self["time"][k-1]
                for k in range(1, len(self["time"]))
            ]
        raise KeyError("no such observable")

    def __init__(self, **kwargs):
        """construct self"""
        self.id = "".join(
            choice(ascii_letters + digits) for _ in range(32)
        )
        super().__init__()
        if "propensity" in kwargs:
            self.propensity = kwargs.pop("propensity")
        if "stoichiometry" in kwargs:
            self.stoichiometry = kwargs.pop("stoichiometry")
        if "max_duration" in kwargs:
            self.max_duration = kwargs.pop("max_duration")
        else:
            self.max_duration = inf
        if "max_steps" in kwargs:
            self.max_steps = kwargs.pop("max_steps")
        else:
            self.max_steps = inf
        if "state" in kwargs:
            self.loader("state", kwargs.pop("state"))
        for k,v in kwargs.items():
            self[k] = v

    def __setitem__(self, key, value):
        """vaildate and return state object"""
        if not isinstance(value, list):
            raise GillespieException("invalid state")
        if len(value) < 1:
            raise GillespieException("invalid state")
        if not isinstance(value[-1], (int, complex, float)):
            raise GillespieException("invalid state")
        if not species_re.fullmatch(k):
            raise GillespieException("invalid state")
        super().__setitem(key, value)

    def equilibriated(self):
        """return True if simulation over else False"""
        if self["time"] >= self.max_duration:
            logger.info("exit on max duration: model_id=%i", self.id)
            return True
        elif len(self["time"]) >= self.max_steps:
            logger.info("exit on max steps: model_id=%i", self.id)
            return True
        elif len(self.events) == 0:
            logger.info(
                "exit on zero cumulative propensity: model_id=%i",
                self.id
            )
            return True
        return False

    def loader(self, obj, data)
        """load/s model object data"""
        if obj == "state":
            if isinstance(state, str):
                state = json_lds(state)
            elif isinstance(state, IOBase):
                state = json_ld(state)
            if not isinstance(state, dict):
                raise GillespyException("invalid state")
            for k,v in state.items():
                self[k] = v
        elif obj in {"propensity", "stoichiometry"}:
            setattr(self, obj, data)
        else:
            raise GillespyException("bad object")


class EfficientModel(Model):
    """model for Gibsonian SSAs"""
    
    invalid_events = None
    sojourn_tree = PriorityQueue()

    @property
    def events(self):
        """return list of possible events"""
        pass

    def equilibriated(self):
        """return True if simulation over else False"""
        pass

    def loader(self, obj, data, serializer):
        """load/s model JSON object"""
        #super().loader(obj, data, serializer)

    def update(self, event):
        """update model, given event"""
        pass
