"""model, and associated objects"""
from io import IOBase
from json import load as json_ld, loads as json_lds
from logging import getLogger
from math import inf
from re import compile, ASCII

from . import zero_propensity, GillespyException


LOGGER = getLogger(__name__)


class Base(dict):
    """base SSA model"""

    species_re = compile(r"(^[a-zA-Z]{1}\w{,19})", flags=ASCII)

    def equilibriated(self):
        """return True if simulation over else False"""
        raise NotImplemented

    def loader(self, obj, data, serializer):
        """load/s model JSON object"""
        if obj == "propensity":
            def process_data(data):
                """return processed propensity object"""
                dependancy_graph = dict()
                for k,v in data.items():
                    dependancy_graph[k] = [
                        s.group(0) for s
                        in self.species_re.findall()
                    ]
                    data[k] = eval(
                        "lambda d: " + self.species_re.sub(
                            lambda mo: "d['" + mo.group(0) + "']",
                            v,
                            flags=ASCII
                        )
                    )
                for event, event_species in dependancy_graph.items():
                    event_dependancies = list()
                    for other_event in dependancy_graph.keys():
                        if set(v).isdisjoint(
                            set(dependancy_graph[other_event])
                        ):
                            continue
                        event_dependancies.append(other_event)
                    dependancy_graph[event] = event_dependancies
                self.dependancy_graph = dependancy_graph
                return data
            setattr(
                self, obj, serializer(data, object_hook=process_data)
            )
        elif obj == "state":
            def process_data(data):
                """return processed state object"""
                for v in data.values():
                    if not isinstance(v, list):
                        raise GillespieException(
                            "state values must be lists"
                        )
                    if len(v) < 1:
                        raise GillespieException(
                            "state must have initial conditions"
                        )
                if "time" not in data:
                    raise GillespieException(
                        "state must have a time parameter"
                    )
                return data
            for k,v in serializer(
                data, object_hook=process_data
            ).items():
                self[k] = v
        elif obj == "stoichiometry":
            process_data = lambda data: data
            setattr(
                self, obj, serializer(data, object_hook=process_data)
            )
        else:
            LOGGER.exception("bad JSON object: obj=%s", obj)
            raise GillespyException("bad JSON object")
        LOGGER.info("validated and loaded: obj=%s", obj)

    def load(self, obj, fp):
        """load model JSON object from file pointer"""
        self.loader(obj, fp, json_ld)

    def loads(self, obj, value):
        """load model JSON object on from string"""
        self.loader(obj, fp, json_lds)


class SimpleEvents:
    """model events descriptor for Gillespian SSAs"""

    def __get__(self, instance, owner):
        """return valid events"""
        if not (
            hasattr(instance, "propensity"),
            or hasattr(instance, "state"),
            or hasattr(instance, "stoichiometry")
        ):
            raise GillespyException("missing model components")
        # ...
        

class Model(Base):
    """model for Gillespian SSAs"""

    events = SimpleEvents()

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
        super().__init__()
        for k,v in {
            "state": kwargs.get("state"),
            "propensity": kwargs.get("propensity"),
            "stoichiometry": kwargs.get("stoichiometry")
        }.items():
            if isinstance(v, dict):
                self.loader(k, v, lambda d: d)
            elif isinstance(v, str):
                self.loader(k, v, json_lds)
            elif isinstance(v, IOBase):
                self.loader(k, v, json_ld)
        self.max_duration = kwargs.get("max_duration", inf)
        self.max_steps = kwargs.get("max_steps", inf)

    def equilibriated(self):
        """return True if simulation over else False"""
        if self["time"] >= self.max_duration:
            LOGGER.info("exit on max duration")
            return True
        elif len(self["time"]) >= self.max_steps:
            LOGGER.info("exit on max steps")
            return True
        elif len(self.events) == 0:
            LOGGER.info("exit on zero cumulative propensity")
            return True
        return False


class EfficientEvents:
    """model events for Gibsonian SSAs"""
    pass


class EfficientModel:
    """model for Gibsonian SSAs"""
    pass
