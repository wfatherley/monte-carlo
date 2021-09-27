"""utilities"""
from re import compile, ASCII


default_seed = None
species_re = compile(r"([a-zA-Z]{1}\w{,31})", flags=ASCII)
zero_propensity = 10**-15


class GillespyException(Exception):
    """something's not right"""
    pass
