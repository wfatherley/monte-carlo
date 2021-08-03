"""gillespy"""
from configparser import ConfigParser
from os import path


class GillespyException(Exception):
    """something's not right"""
    pass


config_seed = None
zero_propensity = 10**-15 # TODO: model/ssa configs