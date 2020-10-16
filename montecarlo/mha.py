"""_Metropolis-Hastings Algorithm (MHA)_

- need loss func like l1 norm
- need density/mass func that loss goes by
- if current loss < loss[-1] accept
- if current loss > loss[-1] then
    * if random() < density(current loss), accept
    * else update with loss[-1]
"""
from random import random


class Base(dict):
    """Container for MHAs"""
    def __init__(self, *args, **kwargs):
        pass