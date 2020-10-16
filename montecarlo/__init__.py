from mha import Base as BaseMHA
from ssa import Base as BaseSSA


class MCException(Exception):
    """Simple monte-carlo exception"""
    pass

# middle square process; mersenne twister