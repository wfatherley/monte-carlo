from mha import MHA, MHAModel
from ssa import SSA, SSAModel


class MCException(Exception):
    """Simple monte-carlo exception"""
    pass

# middle square process; mersenne twister