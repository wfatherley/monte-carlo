from mha import MHA, MHAModel
from ssa import SSA, SSAModel


class MCException(Exception):
    """Simple monte-carlo exception"""
    pass

class MCWarning(Warning):
    """Simple monte-carlo warning"""
    pass

# middle square process; mersenne twister

class Mersenne:
    """Pseudorandom number generating class"""
    def __init__(self, seed=1234):
        """Initialize pseudorandom number generator"""
        self.seed = seed
        self.j = 2**31 - 1
        self.k = 16807
        self.count = 2**30
    def floating(self, interval=None, count=1):
        """Return a pseudorandom float"""
        results = []
        for i in range(count):
            self.seed = (self.k * self.seed) % self.j
            if interval is not None:
                results.append(
                    (
                        interval[1] - interval[0]
                    ) * (self.seed / self.j) + interval[0]
                )
            else: 
                results.append(self.seed / self.j)
        return results
    
    def integer(self, interval, count=1):
        """return a pseudorandom integer"""
        pass

class MiddleSqaure:
    pass