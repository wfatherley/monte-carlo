from warnings import warn

from mha import MHA, MHAModel
from ssa import SSA, SSAModel


class MCException(Exception):
    """Simple monte-carlo exception"""
    pass


class Mersenne:
    """Pseudorandom number generater"""

    def __init__(self, seed=1234):
        """Initialize pseudorandom number generator"""
        self.seed = seed
        self.j = 2**31 - 1
        self.k = 16807
        self.count = 2**30

    def floating(self, interval=None, count=1):
        """
        Return a pseudorandom float. Default is one floating-
        point number between zero and one. Pass in a tuple or
        list, (a,b), to return a floating-point number on
        [a,b]. If count is 1, a single number is returned,
        otherwise a list of numbers is returned.
        """
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
            self.count -= 1
            if self.count == 0:
                warn(
                    "Pseudorandom period nearing!!",
                    category=ResourceWarning
                )
        if count == 1:
            return results.pop()
        else:
            return results
    
    def integer(self, interval=None, count=1):
        """return a pseudorandom integer"""
        results = []
        for i in range(count):
            self.seed = (self.k * self.seed) % self.j
            if interval is not None:
                results.append(
                    int((
                        interval[1] - interval[0] + 1
                    ) * (self.seed / self.j) + interval[0])
                )
            else:
                result = self.seed / self.j
                if result < 0.50:
                    results.append(0)
                else:
                    results.append(1)
            self.count -= 1
            if self.count == 0:
                warn(
                    "Pseudorandom period nearing!!",
                    category=ResourceWarning
                )
        if count == 1:
            return results.pop()
        else:
            return results


class MiddleSqaure:
    pass