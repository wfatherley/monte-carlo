"""gillespy unit tests"""
from random import randint
from unittest import TestCase

from ssa import Base, Direct, FirstReaction, Model


state = {
    'a': [10, 9, 10, 9, 8, 7, 8, 7, 6, 5, 6],
    'b': [0, 1, 0, 1, 2, 3, 2, 3, 4, 5, 4],
    'time': [
        0.0, 0.008992303463353597, 0.5835134815757584,
        0.6194587363159768, 0.6996347430132124,
        0.7998282296135547, 0.9083589167245347,
        1.0768417124506244, 1.1417075818688178,
        1.1520935500793725, 1.314095293338304
    ]
}
propensity = {0: "0.5 * a", 1: "0.5 * b"}
stoichiometry={0: {"a": -1, "b": 1}, 1: {"a": 1, "b": -1}}


class Base(TestCase):
    """base unit test"""

    @classmethod
    def setUpClass(cls):
        """set up class for unit tests"""
        cls.seed = 2
        cls.result = state
        cls.model = Model(
            state=cls.result.copy(),
            propensity=propensity,
            stoichiometry=stoichiometry,
            max_steps=10
        )


class TestModels(Base):
    """test model objects"""a
    
    @classmethod
    def setUpClass(cls):
        """set up class for unit tests"""
        cls.state = state
        cls.propensity = propensity
        cls.stoichiometry = stoichiometry

    def test_Model(self):
        """test Model object"""
        #
        # verify exceptions are raised
        self.assertRaises(Exception, Model)


class TestMethods(Base):
    """test method objects"""

    def test_Base(self):
        """test base method object"""
        #
        # verify exceptions are raised
        self.assertRaises(TypeError, Base)
        self.assertRaises(NotImplemented, Base(self.model).method)
        with self.assertRaises(StopIteration):
            for trajectory in Base(self.model, trajectories=0):
                continue
        #
        # verify constructor alters seed
        iterator = Base(self.model, seed=self.seed)
        self.assertTrue(randint(0, 100) == 7)

    def test_Direct(self):
        """test direct method"""
        #
        # verify method correctly produces seed-constrained results
        trajectory = next(
            Direct(self.model, seed=self.seed, trajectories=1)
        )
        self.assertTrue(trajectory.model["a"] == self.result["a"])
        self.assertTrue(trajectory.model["b"] == self.result["b"])

    def test_FirstReaction(self):
        """test first-recation method"""
        #
        # verify method correctly produces seed-constrained results
        trajectory = next(
            FirstReaction(self.model, seed=self.seed, trajectories=1)
        )
        self.assertTrue(trajectory.model["a"] == self.result["a"])
        self.assertTrue(trajectory.model["b"] == self.result["b"])
        
