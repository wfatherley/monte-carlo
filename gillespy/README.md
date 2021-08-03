# gillespy
Python 3 implementations of the [stochastic simulation algorithm](https://en.wikipedia.org/wiki/Stochastic_simulation), with convenient JSON interface.

## Install
WIP

## Usage
To perform stochastic simulations, `gillespy` exposes varieties of two Python objects,

 - __model__ classes for storing an managing initial and final simulation data;
 - and iterable __SSA__ classes containing the simulation algorithm.

The basic flow is to instantiate one of the model classes, instantiate one of the SSA classes by passing it the model instance, and either use the SSA instance as the iterable in a for-loop or in some other control flow structure.

Instantiating a model requires passing it three JSON objects, including a __state__ object,

```json
{
    "species_a": [4],
    "species_b": [13],
    "time": [0.0]
}
```

whereby each entity is an array with a single value representing that entitys initial condition. The `species_a` entity above can be interpreted as _the initial population of species A is 4._ Although its initial condition can be any numerical value, the state object must always have a `time` entity.

In addition to passing in these arrays of initial conditions via the `state` object, there are the __propensity__ and __stoichiometry__ objects, which specify entities that define _elementary events_ that can happen at each step within the simulation. Suppose the process that is to be simulated with the state object above is the composition of only two elementary events--

 - one of species A transitions to being a B species,
 - and vice versa.

Then the propensity object

```json
{
    "a_2_b": "0.5 * species_a",
    "b_2_a": "1.25 * species_b"
}
```

defines propensity expressions for each of the two possible elementary events. Each propensity expression represents the instructions need to compute the frequency/likelihood of that event, which is required by the SSA object to perform the simulation. During each iteration of a given simulation, the arrays in the state object are updated when the SSA selects an elementary event based on their propensities, and this update occurs via a stoichiometry object that specifies exactly how to update populations--

```json
{
    "a_2_b": {"species_a": -1, "species_b": 1},
    "b_2_a": {"species_a": 1, "species_b": -1}
}
```

There are several ways to pass these three objects to a model class, which are shown below with `gillespy.Model`.

```python
from json import import loads

import gillespy


# read in the JSON files
with (
    open("my_state.json", "r")) as state,
    open("my_propensity.json", "r") as prope,
    open("my_stoichiometry.json", "r") as stoic
):
    # pass in file pointers
    my_model = gillespy.Model(
        state=state, propensity=prope, stoichometry=stoic
    )

    state, prope, stoic = (
        state.read(), prope.read(), stoic.read()
    )

    # or pass strings
    my_duplicate_model = gillespy.Model(
        state=state, propensity=prope, stoichometry=stoic
    )

    # or pass serialized objects
    my_doubly_duplicate_model = gillespy.Model(
        state=loads(state),
        propensity=loads(prope),
        stoichometry=loads(stoic)
    )

    # or instantiate the model without the objects
    another_duplicate_model = gillespy.Model()

    # and read the strings (or file pointers via `Model.load`) in later
    another_duplicate_model.loads("state", state)
    another_duplicate_model.loads("propensity", prope)
    another_duplicate_model.loads("stoichiometry", stoic)
```

See the documentation for more information on selecting and instantiating models.

With an instantiated model, an __SSA__ class can be instantaited, which is an iterator that returns some number of complete simulations, also known as __trajectories__.

```python
# pass model to an SSA class
ssa = gillespy.DirectSSA(my_model)

# get a trajectory (i.e. one simulation)
trajectory = next(ssa)

# or get 100 trajectories
trajectory_count = 100
while trajectory_count > 0:
    trajectory_count -= 1
    print(next(ssa))

# or get inifinity trajectories
for trajectory in ssa:
    print(trajectory)
```

As can be seen, the default behvior of SSA instances is to produce and return an indefinite number of trajectories. See the documentation for more information on selecting and instantiating SSAs, and tuning their behavior.

## References and further reading
JM Hammersley, DC Handscomb (1964). _Monte Carlo Methods_. Methuen and Company LTD.
