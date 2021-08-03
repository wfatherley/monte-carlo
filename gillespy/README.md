# gillespy
Python 3 implementations of the [stochastic simulation algorithm](https://en.wikipedia.org/wiki/Stochastic_simulation), with convenient JSON interface.

## Usage
To perform a simulation, `gillespy` requires three JSON objects, a __state__ object,

```json
{
    "species_a": [4],
    "species_b": [13],
    "time": [0.0]
}
```

whereby each entity is array with a single value representing that entitys initial population. Note that there is also a required `time` entity, and its initial value is the simulation "start time".

In addition to passing in these arrays of initial conditions via the `state` object, there are the __propensity__ and __stoichiometry__ objects, which specify entities that define _elementary events_ that can happen at each step within the simulation.  The propensity object

```json
{
    "a_2_b": "0.5 * species_a",
    "b_2_a": "1.25 * species_b"
}
```

defines the frequency/likelihood of each of the two possible elementary evenets (some number of species A transitioning to B, and vice versa), and the stoichiometry object specifies exactly the "some number"--

```json
{
    "a_2_b": {"species_a": -1, "species_b": 1},
    "b_2_a": {"species_a": 1, "species_b": -1}
}
```

The "get-started" usage for simulating the system specified by the above three objects occurs in two steps. The first is creating a _model_ instance that handles all logic related to insuring the validity of simulation data. There are several different model classes, but they all accept the above three JSON objects in the same multiples of ways:

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

    # or instantiate the model without JSON
    another_duplicate_model = gillespy.Model()

    # and read the strings in later
    # (there is also a `load` method for filelikes)
    another_duplicate_model.loads("state", state)
    another_duplicate_model.loads("propensity", prope)
    another_duplicate_model.loads("stoichiometry", stoic)
```

See the documentation for more information on instantiating models.

With an instantiated model, an __SSA (stochastic simulation algorithm)__ class can be instantaited, which is an iterator that returns some number of complete simulation __trajectories__.


`gillespy` has several classes that can be used to perform a simulation with the above three objects, and they all differ mainly in efficiency. 

## References and further reading
JM Hammersley, DC Handscomb (1964). _Monte Carlo Methods_. Methuen and Company LTD.
