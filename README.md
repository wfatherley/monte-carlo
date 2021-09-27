# gillespy
Python 3 implementations of the [stochastic simulation algorithm (SSA)](https://en.wikipedia.org/wiki/Stochastic_simulation), with convenient JSON interface.

## Install
WIP

## Usage
To perform stochastic simulations, pass model objects and desired parameters to `gillespy.Model`:

```python
import gillespy


# the "state" model object
state = {"a": [5], "b": [3], "time": [0.0]}

# the "propensity" model object
propensity = {
    "a_to_b": "3.0 * a",
    "b_to_a": "0.5 * b"
}

# the "stoichiometry" model object
stoichiometry = {
    "a_to_b": {"a": -1, "b": 1},
    "b_to_a": {"a": 1, "b": -1}
}


my_model = gillespy.Model(
    state=state,
    propensity=propensity,
    stoichiometry=stoichiometry,
    steps=5
)
```

Then, pass this model instance to one of the simulation algorithm iterators:

```python
# simulate with the "direct method"
for simulation in gillespy.Direct(my_model):
    print(simulation)
    # prints, e.g.:
    # {
    #   "a": [5, 6, 5, 5, 4, 3],
    #   "b": [3, 2, 3, 3, 4, 5],
    #   "time": [0.0, 0.7, 1.9, 2.0, 2.4, 2.7]
    # }

# simulate with the "first-family method"
simulator = 
```

For detailed specification and usage information, [read the docs]().

## References and further reading
JM Hammersley, DC Handscomb (1964). _Monte Carlo Methods_. Methuen and Company LTD.
