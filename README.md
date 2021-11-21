# Critical Behavior of Coupled Scalar Fields with Different Dynamics

---

This project was created for my bachelor's thesis.
The goal of this thesis is to investigate properties of two coupled fields near the critical point.


## Compiliation and Execution

A recent version of the nvcc and the appropriate nvdida driver is needed.

```shell
./cmake .
./make coupledDynamics

# Execution
./coupledDynamics --help
```


## Output

The output is stored in a JSON-like Yaml-file. The default name is measurements.yaml.
Every simulation is stored as a json-object in a single line.
Therefore, the same file can be used for multiple simulations.
Filtering of the results can be done with the database.py script.

The json-object contains every parameter as well as every measurement.
Schema of the measurement-object:

```javascript
{
  "N": Integer, //lattice-sites
  "T": Float, //temperature
  "C": Float, //coupling constant
  "D": Float, //dimension
  "t": Float, //step-width
  "msq_A": Float,
  "lambda_A": Float,
  "gamma_A": Float,
  "msq_B": Float,
  "gamma_B": Float,
  "J": Float, //external field
  "mu": Float,
  "fastTherm": Boolean,
  "thermTime1": Integer,
  "thermTime2": Integer,
  "thermTime3": Integer,
  "dynamicChangeMode": ENUM("ReInit","Zero","Subtract"),
  "measure_time": Integer,
  "p": Boolean, //print every step and output to additional file
  "seed": Integer, 
  "A": Integer, // amount of simulations
  "S": Float, // external thermalization field
  // mesaurement results:
  "sigma": [Float],
  "n": [Float],
  "sigma_sq": [Float],
  "n_sq": [Float]
},
```

If the observables are calculated for every step, the measurement also contains measurements for the momentum fields.

## Contributing
By creating a pull request you accept that the code is licenced under the conditions specified in the LICENCE.txt file.