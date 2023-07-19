# Comparison of EM Algorithm to Hierarchical Estimator

### EM Algorithm

The EM Algorithm is an iterative method to calculate the MLE
of observed data when there is a latent variable. This
algorithm is implemented in Python in `emAlgorithm.py`.

### Hierarchical Estimator

The Hierarchical Estimator uses Markov Chain Composite
Likelihood (MCCL) to approximate the MLE by assuming the
observed data follow a Markov Chain. This algorithm is
implemented in R. See `DOCUMENTATION.txt` for details.


## How To Use

### To run 1 data set

Create a file of the binarized DNA sequences, with one observed
sequence per line. 

To use the EM algorithm: in `emAlgorithm.py`, run main with parameters `L` 
and `q` representing the length and probability of recombination,
respectively. The third argument is the path to the file with the
data.

To use the Hierarchical estimator: *will be added later*

### To generate data from a known pi

Use `main()` in  `seqGenerator.py` with a given `L`, `q`, `true_pi`,
and `n`. Can also include a number of data sets to generate (`runs`)
and a file path (`file_name_info`) within `simData` to save to.
Data will be saved to `simData/{file_name_info}/s{num}.txt` for use
in simulation.

### To run a simulation

Create all data files with the format `simData/{label}/s{num}.txt`,
where `label` groups data with common L and q values and `num`
goes from 1 to the number of data sets. Adjust the `L` and `q` in
`simulator.py` and `simulator.R` as needed, then run both.

Use `format_all()` in `dataInterpreter.py` to get a read-out of the
results from the simulation. Be sure to set the `L`, `q`, and `pi`
values that the data were generated with.