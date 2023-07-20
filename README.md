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

### To run a simulation

In `main.py`, set the parameters `test_l`, `test_q`, and `test_pi`.
Generate the data sets with `simulator.generate_all()`. By default,
this will generate 50 data sets of 100 DNA sequences for each choice
of l, q, and pi and save all the data in a `simData` directory. Then,
run `simulator.run_all_sims()` to run the EM algorithm on each
data set. Then, run the entire `simulator.R` file to run the
hierarchical estimator on each data set. Finally, run
`dataInterpreter.format_all()` to get a summary of the results. 
This will output the parameter estimates, the MSE, and the computation
time for each method.