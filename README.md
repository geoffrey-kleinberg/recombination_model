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

To use the EM algorithm: in `emAlgorithm.py`, run `main()` with 
the following parameters:
- `L` is an integer with the length of the DNA sequences
- `q` is a float representing the probability of recombination
- `file_name` is a path to the file containing the data

To use the Hierarchical estimator: In `simulator.R`, run 
`get_and_write_estimates` with the following parameters:
- `L` is an integer with the length of the DNA sequences
- `q_val` is a float representing the probability of recombination
- `m` is an integer representing the margin to use, which must be 
between 1 and `L` - 2, and can't be greater than 3.
- `n` is the number of DNA sequences in the data set
- `file_name` is a path to the file containing the data

### To run a simulation

In `main.py`, set the parameters `test_l`, `test_q`, `test_pi`, and
`test_pi`. Then, call `run_full_simulation()` to do the simulation. 
- `test_l` is an array of integers representing the lengths you
want to test
- `test_q` is an array of floats representing the recombination
probabilities you want to test
- `test_pi` is a dictionary with a key for each element, *L* of 
`test_l`. The value associated with each key is an array of length
2^*L* representing the true ancestral distribution.
- `test_m` is a dictionary with a key for each element, *L* of 
`test_l`. The value associated with each key is an array with integers
ranging from 1 to *L -* 2, with a maximum of 4.
- `generate_new_data` is a Boolean. It is False by default, but can be 
True if you don't have existing data or want to generate new data.