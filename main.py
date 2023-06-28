import emAlgorithm
import seqGenerator
import numpy as np
import scipy.stats


if __name__ == '__main__':
    # set these parameters before running
    l = 2
    q = 0.05
    pi = [0.25, 0.25, 0.25, 0.25]
    n = 100
    runs = 50
    file_name = "l2q005pi1"
    dataAlreadyGenerated = True

    print("generating data")
    # makes the data files
    if not dataAlreadyGenerated:
        seqGenerator.main(l, q, pi, n, runs=runs, file_name_info=file_name)

    print("running EM algorithm")
    # does the EM algorithm on all of them (records results in array)
    calculated_pis = np.zeros([1, 4])
    for run in range(runs):
        file = f'simData/{file_name}s{run + 1}.txt'
        calculated_pis = np.vstack((calculated_pis, emAlgorithm.main(l, q, file)))
        if run % 10 == 9:
            print(f'completed run {run + 1}/50')

    calculated_pis = calculated_pis[1:]
    # calculates bias for each
    biases = np.average(calculated_pis, axis=0) - pi
    print("Biases: " + str(biases))
    # calculates standard error
    standard_errors = scipy.stats.sem(calculated_pis)
    print("Standard error: " + str(standard_errors))
    # calculates MSE
    mse = biases ** 2 + standard_errors ** 2
    print("MSE: " + str(mse))24

