import emAlgorithm
import seqGenerator
import numpy as np
import scipy.stats
import os
import shutil


def run_one_sim(l, q, pi, folder_name, n=100, runs=50, override=False):
    print("generating data")
    if override:
        shutil.rmtree(f'simData/{folder_name}')
    # makes the data files
    if not os.path.isdir(f'simData/{folder_name}'):
        os.makedirs(f'simData/{folder_name}')
        seqGenerator.main(l, q, pi, n, runs=runs, file_name_info=folder_name)

    print("running EM algorithm")
    # does the EM algorithm on all of them (records results in array)
    calculated_pis = np.zeros([1, 2 ** l])
    for run in range(runs):
        file = f'simData/{folder_name}/s{run + 1}.txt'
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
    print("MSE: " + str(mse))


if __name__ == '__main__':
    # iterate over all l
    # test_l = [2, 3, 4, 5]
    # then over all q
    # test_q = [0.001, 0.01, 0.05, 0.1]
    # then over all pi
    # test_pi = {
    # 2: [[], [], []],
    # 3: [[], [], []], ...
    # }
    # for i in range(len(test_l)):
    #    l = test_l[i]
    #    for j in range(len(test_q)):
    #        for k in range(len(test_pi[l])):
    #            folder_name = f'l{i + 2}q{j + 1}pi{k + 1}'
    #            run_one_sim(test_l[i], test_q[j], test_pi[l][k], folder_name)


    # set these parameters before running
    run_one_sim(3, 0.05, [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125], 'test', override=True)

