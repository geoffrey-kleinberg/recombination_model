import time

import emAlgorithm
import seqGenerator
import numpy as np
import scipy.stats
import os
import shutil
import simulator
import dataInterpreter


def run_one_sim(l, q, pi, folder_name, n=100, runs=50, override=False, write_to=''):
    print("generating data")
    if override:
        shutil.rmtree(f'simData/{folder_name}')
    # makes the data files
    if not os.path.isdir(f'simData/{folder_name}'):
        os.makedirs(f'simData/{folder_name}')
        seqGenerator.main(l, q, pi, n, runs=runs, file_name_info=folder_name)

    print("running EM algorithm")
    t = time.perf_counter()
    # does the EM algorithm on all of them (records results in array)
    calculated_pis = np.zeros([1, 2 ** l])
    for run in range(runs):
        file = f'simData/{folder_name}/s{run + 1}.txt'
        calculated_pis = np.vstack((calculated_pis, emAlgorithm.main(l, q, file)[0]))
        if run % 10 == 9:
            print(f'completed run {run + 1}/50')

    run_time = time.perf_counter() - t

    calculated_pis = calculated_pis[1:]

    # calculates bias for each
    biases = np.average(calculated_pis - pi, axis=0)
    # calculates standard error
    standard_errors = scipy.stats.sem(calculated_pis)
    # calculates MSE
    mse = biases ** 2 + standard_errors ** 2

    if write_to != '':
        with open(write_to, 'a') as file:
            # file.write(f'L = {l}, q = {q}, pi = {pi}\n')
            tab = '\t'
            file.write(f'{tab.join([str(i) for i in biases])}\t')
            file.write(f'{tab.join([str(i) for i in standard_errors])}\t')
            file.write(f'{tab.join([str(i) for i in mse])}\t')
            file.write(f'{run_time}\n')
    else:
        print("Biases: " + str(biases))
        print("Standard error: " + str(standard_errors))
        print("MSE: " + str(mse))


def run_all():
    # iterate over all l
    test_l = [2, 3, 4, 5]
    # then over all q
    test_q = [0.001, 0.01, 0.05, 0.1]
    # then over all pi
    test_pi = {
        2: [[0.25, 0.25, 0.25, 0.25], [0.1, 0.2, 0.3, 0.4], [0, 0.1, 0.2, 0.7]],
        3: [
            [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125],
            [0.04, 0.06, 0.08, 0.12, 0.14, 0.16, 0.18, 0.22],
            [0, 0.02, 0.03, 0.05, 0.08, 0.12, 0.3, 0.4]
        ],
        4: [
            [0.0625 for _ in range(16)],
            [0.01, 0.02, 0.03, 0.04, 0.05, 0.03, 0.05, 0.07, 0.06, 0.08, 0.07, 0.09, 0.08, 0.1, 0.1, 0.12],
            [0, 0.02, 0.03, 0.05, 0, 0.06, 0.05, 0.09, 0, 0.08, 0.07, 0, 0.1, 0.1, 0.15, 0.2]
        ],
        5: [
            [0.03125 for _ in range(32)],
            [0.004, 0.006, 0.008, 0.012, 0.014, 0.016, 0.018, 0.022, 0.023, 0.024, 0.013, 0.017, 0.026, 0.027, 0.032,
             0.038, 0.028, 0.032, 0.036, 0.044, 0.033, 0.037, 0.042, 0.048, 0.038, 0.042, 0.046, 0.048, 0.052, 0.054,
             0.056, 0.064],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0.03, 0, 0.05, 0, 0.07, 0, 0.06, 0, 0.08, 0, 0.07, 0, 0.09, 0.01, 0.08,
             0.02, 0.1, 0.03, 0.1, 0.04, 0.12]
        ]
    }

    for i in range(len(test_l)):
        l = test_l[i]
        for j in range(len(test_q)):
            for k in range(len(test_pi[l])):
                folder_name = f'l{l}q{j + 1}pi{k + 1}'
                print(f'L = {l}, q = {test_q[j]}, pi = {test_pi[l][k]}')
                run_one_sim(test_l[i], test_q[j], test_pi[l][k], folder_name, write_to=f'simResults/out{l}EM.txt')
                print()
                print()


if __name__ == '__main__':
    #t = time.perf_counter()
    simulator.run_all_sims()
    dataInterpreter.format_all()
    #print(time.perf_counter() - t)

    # set these parameters before running
    #run_one_sim(3, 0.05, [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125], 'test', override=False)

