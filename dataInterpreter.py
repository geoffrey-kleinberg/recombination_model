import numpy as np
import scipy.stats
import os


# turns one set of raw data into formatted data
def read_one_set(folder_name, file_name, pis, true_pi):
    data = []
    with open(f'{folder_name}/{file_name}.txt', 'r') as f:
        for line in f:
            data.append([float(i) for i in line.split('\t')])
    # gets the estimates from each run of simulation
    # other is data such as time, num iterations, etc
    estimates, other = np.hsplit(np.array(data), np.arange(pis, pis + 1))
    # take the average of each
    averages = np.average(estimates, axis=0)
    other_averages = np.atleast_2d(np.average(other, axis=0)).T
    # calculate bias, standard error, and MSE
    biases = averages - true_pi
    se = scipy.stats.sem(estimates)
    mse = biases ** 2 + se ** 2

    return averages, mse, other_averages


# compiles all raw simulation results into formatted output for one set of samples
def format_file(file_name, seq_len, true_pi):
    pis = 2 ** seq_len
    max_m = seq_len - 2
    if max_m > 3:
        max_m = 3
        print("max_m is greater than 3")
    # max_m = 3
    # max_m = 0

    estimator_data = [np.array(true_pi)]
    single_data = []

    # organizes the EM algorithm data
    em_averages, em_mse, em_other_averages = read_one_set('rawResultsEM', file_name, pis, true_pi)
    estimator_data += [em_averages, em_mse]
    single_data += [em_other_averages]

    # organizes the hierarchical estimator data
    for i in range(1, max_m + 1):
        folder = f'rawResultsMCCLm{i}'
        averages, mse, other_averages = read_one_set(folder, file_name, pis, true_pi)
        estimator_data = estimator_data[:(i + 1)] + [averages] + estimator_data[(i + 1):]
        estimator_data += [mse]
        single_data += [other_averages]

    # adds the non-estimate data to arrays
    other_data = np.vstack(single_data)
    other_data = np.hstack([other_data, np.zeros((3 + max_m, 2 * (max_m + 1)))])

    out = np.hstack([np.atleast_2d(i).T for i in estimator_data])

    out = np.vstack([out, other_data])

    if not os.path.isdir(f'simResults'):
        os.mkdir(f'simResults')

    # writes all data to file
    with open(f'simResults/{file_name}.txt', 'w') as f:
        f.write('')

    with open(f'simResults/{file_name}.txt', 'a') as f:
        new_line = '\n'
        tab = '\t'
        f.write(new_line.join([tab.join([str(j) for j in i]) for i in out]))
        f.write(new_line)


# formats data for all simulation results
def format_all(test_l, test_q, test_pi):
    print('formatting all data')
    # calls format_file() for each choice of l, q, and pi
    for seq_len in test_l:
        for j in range(len(test_q)):
            for k in range(len(test_pi[seq_len])):
                file_name = f'l{seq_len}q{j + 1}pi{k + 1}'
                print(file_name)
                format_file(file_name, seq_len, test_pi[seq_len][k])


if __name__ == '__main__':
    format_all()
