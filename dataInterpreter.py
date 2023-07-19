import numpy as np
import scipy.stats


def read_one_set(folder_name, file_name, pis, true_pi):
    data = []
    with open(f'{folder_name}/{file_name}.txt', 'r') as f:
        for line in f:
            data.append([float(i) for i in line.split('\t')])
    estimates, other = np.hsplit(np.array(data), np.arange(pis, pis + 1))
    averages = np.average(estimates, axis=0)
    other_averages = np.atleast_2d(np.average(other, axis=0)).T
    biases = averages - true_pi
    se = scipy.stats.sem(estimates)
    mse = biases ** 2 + se ** 2

    return averages, mse, other_averages


def format_file(file_name, seq_len, true_pi):
    pis = 2 ** seq_len
    max_m = seq_len - 2

    estimator_data = [np.array(true_pi)]
    single_data = []

    em_averages, em_mse, em_other_averages = read_one_set('rawResultsEM', file_name, pis, true_pi)
    estimator_data += [em_averages, em_mse]
    single_data += [em_other_averages]

    for i in range(1, max_m + 1):
        folder = f'rawResultsMCCLm{i}'
        averages, mse, other_averages = read_one_set(folder, file_name, pis, true_pi)
        estimator_data = estimator_data[:(i + 1)] + [averages] + estimator_data[(i + 1):]
        estimator_data += [mse]
        single_data += [other_averages]

    other_data = np.vstack(single_data)
    other_data = np.hstack([other_data, np.zeros((3 + max_m, 2 * (max_m + 1)))])

    out = np.hstack([np.atleast_2d(i).T for i in estimator_data])

    out = np.vstack([out, other_data])

    with open(f'simResults/{file_name}.txt', 'w') as f:
        f.write('')

    with open(f'simResults/{file_name}.txt', 'a') as f:
        new_line = '\n'
        tab = '\t'
        f.write(new_line.join([tab.join([str(j) for j in i]) for i in out]))
        f.write(new_line)


def format_all():
    print('formatting all data')
    test_l = [3, 4, 5]
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

    for seq_len in test_l:
        for j in range(4):
            for k in range(3):
                file_name = f'l{seq_len}q{j + 1}pi{k + 1}'
                print(file_name)
                format_file(file_name, seq_len, test_pi[seq_len][k])


if __name__ == '__main__':
    format_all()
