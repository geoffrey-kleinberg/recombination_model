import random
import utilities
import os


# main function to generate data
def generate_one(seq_len, q, true_pi, n=100, runs=1, file_name_info=''):
    values = utilities.make_i_array(seq_len)
    all_k = utilities.make_k_array(seq_len)
    qk = utilities.make_q_array(q, all_k)
    fk = utilities.make_f_array(all_k)
    # set the weights of each descendant sequence according to recombination model
    weights = [utilities.f(i, fk, qk, true_pi, all_k) for i in values]

    if not os.path.isdir(f'simData/{file_name_info}'):
        os.makedirs(f'simData/{file_name_info}')

    for run in range(runs):
        # choose the correct number of data points from distribution
        data = random.choices(population=values, weights=weights, k=n)
        file_name = f'simData/{file_name_info}/s{run + 1}.txt'

        # writes all data to the file
        with open(file_name, 'w') as f:
            for d in data:
                f.write(d + "\n")


def generate_all(test_l, test_q, test_pi, n=100, runs=50):
    print('generating data')
    for seq_len in test_l:
        for q in range(len(test_q)):
            for pi in range(len(test_pi[seq_len])):
                generate_one(seq_len, test_q[q], test_pi[seq_len][pi], n=n, runs=runs,
                             file_name_info=f'l{seq_len}q{q + 1}pi{pi + 1}')


if __name__ == '__main__':
    generate_one(2, 0.05, [0.1, 0.2, 0.3, 0.4], n=100, runs=50, file_name_info='test')
