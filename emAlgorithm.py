import math
import utilities

import numpy as np
from scipy.optimize import minimize


# overall function for Q(pi | pi^t)
def q_t(x, pi_old, counts, qk, fk, all_i, all_k, pi_old_memo=None, f_memo=None, t_memo=None):
    if t_memo is None:
        t_memo = {}
    if f_memo is None:
        f_memo = {}
    if pi_old_memo is None:
        pi_old_memo = {}
    pi = np.append(x, 1 - sum(x))

    pi_memo = {}

    total = 0
    # for all possible sequences
    for i in all_i:
        # if there are none we can skip it
        if counts[utilities.get_ind(i)] == 0:
            continue
        # for all possible recombinations
        for k in all_k:
            # add to the total likelihood according to equation
            total += utilities.t(i, k, fk, qk, pi_old, all_k, pi_old_memo, f_memo, t_memo) * \
                     math.log(fk[utilities.get_ind(k)](i, pi, pi_memo)) \
                     * counts[utilities.get_ind(i)]

    # multiply by -1 because the optimization function minimizes
    return total * -1


def main(seq_len, q, file_name=''):
    if file_name == '':
        file_name = 'oldData/sequences' + str(seq_len) + '.txt'
    # t1 = time.perf_counter()
    # gets the observed data
    counts = np.array(utilities.get_counts(file_name, seq_len))
    # our initial guess is just observed data
    guess = counts / sum(counts)

    # gets all possible DNA sequences
    all_i = utilities.make_i_array(seq_len)
    # gets all possible recombination locations
    all_k = utilities.make_k_array(seq_len)

    # sets qk values
    qk = utilities.make_q_array(q, all_k)
    # sets all the fk functions
    fk = utilities.make_f_array(all_k)

    # print(guess)
    prev_ll = utilities.total_log_likelihood(counts, guess, all_i, all_k, fk, qk)
    for run in range(40):
        # t2 = time.perf_counter()
        # minimizes the function in (2 ** seq_len) - 1 dimensions in order to guarantee sum is 1
        d = minimize(q_t, guess[:-1], method="Nelder-Mead", bounds=[(0, 1) for _ in range(2 ** seq_len - 1)],
                     args=(guess, counts, qk, fk, all_i, all_k, {}, {}, {}))

        # the minimum is the new guess
        guess = np.append(d.x, 1 - sum(d.x))

        # print(guess)
        new_ll = utilities.total_log_likelihood(counts, guess, all_i, all_k, fk, qk)
        improvement = new_ll - prev_ll
        # print(improvement)
        if improvement <= 10 ** -6:
            break

        prev_ll = new_ll

    return guess
    # print(time.perf_counter() - t2)

    # print(time.perf_counter() - t2)
    # print(time.perf_counter() - t1)


if __name__ == '__main__':
    print(main(2, 0.05, file_name='oldData/sequences2.txt'))
