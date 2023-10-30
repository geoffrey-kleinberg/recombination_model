import math
import time

import utilities

import numpy as np
from scipy.optimize import minimize


# overall function for Q(pi | pi^t)
def q_t(x, pi_old, counts, qk, fk, all_i, all_k, removed, pi_old_memo=None, f_memo=None, t_memo=None):
    if t_memo is None:
        t_memo = {}
    if f_memo is None:
        f_memo = {}
    if pi_old_memo is None:
        pi_old_memo = {}

    pi = np.insert(x, removed, 1 - sum(x))

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


def main(seq_len, q, file_name='', view_print=False):
    if file_name == '':
        file_name = 'oldData/sequences' + str(seq_len) + '.txt'

    t = time.perf_counter()

    # gets the observed data
    counts = np.array(utilities.get_counts(file_name, seq_len))
    # our initial guess is just observed data
    guess = counts / sum(counts)

    # THERE IS AN ERROR WHEN THE INDEX OF GUESS NOT PASSED INTO Q HAS VALUE 0
    # THIS GUARANTEES THAT WON'T HAPPEN
    non_zero_location = np.where(np.isclose(guess, max(guess)))[0][0]

    # gets all possible DNA sequences
    all_i = utilities.make_i_array(seq_len)

    # gets all possible recombination locations
    all_k = utilities.make_k_array(seq_len)

    # sets qk values
    qk = utilities.make_q_array(q, all_k)

    # sets all the fk functions
    fk = utilities.make_f_array(all_k)

    # gets the likelihood
    prev_ll = utilities.total_likelihood(counts, guess, all_i, all_k, fk, qk)
    if view_print:
        print(prev_ll)
    # records the number of runs and optimization iterations
    runs = 0
    opt_it = 0
    opt_time = 0
    while True:
        runs += 1
        # minimizes the function in (2 ** seq_len) - 1 dimensions in order to guarantee sum is 1
        guess_arg = np.concatenate((guess[:non_zero_location], guess[(non_zero_location + 1):]))
        t1 = time.perf_counter()

        d = minimize(q_t, guess_arg, method="Nelder-Mead", bounds=[(0, 1) for _ in range(2 ** seq_len - 1)],
                     args=(guess, counts, qk, fk, all_i, all_k, non_zero_location, {}, {}, {}))
        opt_time += time.perf_counter() - t1

        opt_it += d.nit

        # the minimum is the new guess
        guess = np.insert(d.x, non_zero_location, 1 - sum(d.x))

        # calculate the new likelihood
        new_ll = utilities.total_likelihood(counts, guess, all_i, all_k, fk, qk)
        # if the improvement is small enough, end the algorithm
        improvement = (new_ll - prev_ll) / prev_ll
        if view_print:
            print("improvement ", improvement)
        if improvement <= 10 ** -6:
            break

        prev_ll = new_ll

        if view_print:
            print("time:", time.perf_counter() - t)
            print("guess:", guess)

    total_time = time.perf_counter() - t

    if view_print:
        print("time:", opt_time)
        print("guess:")
        for i in guess:
            print(i)

    return [guess, runs, opt_it, total_time]


if __name__ == '__main__':
    # main(3, 0.01, file_name='oldData/sequences3.txt', view_print=True)
    main(10, 0.0003, file_name='real_data.txt', view_print=True)
