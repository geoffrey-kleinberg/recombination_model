import math
import time

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import distance


def total_log_likelihood(counts, pi, all_i, all_k, fk, qk):
    likelihood = 0
    pi_memo = {}
    for seq in all_i:
        likelihood += counts[get_ind(seq)] * math.log(f(seq, fk, qk, pi, all_k, pi_memo, {}))

    return likelihood


# gets the proportion of sequences that have string match from start to stop
def pi_seg(start, stop, match, pi, seq_len, memo={}):
    # for memoization
    address = str(start) + str(stop) + str(match)
    if address in memo.keys():
        return memo[address]

    if start >= stop:
        return 1

    # sums along pi if the sequences match in the desired location
    total = 0
    for i in range(len(pi)):
        if ind_to_seq(i, seq_len)[start:stop] == match:
            total += pi[i]

    # memoization and return
    memo[address] = total
    return total


def get_ind(seq):
    # helper function to get the index corresponding to sequence seq
    return int(seq, 2)


def ind_to_seq(ind, seq_len):
    # takes an index and returns the corresponding sequence
    return bin(ind)[2:].zfill(seq_len)


# the f function, the complete likelihood for a given sequence
def f(seq, fk, qk, pi, possible_k, pi_memo={}, f_memo={}):
    if seq in f_memo.keys():
        return f_memo[seq]

    # sums f_k * q_k along all k
    a = 0
    for k in possible_k:
        a += fk[get_ind(k)](seq, pi, pi_memo) * qk[get_ind(k)]

    f_memo[seq] = a
    return a


# the t function that calculates the expectation of r
def t(seq, k, fk, qk, pi_old, possible_k, pi_old_memo={}, f_memo={}, t_memo={}):
    # for memoization
    address = str(seq) + str(k)
    if address in t_memo.keys():
        return t_memo[address]

    # the T function we calculated
    t_memo[address] = fk[get_ind(k)](seq, pi_old, pi_old_memo) \
                      * qk[get_ind(k)] / f(seq, fk, qk, pi_old, possible_k, pi_old_memo, f_memo)

    return t_memo[address]


def get_counts(file, seq_len):
    # counts the number of appearances of each sequence
    fi = open(file, 'r')

    counts = [0 for _ in range(2 ** seq_len)]

    for line in fi:
        counts[get_ind(line.strip())] += 1

    fi.close()

    return counts


# makes an array of all q_k based on q and number of recombinations
def make_q_array(q, all_k):
    arr = []
    for k in all_k:
        # 0 is recombination and occurs with probability q
        # 1 is no recombination, occures with probability 1 - q
        arr.append(q ** k.count('0') * (1 - q) ** k.count('1'))

    return arr


def prod_creator(start, stop):
    # returns a function that calculates one of the factors in the f_k
    return lambda seq, pi, memo: pi_seg(start, stop, seq[start:stop], pi, len(seq), memo)


def f_creator(prods):
    # returns a function that calculates a given f_k
    return lambda seq, pi, memo: math.prod([s(seq, pi, memo) for s in prods])


# makes all the f_k functions
def make_f_array(all_k):
    arr = []

    # for every possible permutation of recombination we need an f_k
    for k in all_k:
        # stores the locations where recombination occurs
        recombination_locs = [0]
        for l in range(len(k)):
            if k[l] == '0':
                recombination_locs.append(l + 1)

        # array has 0 and len(k) + 1 at start and end
        recombination_locs.append(len(k) + 1)

        # generates the factors in f_k based ont the recombination locations
        prods = []
        for l in range(len(recombination_locs) - 1):
            curr = recombination_locs[l]
            second = recombination_locs[l + 1]
            prods.append(prod_creator(curr, second))

        # generates the final f_k function
        arr.append(f_creator(prods))

    return arr


# generates all possible sequences of length seq_len
def make_i_array(seq_len):
    arr = []
    for i in range(2 ** seq_len):
        arr.append(ind_to_seq(i, seq_len))

    return arr


# generates all possible recombination status for DNA of length seq_len
def make_k_array(seq_len):
    arr = []
    for i in range(2 ** (seq_len - 1)):
        arr.append(ind_to_seq(i, seq_len - 1))

    return arr


# overall function for Q(pi | pi^t)
def q_t(x, pi_old, counts, qk, fk, all_i, all_k, pi_old_memo={}, f_memo={}, t_memo={}):
    pi = np.append(x, 1 - sum(x))

    pi_memo = {}

    total = 0
    # for all possible sequences
    for i in all_i:
        # if there are none we can skip it
        if counts[get_ind(i)] == 0:
            continue
        # for all possible recombinations
        for k in all_k:
            # add to the total likelihood according to equation
            total += t(i, k, fk, qk, pi_old, all_k, pi_old_memo, f_memo, t_memo) * math.log(fk[get_ind(k)](i, pi, pi_memo)) \
                     * counts[get_ind(i)]

    # multiply by -1 because the optimization function minimizes
    return total * -1


def main(seq_len, q, file_name=''):
    if file_name == '':
        file_name = 'oldData/sequences' + str(seq_len) +'.txt'
    # t1 = time.perf_counter()
    # gets the observed data
    counts = np.array(get_counts(file_name, seq_len))
    # our initial guess is just observed data
    guess = counts / sum(counts)

    # gets all possible DNA sequences
    all_i = make_i_array(seq_len)
    # gets all possible recombination locations
    all_k = make_k_array(seq_len)

    # sets qk values
    qk = make_q_array(q, all_k)
    # sets all the fk functions
    fk = make_f_array(all_k)

    #print(guess)
    prev_ll = total_log_likelihood(counts, guess, all_i, all_k, fk, qk)
    for run in range(40):
        # t2 = time.perf_counter()
        # minimizes the function in (2 ** seq_len) - 1 dimensions in order to guarantee sum is 1
        d = minimize(q_t, guess[:-1], method="Nelder-Mead", bounds=[(0, 1) for _ in range(2 ** seq_len - 1)],
                     args=(guess, counts, qk, fk, all_i, all_k, {}, {}, {}))


        # the minimum is the new guess
        guess = np.append(d.x, 1 - sum(d.x))

        #print(guess)
        new_ll = total_log_likelihood(counts, guess, all_i, all_k, fk, qk)
        improvement = new_ll - prev_ll
        #print(improvement)
        if improvement <= 10 ** -6:
            break

        prev_ll = new_ll

    return guess
        # print(time.perf_counter() - t2)

    # print(time.perf_counter() - t2)
    # print(time.perf_counter() - t1)


if __name__ == '__main__':
    main(2, 0.05, file_name='oldData/sequences2.txt')
