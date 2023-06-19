import math
import time

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import distance


def pi_seg(start, stop, match, pi, seq_len, memo):
    address = str(start) + str(stop) + str(match)
    if address in memo.keys():
        return memo[address]

    if start >= stop:
        return 1

    total = 0
    for i in range(len(pi)):
        if ind_to_seq(i, seq_len)[start:stop] == match:
            total += pi[i]

    memo[address] = total
    return total


def get_ind(seq):
    # helper function to get the index corresponding to sequence seq
    return int(seq, 2)


def ind_to_seq(ind, seq_len):
    return bin(ind)[2:].zfill(seq_len)


def f(seq, fk, qk, pi, possible_k, pi_memo, f_memo):
    if seq in f_memo.keys():
        return f_memo[seq]

    a = 0
    for k in possible_k:
        a += fk[get_ind(k)](seq, pi, pi_memo) * qk[get_ind(k)]

    f_memo[seq] = a
    return a


def t(seq, k, fk, qk, pi_old, possible_k, pi_old_memo, f_memo, t_memo):
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


def make_q_array(q, all_k):
    arr = []
    for k in all_k:
        arr.append(q ** k.count('0') * (1 - q) ** k.count('1'))

    return arr


def prod_creator(start, stop):
    return lambda seq, pi, memo: pi_seg(start, stop, seq[start:stop], pi, len(seq), memo)


def f_creator(prods):
    return lambda seq, pi, memo: math.prod([s(seq, pi, memo) for s in prods])


def make_f_array(all_k):
    arr = []

    for k in all_k:
        recombination_locs = [0]
        for l in range(len(k)):
            if k[l] == '0':
                recombination_locs.append(l + 1)
        
        recombination_locs.append(len(k) + 1)

        prods = []
        for l in range(len(recombination_locs) - 1):
            curr = recombination_locs[l]
            second = recombination_locs[l + 1]
            prods.append(prod_creator(curr, second))

        arr.append(f_creator(prods))

    return arr


def make_i_array(seq_len):
    arr = []
    for i in range(2 ** seq_len):
        arr.append(ind_to_seq(i, seq_len))

    return arr


def make_k_array(seq_len):
    arr = []
    for i in range(2 ** (seq_len - 1)):
        arr.append(ind_to_seq(i, seq_len - 1))

    return arr


def q_t(x, pi_old, counts, qk, fk, all_i, all_k, pi_old_memo, f_memo, t_memo):
    pi = np.append(x, 1 - sum(x))

    pi_memo = {}

    total = 0
    for i in all_i:
        if counts[get_ind(i)] == 0:
            continue
        for k in all_k:
            total += t(i, k, fk, qk, pi_old, all_k, pi_old_memo, f_memo, t_memo) * math.log(fk[get_ind(k)](i, pi, pi_memo)) \
                     * counts[get_ind(i)]

    # multiply by -1 because the optimization function minimizes
    return total * -1


def main(seq_len, q):
    t1 = time.perf_counter()
    # gets the observed data
    counts = np.array(get_counts('sequences' + str(seq_len) + '.txt', seq_len))
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

    print(guess)
    for run in range(40):
        t2 = time.perf_counter()
        # minimizes the function in (2 ** seq_len) - 1 dimensions in order to guarantee sum is 1
        d = minimize(q_t, guess[:-1], method="Nelder-Mead", bounds=[(0, 1) for _ in range(2 ** seq_len - 1)],
                     args=(guess, counts, qk, fk, all_i, all_k, {}, {}, {}))

        # if there is no significant change
        if distance.euclidean(np.append(d.x, 1 - sum(d.x)), guess) < 10 ** -7:
            break

        # the minimum is the new guess
        guess = np.append(d.x, 1 - sum(d.x))
        print(guess)
        print(time.perf_counter() - t2)

    print(guess)
    print(time.perf_counter() - t2)
    print(time.perf_counter() - t1)


if __name__ == '__main__':
    main(6, 0.05)
