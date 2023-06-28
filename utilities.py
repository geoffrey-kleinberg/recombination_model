import math
from itertools import combinations


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
def f(seq, fk, qk, pi, possible_k, pi_memo=None, f_memo=None):
    if f_memo is None:
        f_memo = {}
    if pi_memo is None:
        pi_memo = {}
    if seq in f_memo.keys():
        return f_memo[seq]

    # sums f_k * q_k along all k
    a = 0
    for k in possible_k:
        a += fk[get_ind(k)](seq, pi, pi_memo) * qk[get_ind(k)]

    f_memo[seq] = a
    return a


# the t function that calculates the expectation of r
def t(seq, k, fk, qk, pi_old, possible_k, pi_old_memo=None, f_memo=None, t_memo=None):
    # for memoization
    if t_memo is None:
        t_memo = {}
    if f_memo is None:
        f_memo = {}
    if pi_old_memo is None:
        pi_old_memo = {}
    address = str(seq) + str(k)
    if address in t_memo.keys():
        return t_memo[address]

    # the T function we calculated
    t_memo[address] = fk[get_ind(k)](seq, pi_old, pi_old_memo) * \
                      qk[get_ind(k)] / f(seq, fk, qk, pi_old, possible_k, pi_old_memo, f_memo)

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
        # 1 is no recombination, occurs with probability 1 - q
        arr.append(q ** k.count('0') * (1 - q) ** k.count('1'))

    return arr


def prod_creator(start, stop):
    # returns a function that calculates one of the factors in the f_k
    return lambda seq, pi, memo: pi_seg(start, stop, seq[start:stop], pi, len(seq), memo)


def f_creator(prods):
    # returns a function that calculates a given f_k
    return lambda seq, pi, memo: math.prod([s(seq, pi, memo) for s in prods])


def f_sum(to_sum):
    return lambda seq, pi, memo: sum([f_a(seq, pi, memo) for f_a in to_sum])


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

        print(recombination_locs)
        # generates the factors in f_k based on the recombination locations
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


def make_k_array_2(seq_len):
    arr = []
    bin_len = math.ceil(math.log2(seq_len))
    for i in range(seq_len):
        arr.append(ind_to_seq(i, bin_len))

    return arr


def make_q_array_2(q, all_k):
    arr = []
    highest = get_ind(all_k[-1])
    for i in all_k:
        num = get_ind(i)
        arr.append(q ** (highest - num) * (1 - q) ** num)

    return arr


def make_f_array_2(all_k):
    arr = []
    num_locations = get_ind(all_k[-1])

    for i in all_k:
        num_recombinations = num_locations - get_ind(i)
        possible_recombination_locations = list(combinations(range(1, num_locations + 1), num_recombinations))
        to_sum = []
        for r in possible_recombination_locations:
            locations = list(r)
            locations.insert(0, 0)
            locations.append(num_locations + 1)

            prods = []
            for l in range(len(locations) - 1):
                curr = locations[l]
                second = locations[l + 1]
                prods.append(prod_creator(curr, second))

            to_sum.append(f_creator(prods))

        arr.append(f_sum(to_sum))

    return arr


if __name__ == '__main__':
    pass
