import math
import numpy as np
from scipy.optimize import minimize


def opp(i):
    # returns the opposite of i
    return 0 if i == 1 else 1


def get_ind(i, j):
    # helper function to get the index corresponding to sequence ij
    return int(str(i) + str(j), 2)


def f1(i, j, pi):
    # f1 is probability with recombination
    return (pi[get_ind(i, j)] + pi[get_ind(i, opp(j))]) * (pi[get_ind(i, j)] + pi[get_ind(opp(i), j)])


def f2(i, j, pi):
    # f2 is probability with no recombination
    return pi[get_ind(i, j)]


def f(i, j, qk, pi, fk):
    # f is the incomplete likelihood
    a = 0
    for k in range(2):
        a += fk[k](i, j, pi) * qk[k]

    return a


def t(i, j, k, fk, qk, pi_old):
    # the T function we calculated
    return fk[k](i, j, pi_old) * qk[k] / f(i, j, qk, pi_old, fk)


def get_counts(file):
    # counts the number of appearances of each sequence
    fi = open(file, 'r')

    counts = [0, 0, 0, 0]

    for line in fi:
        counts[int(line.strip(), 2)] += 1

    fi.close()

    return counts


def q_t(x, pi_old, counts, qk, fk):
    # adds a probability for 11
    pi = np.append(x, 1 - sum(x))

    # evaluates the sum that we calculated by hand
    total = 0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                total += t(i, j, k, fk, qk, pi_old) * math.log(fk[k](i, j, pi)) * counts[get_ind(i, j)]

    # multiply by -1 because the optimization function minimizes
    return total * -1


def main(file, q):
    # gets the observed data
    counts = np.array(get_counts(file))
    # our initial guess is just observed data
    guess = counts / sum(counts)
    # sets q1, q2
    qk = [q, 1 - q]
    # array to easily reference f1, f2
    fk = [f1, f2]
    print(guess)
    for run in range(40):
        # minimizes the function in 3 dimensions in order to guarantee sum is 1
        d = minimize(q_t, guess[:-1], method="Nelder-Mead", bounds=[(0, 1), (0, 1), (0, 1)], args=(guess, counts, qk, fk))
        # the minimum is the new guess
        guess = np.append(d.x, 1 - sum(d.x))
        print(guess)

    print(guess)


if __name__ == '__main__':
    main('oldData/sequences2.txt', 0.05)
