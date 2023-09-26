import utilities


def convert_to_markov(true_pi, seq_len, m):

    sequences = utilities.make_i_array(seq_len)

    markov_likelihoods = []

    for seq in sequences:
        markov_likelihoods.append(get_one_markov(true_pi, seq, seq_len, m))

    return markov_likelihoods


def get_one_markov(true_pi, seq, seq_len, m):

    prob = 1
    for i in range(0, seq_len - m):
        prob *= utilities.pi_seg(i, i + m + 1, seq[i:i + m + 1], true_pi, seq_len)
        if i != 0:
            prob /= utilities.pi_seg(i, i + m, seq[i:i + m], true_pi, seq_len)

    return prob


if __name__ == '__main__':
    print(sum([0.2, 0.1, 0.2, 0.05, 0.1, 0.15, 0.12, 0.08]))
    print(convert_to_markov([0.2, 0.1, 0.2, 0.05, 0.1, 0.15, 0.12, 0.08], 3, 1))
