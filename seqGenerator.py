import random
import emAlgorithm


def main(l, q, true_pi, n, runs=1, file_name_info=''):
    values = emAlgorithm.make_i_array(l)
    all_k = emAlgorithm.make_k_array(l)
    qk = emAlgorithm.make_q_array(q, all_k)
    fk = emAlgorithm.make_f_array(all_k)
    weights = [emAlgorithm.f(i, fk, qk, true_pi, all_k) for i in values]

    for run in range(runs):

        data = random.choices(population=values, weights=weights, k=n)
        file_name = f'simData/{file_name_info}s{run + 1}.txt'

        with open(file_name, 'w') as f:
            for d in data:
                f.write(d + "\n")


if __name__ == '__main__':
    main(2, 0.05, [0.1, 0.2, 0.3, 0.4], 100, runs=2, file_name_info='pi1q005l2')
