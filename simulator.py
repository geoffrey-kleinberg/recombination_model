import emAlgorithm


def run_sim(sim_label, seq_len, q):
    with open(f'rawResultsEM/{sim_label}.txt', 'w') as f:
        f.write('')

    for i in range(1, 51):
        this_result = emAlgorithm.main(seq_len, q, f'simData/{sim_label}/s{i}.txt')
        tab = '\t'
        this_result[0] = tab.join([str(j) for j in this_result[0]])
        with open(f'rawResultsEM/{sim_label}.txt', 'a') as f:
            f.write(tab.join([str(j) for j in this_result]))
            f.write('\n')

        if i % 10 == 0:
            print(f'{i} out of 50')


def run_all_sims():
    test_l = [3, 4, 5]
    test_q = [0.001, 0.01, 0.05, 0.1]

    for seq_len in test_l:
        for j in range(len(test_q)):
            for k in range(1, 4):
                file_name = f'l{seq_len}q{j + 1}pi{k}'
                print(file_name)
                run_sim(file_name, seq_len, test_q[j])


if __name__ == '__main__':
    run_all_sims()
