import emAlgorithm


# runs a simulation (EM algorithm) for all data sets of one l, q, and pi
def run_sim(sim_label, seq_len, q, num=50):
    # clear the results
    with open(f'rawResultsEM/{sim_label}.txt', 'w') as f:
        f.write('')

    # for each data set
    for i in range(1, num + 1):
        this_result = emAlgorithm.main(seq_len, q, f'simData/{sim_label}/s{i}.txt')
        tab = '\t'
        # write all the data to the file
        this_result[0] = tab.join([str(j) for j in this_result[0]])
        with open(f'rawResultsEM/{sim_label}.txt', 'a') as f:
            f.write(tab.join([str(j) for j in this_result]))
            f.write('\n')

        # counter to see progress
        if i % 10 == 0:
            print(f'{i} out of {num}')


def run_all_sims(test_l, test_q, test_pi):
    # runs all simulations over a specified range of l and q
    for seq_len in test_l:
        for j in range(len(test_q)):
            for k in range(len(test_pi[seq_len])):
                file_name = f'l{seq_len}q{j + 1}pi{k + 1}'
                print(file_name)
                run_sim(file_name, seq_len, test_q[j])
