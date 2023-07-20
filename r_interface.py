import subprocess
import os


def run_all(test_l, test_q, test_pi, test_m):
    print('calculating hierarchical estimator')
    for seq_len in test_l:
        for j in range(len(test_q)):
            for k in range(len(test_pi[seq_len])):
                file_name = f'l{seq_len}q{j + 1}pi{k + 1}'
                print(file_name)
                for m in test_m[seq_len]:
                    if not os.path.isdir(f'rawResultsMCCLm{m}'):
                        os.mkdir(f'rawResultsMCCLm{m}')
                    run_one(seq_len, test_q[j], m, file_name)


def run_one(seq_len, q, m, file_name, n=100):
    command = 'Rscript'
    file = 'sim_interface.R'
    args = [seq_len, q, m, n, file_name]

    cmd = [command, file] + [str(i) for i in args]

    subprocess.run(cmd)


if __name__ == '__main__':
    pass
