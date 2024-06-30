import simulator
import dataInterpreter
import seqGenerator
import r_interface
import markovChain
import emAlgorithm
import utilities


def run_full_simulation(sim_l, sim_q, sim_pi, sim_m, generate_new_data=False):
    if generate_new_data:
        seqGenerator.generate_all(sim_l, sim_q, sim_pi)

    simulator.run_all_sims(sim_l, sim_q, sim_pi)
    r_interface.run_all(sim_l, sim_q, sim_pi, sim_m)
    dataInterpreter.format_all(sim_l, sim_q, sim_pi)


if __name__ == '__main__':
    # TO RUN SIMULATIONS
    test_l = [3, 4, 5]
    test_q = [0.001, 0.01, 0.05, 0.1]
    three = [0.02, 0.08, 0.18, 0.22, 0.22, 0.18, 0.08, 0.02]
    four = [0.005, 0.015, 0.03, 0.05, 0.08, 0.1, 0.105, 0.115, 0.115, 0.105, 0.1, 0.08, 0.05, 0.03, 0.015, 0.005]
    five = [0.001, 0.002, 0.004, 0.006, 0.01, 0.014, 0.018, 0.02, 0.025, 0.035, 0.04, 0.055, 0.06, 0.065, 0.07, 0.075,
            0.075, 0.07, 0.065, 0.06, 0.055, 0.04, 0.035, 0.025, 0.02, 0.018, 0.014, 0.01, 0.006, 0.004, 0.002, 0.001]
    test_pi = {
        3: [
            three,
            markovChain.convert_to_markov(three, 3, 1)
        ],
        4: [
            four,
            markovChain.convert_to_markov(four, 4, 1),
            markovChain.convert_to_markov(four, 4, 2)
        ],
        5: [
            five,
            markovChain.convert_to_markov(five, 5, 1),
            markovChain.convert_to_markov(five, 5, 2),
            markovChain.convert_to_markov(five, 5, 3)
        ]
    }
    test_m = {
        3: [1],
        4: [1, 2],
        5: [1, 2, 3]
    }

    run_full_simulation(test_l, test_q, test_pi, test_m, generate_new_data=False)

    # TO RUN ONE DATA SET
    # print(emAlgorithm.main(3, 0.05, 'file.txt'))

