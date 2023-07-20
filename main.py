import simulator
import dataInterpreter
import seqGenerator
import r_interface
import emAlgorithm


if __name__ == '__main__':
    # TO RUN SIMULATIONS
    test_l = [3, 4, 5]
    test_q = [0.001, 0.01, 0.05, 0.1]
    test_pi = {
        3: [
            [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125],
            [0.04, 0.06, 0.08, 0.12, 0.14, 0.16, 0.18, 0.22],
            [0, 0.02, 0.03, 0.05, 0.08, 0.12, 0.3, 0.4]
        ],
        4: [
            [0.0625 for _ in range(16)],
            [0.01, 0.02, 0.03, 0.04, 0.05, 0.03, 0.05, 0.07, 0.06, 0.08, 0.07, 0.09, 0.08, 0.1, 0.1, 0.12],
            [0, 0.02, 0.03, 0.05, 0, 0.06, 0.05, 0.09, 0, 0.08, 0.07, 0, 0.1, 0.1, 0.15, 0.2]
        ],
        5: [
            [0.03125 for _ in range(32)],
            [0.004, 0.006, 0.008, 0.012, 0.014, 0.016, 0.018, 0.022, 0.023, 0.024, 0.013, 0.017, 0.026, 0.027, 0.032,
             0.038, 0.028, 0.032, 0.036, 0.044, 0.033, 0.037, 0.042, 0.048, 0.038, 0.042, 0.046, 0.048, 0.052, 0.054,
             0.056, 0.064],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0.03, 0, 0.05, 0, 0.07, 0, 0.06, 0, 0.08, 0, 0.07, 0, 0.09, 0.01, 0.08,
             0.02, 0.1, 0.03, 0.1, 0.04, 0.12]
        ]
    }
    test_m = {
        3: [1],
        4: [1, 2],
        5: [1, 2, 3]
    }
    # seqGenerator.generate_all(test_l, test_q, test_pi)
    simulator.run_all_sims(test_l, test_q, test_pi)
    r_interface.run_all(test_l, test_q, test_pi, test_m)
    dataInterpreter.format_all(test_l, test_q, test_pi)

    # TO RUN ONE DATA SET
    # print(emAlgorithm.main(3, 0.05, 'file.txt'))
    pass

