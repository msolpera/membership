
import numpy as np

# def probCnvrg(probs_all, probs_old, perc_prob_cnvrg, runs_old, min_runs=3):
#     """
#     Check if probabilities converged to within 'perc_prob_cnvrg'%.

#     Break only if a minimum of 'min_runs' consecutive runs have been processed.
#     """
#     cnvrg_flag = False

#     if not probs_all:
#         return cnvrg_flag, probs_old, runs_old

#     # Average all probabilities.
#     prob_avrg = np.mean(probs_all, 0)

#     # if runs_old > min_runs:
#     #     for _ in np.linspace(0.01, .5, 50):
#     #         if np.allclose(probs_old, prob_avrg, _):
#     #             print("Relative tolerance for probabilities: {:.0f}%".format(
#     #                 100. * _))
#     #             break

#     # Check if probabilities changed less than 'perc_prob_cnvrg'% with
#     # respect to the previous iteration.
#     if np.allclose(probs_old, prob_avrg, perc_prob_cnvrg):
#         if runs_old > min_runs:
#             # Arrays are equal within tolerance
#             cnvrg_flag = True
#         runs_old += 1
#     else:
#         # Store new array in old one and proceed to new iteration.
#         probs_old = prob_avrg
#         # Reset
#         runs_old = 0

#     return cnvrg_flag, probs_old, runs_old
