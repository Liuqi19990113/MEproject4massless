# %%
import numpy as np
import os
import itertools as it
from multiprocessing.pool import Pool

# %%
def generate_data(parameter_list):
    event_order = parameter_list[0]
    lambda_e = parameter_list[1]
    lambda_pi = parameter_list[2]
    gamma1 = parameter_list[3]
    gamma2 = parameter_list[4]
    os.system('./run.exe {} {} {} {} > event_{}.dat'.format(lambda_e, lambda_pi, gamma1, gamma2,event_order))

# %%
def generate_config():
    lambda_e_max = 1.2
    lambda_e_min = 0.05
    lambda_e_array = np.linspace(lambda_e_min, lambda_e_max, 50)
    final_config = []
    config_count = 0
    for lambda_e in lambda_e_array:
        lambda_pi_max = 0.25*lambda_e
        lambda_pi_min = -0.25*lambda_e
        lambda_pi_array = np.linspace(lambda_pi_min, lambda_pi_max, 20)
        gamma1_max = 0.25*lambda_e
        gamma1_min = -0.25*lambda_e
        gamma1_array = np.linspace(gamma1_min, gamma1_max, 20)
        gamma2_max = 0.25*lambda_e
        gamma2_min = -0.25*lambda_e
        gamma2_array = np.linspace(gamma2_min, gamma2_max, 20)
        tem_config = it.product([lambda_e], lambda_pi_array, gamma1_array, gamma2_array)
        for thing in tem_config:
            config_count += 1
            final_config.append([config_count]+list(thing))
    return final_config

# %%
if __name__ == '__main__':
    final_config = generate_config()
    with Pool() as pool:
        pool.map(generate_data, final_config)


