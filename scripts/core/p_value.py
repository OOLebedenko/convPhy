import numpy as np
import pandas as pd
from tqdm import tqdm


def get_p_value(R_S, number_of_choices=10000):
    resistance_branches = R_S["resistant"]
    p_values_out = []
    sampling = np.random.choice(resistance_branches, number_of_choices, replace=True)

    for resistant in tqdm(resistance_branches, desc="run p_value)"):
        p_value = (sampling >= resistant).sum()/ number_of_choices
        p_values_out.append(p_value)
    R_S['p_value'] = p_values_out
    return R_S[R_S['p_value'] < 0.05]
