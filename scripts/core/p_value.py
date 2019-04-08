import numpy as np
from tqdm import tqdm


def get_p_value(R_S, info_pos, n=10000):
    p_value = []
    with open(info_pos) as file:
        info_pos = [line.strip() for line in file.readlines()]
    for ind, real in tqdm(enumerate(R_S), desc="run p_value)"):
        counter = 0
        for i in range(n):
            point = np.random.randint(len(R_S))
            sampling = R_S[point]
            if int(sampling[1]) >= int(real[1]):
                counter += 1
        p = float(counter) / float(n)
        if p < 0.05:
            index = real[0]
            p_value.append(str(index) + "\t" + str(info_pos[int(index)]) + "\t" + str(p) + "\n")
    return (p_value)

