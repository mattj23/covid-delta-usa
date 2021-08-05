"""
    Tool to help estimate the probability distributions from the literature
"""
import math
from typing import List

import numpy as np
from scipy import stats
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

values = [
    (0.5069958847736888, 0.007294117647063891),
    (0.8296296296296539, 0.04011764705882778),
    (1.3827160493827364, 0.10894117647059068),
    (2.0510288065843785, 0.16470588235294217),
    (2.3967078189300555, 0.17600000000000066),
    (2.7193415637860205, 0.17882352941176533),
    (3.2032921810699673, 0.17352941176470663),
    (3.8485596707818974, 0.1555294117647071),
    (5.669135802469128, 0.09588235294117925),
    (9.033744855967049, 0.03376470588235736),
    (12.697942386831222, 0.011882352941181451),
    (17.652674897119258, 0.0034117647058875544),
    (19.86502057613159, 0.0023529411764758312),
]


def weib(a, c, loc, s, x: float) -> float:
    return stats.exponweib.pdf(x, a, c, loc, s)


def weib2(a, b, x) -> float:
    return (a/b) * (x/b)**(a-1) * np.exp(-(x/b)**a)


def opt2(x) -> np.array:
    residuals = []
    a, b, s = x
    for x_, y_ in values:
        v = weib2(a, b, x_/s) / s
        residuals.append(v - y_)
    return residuals


def opt(x, *args, **kwargs) -> np.array:
    residuals = []
    a, c, s = x
    for x_, y_ in values:
        residuals.append(weib(a, c, 0, s, x_) - y_)
    return residuals


def main():
    result = least_squares(opt, [1, 1, 1])
    print(result)
    a, c, s = result.x
    print(a, c, s)

    x_ = np.linspace(0, 21)
    # y = weib2(a, b, x_/s) / s
    y = stats.exponweib.pdf(x_, a, c, 0, s)
    plt.plot(*zip(*values), "o")
    plt.plot(x_, y)

    # z = stats.exponweib.pdf(x_, a, c, 0, s)
    # plt.plot(x_, z)

    plt.show()

    # Generate an inverse cdf table
    output = []
    for i in range(20):
        cdf = stats.exponweib.cdf(i+1, a, c, 0, s)
        output.append(cdf)
    print("{ " + ", ".join(f"{v:0.3f}" for v in output) + "}")



if __name__ == '__main__':
    main()
