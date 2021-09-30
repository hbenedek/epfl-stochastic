from statistics import NormalDist
import numpy as np
import matplotlib.pyplot as plt


def pos_side(x: float) -> float:
    if x > 0: return x
    return 0

def simulate(time: int, account: float, K: float, number_of_puts: int, mu: float, sigma2: float):
    history = [0]
    prices = NormalDist(mu, sigma2).samples(time)
    cdf = NormalDist(mu, sigma2).cdf
    pdf = NormalDist(mu, sigma2).pdf
    P = (K - mu) * cdf(K) + np.sqrt(sigma2) * pdf(K)
    date = 0
    while account >= 0 and date < time:
        account += number_of_puts * (P - pos_side(K-prices[date]))
        history.append(account)
        date = date + 1
    return history


simulations= []
for i in range (100):
    history = simulate(1000, 0, 10, 1, 0, 1)
    simulations.append(history)

color = []
g = 0
r = 0
for s in simulations:
    if s[-1]<0: 
    #    color='red' 
       r = r + 1
    else: 
    #    color='green'
        g = g + 1
    plt.plot(np.arange(len(s)), s)

print(f'green{g}, red: {r}')
plt.show()

