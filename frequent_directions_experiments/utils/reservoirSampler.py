from __future__ import absolute_import
from __future__ import print_function
from numpy.random import binomial
from random import sample


class ReservoirSampler:
    def __init__(self, t_paralel_sampleres=1):
        self.t = t_paralel_sampleres
        self.t_range = list(range(self.t))
        self.items = [None] * self.t
        self.items_weights = [0.0] * self.t

        self.item_probability = 0.0
        self.sum_w = 0.0
        self.machine_precision = 1e-10

    def add(self, item, w=1):

        w = float(w)
        if w <= 0.0:
            return
        self.sum_w += w
        p = w / max(self.sum_w, self.machine_precision)

        num_items_to_update = binomial(self.t, p)
        items_to_update = sample(self.t_range, num_items_to_update)

        for i in items_to_update:
            self.items[i] = item
            self.items_weights[i] = w

    def get(self, with_probabilities=False):
        if with_probabilities:
            probs = [w / self.sum_w for w in self.items_weights]
            return list(zip(self.items, probs))
        else:
            return self.items


if __name__ == "__main__":
    n = 1000
    items = list(range(n))
    weights = list(range(n))

    rs = ReservoirSampler(1000)

    for i in range(n):
        rs.add(items[i], weights[i])

    print(sorted(rs.get()))
