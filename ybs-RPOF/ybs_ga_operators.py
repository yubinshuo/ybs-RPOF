import numpy as np
#---------- crossover-------#
class crossover:
    def __init__(self,Chrom, size_pop, len_chrom):
        self.Chrom, self.size_pop, self.len_chrom = Chrom, size_pop, len_chrom
    def crossover_1point(self):
        Chrom, size_pop, len_chrom = self.Chrom, self.size_pop, self.len_chrom
        for i in range(0, size_pop, 2):
            n = np.random.randint(0, self.len_chrom)
            seg1, seg2 = self.Chrom[i, n:].copy(), self.Chrom[i + 1, n:].copy()
            self.Chrom[i, n:], self.Chrom[i + 1, n:] = seg2, seg1
        return self.Chrom
    def crossover_2point(self):
        Chrom, size_pop, len_chrom = self.Chrom, self.size_pop, self.len_chrom
        for i in range(0, size_pop, 2):
            n1, n2 = np.random.randint(0, self.len_chrom, 2)
            if n1 > n2:
                n1, n2 = n2, n1
            seg1, seg2 = self.Chrom[i, n1:n2].copy(), self.Chrom[i + 1, n1:n2].copy()
            self.Chrom[i, n1:n2], self.Chrom[i + 1, n1:n2] = seg2, seg1
        return self.Chrom
    def crossover_2point_bit(self):
        Chrom, size_pop, len_chrom = self.Chrom, self.size_pop, self.len_chrom
        half_size_pop = int(size_pop / 2)
        Chrom1, Chrom2 = Chrom[:half_size_pop], Chrom[half_size_pop:]
        mask = np.zeros(shape=(half_size_pop, len_chrom), dtype=int)
        for i in range(half_size_pop):
            n1, n2 = np.random.randint(0, self.len_chrom, 2)
            if n1 > n2:
                n1, n2 = n2, n1
            mask[i, n1:n2] = 1
        mask2 = (Chrom1 ^ Chrom2) & mask
        Chrom1 ^= mask2
        Chrom2 ^= mask2
        return self.Chrom
    def crossover_pmx(self):
        Chrom, size_pop, len_chrom = self.Chrom, self.size_pop, self.len_chrom
        for i in range(0, size_pop, 2):
            Chrom1, Chrom2 = self.Chrom[i], self.Chrom[i + 1]
            cxpoint1, cxpoint2 = np.random.randint(0, self.len_chrom - 1, 2)
            if cxpoint1 >= cxpoint2:
                cxpoint1, cxpoint2 = cxpoint2, cxpoint1 + 1
            pos1_recorder = {value: idx for idx, value in enumerate(Chrom1)}
            pos2_recorder = {value: idx for idx, value in enumerate(Chrom2)}
            for j in range(cxpoint1, cxpoint2):
                value1, value2 = Chrom1[j], Chrom2[j]
                pos1, pos2 = pos1_recorder[value2], pos2_recorder[value1]
                Chrom1[j], Chrom1[pos1] = Chrom1[pos1], Chrom1[j]
                Chrom2[j], Chrom2[pos2] = Chrom2[pos2], Chrom2[j]
                pos1_recorder[value1], pos1_recorder[value2] = pos1, j
                pos2_recorder[value1], pos2_recorder[value2] = j, pos2
            self.Chrom[i], self.Chrom[i + 1] = Chrom1, Chrom2
        return self.Chrom
#----------- mutation ---------------------------#
class  mutation:
    def __init__(self, Chrom, size_pop, len_chrom,prob_mut):
        self.Chrom, self.size_pop, self.len_chrom, self.prob_mut = Chrom, size_pop, len_chrom, prob_mut
    def mutation(self):
        mask = (np.random.rand(self.size_pop, self.len_chrom) < self.prob_mut)
        self.Chrom ^= mask
        return self.Chrom
    def mutation_TSP_1(self):
        for i in range(self.size_pop):
            for j in range(self.n_dim):
                if np.random.rand() < self.prob_mut:
                    n = np.random.randint(0, self.len_chrom, 1)
                    self.Chrom[i, j], self.Chrom[i, n] = self.Chrom[i, n], self.Chrom[i, j]
        return self.Chrom
    def swap(individual):
        n1, n2 = np.random.randint(0, individual.shape[0] - 1, 2)
        if n1 >= n2:
            n1, n2 = n2, n1 + 1
        individual[n1], individual[n2] = individual[n2], individual[n1]
        return individual
    def reverse(individual):
        n1, n2 = np.random.randint(0, individual.shape[0] - 1, 2)
        if n1 >= n2:
            n1, n2 = n2, n1 + 1
        individual[n1:n2] = individual[n1:n2][::-1]
        return individual
    def transpose(individual):
        n1, n2, n3 = sorted(np.random.randint(0, individual.shape[0] - 2, 3))
        n2 += 1
        n3 += 2
        slice1, slice2, slice3, slice4 = individual[0:n1], individual[n1:n2], individual[n2:n3 + 1], individual[n3 + 1:]
        individual = np.concatenate([slice1, slice3, slice2, slice4])
        return individual
    def mutation_reverse(self):
        for i in range(self.size_pop):
            if np.random.rand() < self.prob_mut:
                self.Chrom[i] = reverse(self.Chrom[i])
        return self.Chrom
    def mutation_swap(self):
        for i in range(self.size_pop):
            if np.random.rand() < self.prob_mut:
                self.Chrom[i] = swap(self.Chrom[i])
        return self.Chrom
# ------------------- ranking -------------------- #
class  ranking:
    def __init__(self, Y):
        self.Y = Y
    def ranking(self):
        self.FitV = -self.Y
    def ranking_linear(self):
        self.FitV = np.argsort(np.argsort(-self.Y))
        return self.FitV
# ------------------- selection --------------------#
class  selection:
    def __init__(self, Chrom, size_pop, len_chrom,FitV):
        self.Chrom, self.size_pop, self.len_chrom,self.FitV  = Chrom, size_pop, len_chrom,FitV
    def selection_tournament(self, tourn_size=3):
        FitV = self.FitV
        sel_index = []
        for i in range(self.size_pop):
            aspirants_index = np.random.randint(self.size_pop, size=tourn_size)
            sel_index.append(max(aspirants_index, key=lambda i: FitV[i]))
        self.Chrom = self.Chrom[sel_index, :]
        return self.Chrom
    def selection_tournament_faster(self, tourn_size=3):
        aspirants_idx = np.random.randint(self.size_pop, size=(self.size_pop, tourn_size))
        aspirants_values = self.FitV[aspirants_idx]
        winner = aspirants_values.argmax(axis=1)
        sel_index = [aspirants_idx[i, j] for i, j in enumerate(winner)]
        self.Chrom = self.Chrom[sel_index, :]

        return self.Chrom
    def selection_roulette_1(self):
        FitV = self.FitV
        FitV = FitV - FitV.min() + 1e-10
        sel_prob = FitV / FitV.sum()
        sel_index = np.random.choice(range(self.size_pop), size=self.size_pop, p=sel_prob)
        self.Chrom = self.Chrom[sel_index, :]
        return self.Chrom
    def selection_roulette_2(self):
        FitV = self.FitV
        FitV = (FitV - FitV.min()) / (FitV.max() - FitV.min() + 1e-10) + 0.2
        sel_prob = FitV / FitV.sum()
        sel_index = np.random.choice(range(self.size_pop), size=self.size_pop, p=sel_prob)
        self.Chrom = self.Chrom[sel_index, :]
        return self.Chrom
if __name__ == '__main__':
    pass