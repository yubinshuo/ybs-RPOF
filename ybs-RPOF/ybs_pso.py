
import numpy as np
import pandas as pd
import os
from collections import Iterable, Counter

class PSO:

    def __init__(self, func, bound, POP_SIZE,x_start,w=1, c1=0.2, c2=0.2, v_max=0.05, *, var_name=None):
        bounds = Counter([isinstance(a, Iterable) for a in bound])[True]
        Var_size = int(np.ceil(POP_SIZE ** (1 / bounds)))
        vals = [np.linspace(var[0], var[1], Var_size) if isinstance(var, Iterable) else np.array([var]) for var in
                bound]
        for i in range(len(x_start)): vals[i][0] = x_start[i]
        for i in range(len(x_start)): vals[i][-1] = vals[i][-1] - 0.0005
        vals = np.array(list(map(lambda var: var.flatten(), np.meshgrid(*vals))))
        self.var_quantity, self.POP_SIZE = vals.shape
        self.func = func
        self.bound = bound
        self.w = w
        self.c1 = c1
        self.c2 = c2
        self.v_max = v_max
        self.var_name = var_name
        #
        self.particles = np.array(list(zip(*vals)))
        self.velocity = np.random.rand(*self.particles.shape) * v_max
        #
        self.person_best = self.particles.copy()
        #
        self.global_best = max(self.person_best, key=lambda particle: self.func(*particle)).copy()
        #
        self.global_best_fitness = 9999999999.9999

    def get_fitness(self):
        return np.array([self.func(*particle) for particle in self.particles])

    def update_position(self):
        for index, particle in enumerate(self.particles):
            V_k_plus_1 = self.w * self.velocity[index] \
                         + self.c1 * np.random.rand() * (self.person_best[index] - particle) \
                         + self.c2 * np.random.rand() * (self.global_best - particle)
            self.particles[index] = self.particles[index] + V_k_plus_1
            self.velocity[index] = V_k_plus_1
            for i, var in enumerate(particle):
                if isinstance(self.bound[i], Iterable):
                    if var < self.bound[i][0]:
                        self.particles[index][i] = self.bound[i][0]
                    elif var > self.bound[i][1]:
                        self.particles[index][i] = self.bound[i][1]
                elif var != self.bound[i]:
                    self.particles[index][i] = self.bound[i]

    def update_best(self):
        self.global_best_fitness = self.func(*self.global_best)
        person_best_value = np.array([self.func(*particle) for particle in self.person_best])
        for index, particle in enumerate(self.particles):
            current_particle_fitness = self.func(*particle)
            if current_particle_fitness > person_best_value[index]:
                person_best_value[index] = current_particle_fitness
                self.person_best[index] = particle
            if current_particle_fitness > self.global_best_fitness:
                self.global_best_fitness = current_particle_fitness
                self.global_best = particle
        with open(os.path.join(os.getcwd() , "pso_result"),'a') as ww:
            ww.write(str(self.global_best) + "  " + str(-self.global_best_fitness) + '\n')
    def pso(self):
        self.update_position()
        self.update_best()
        #print(self.global_best)
        return self.global_best, -self.global_best_fitness

    def info(self):
        result = pd.DataFrame(self.particles)
        if self.var_name == None:
            result.columns = [f'x{i}' for i in range(len(self.bound))]
        else:
            result.columns = self.var_name
        result['fitness'] = self.get_fitness()
        return result

