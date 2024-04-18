"""Provide the primary functions."""

import numpy as np
import math

class BitString:
    """
    Simple class to implement a config of bits
    """
    def __init__(self, N):
        self.N = N
        self.config = np.zeros(N, dtype=int) 
        


    def __repr__(self):
        ans = ""
        for i in self.config:
            ans += str(i)
        return ans;

    def __eq__(self, other):        
        return self.int() == other.int()
    
    def __len__(self):
        return self.N

    def on(self):
        count = 0
        for i in self.config:
            if i == 1:
                count+=1
        return count

    
    def off(self):
        count = 0
        for i in self.config:
            if i == 0:
                count+=1
        return count
    
    def flip_site(self,i):
        self.config[i] = (self.config[i] - 1) * -1
    
    def int(self):
        ans = 0
        for i in range(len(self.config)):
            ans += (self.config[len(self.config) - i - 1] * (pow(2, i)))
        return ans

 

    def set_config(self, s):
        self.config = np.array(s)
        
    def set_int_config(self, dec:int):
        self.config = np.zeros(self.N, dtype=int)
        num = dec
        index = len(self.config) - 1
        while num != 0:
            self.config[index] = num % 2
            num = num // 2
            index -= 1

class IsingHamiltonian:
    def __init__(self, J=[[()]], mu=np.zeros(1)):
        """Constructor

        Parameters
        ----------
        J: list of lists of tuples, optional
            Strength of coupling, e.g,
            [(4, -1.1), (6, -.1)]
            [(5, -1.1), (7, -.1)]
        mu: vector, optional
            local fields
        """
        self.J = J
        self.mu = mu
        self.N = len(self.J)
        # self.nodes = []
        # self.js = []

        # for i in range(len(self.J)):
        #     self.nodes.append(np.zeros(len(self.J[i]), dtype=int))
        #     self.js.append(np.zeros(len(self.J[i])))
        #     for jidx, j in enumerate(self.J[i]):
        #         self.nodes[i][jidx] = j[0]
        #         self.js[i][jidx] = j[1]
        # self.mu = np.array([i for i in self.mu])
        # self.N = len(self.J)

    def energy(self, config):
        """Compute energy of configuration, `config`

            .. math::
                E = \\left<\\hat{H}\\right>

        Parameters
        ----------
        config   : BitString
            input configuration

        Returns
        -------
        energy  : float
            Energy of the input configuration
        """
        if len(config.config) != len(self.J):
            raise ValueError("wrong dimension")

        e = 0.0
        for i in range(config.N):
            # print()
            # print(i)
            for j in self.J[i]:
                if j[0] < i:
                    continue
                # print(j)
                if config.config[i] == config.config[j[0]]:
                    e += j[1]
                else:
                    e -= j[1]

        e += np.dot(self.mu, 2 * config.config - 1)
        return e
    def compute_average_values(self, T):


        bs = BitString(self.N)
                #Boltzman constant
        k =1# 1.38064852 * math.pow(10, -23)
        beta = 1/(k * T)

        #Calculate normalization constant Z
        Z = 0
        for config in range( 0,  2 ** bs.N ):
            bs.set_int_config(config)
            Z += math.exp(-1 * beta * self.energy(bs))

        E = 0
        M = 0
        E_square = 0
        M_square = 0
        
        for config in range( 0,  2 ** bs.N ):
            bs.set_int_config(config)
            e = self.energy(bs)
            P = math.exp(-1 * beta * e) / Z
            E += e * P
            E_square += e * e * P
            m = bs.config.tolist().count(1) - bs.config.tolist().count(0)
            M += m * P
            M_square += m * m * P

        HC = (E_square - (E * E)) / (T * T)
        MS = (M_square - (M * M)) / T

        
        return E, M, HC, MS


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print()
