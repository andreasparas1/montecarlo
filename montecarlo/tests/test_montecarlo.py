"""
Unit and regression test for the montecarlo package.
"""

# Import package, test suite, and other packages as needed
import montecarlo
import pytest
import sys
import random
import numpy as np
import copy as cp

def test_montecarlo_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "montecarlo" in sys.modules

def test_BS_1():
    my_bs = montecarlo.BitString(8)
    my_bs.flip_site(2)
    my_bs.flip_site(2)
    
    assert(str(my_bs) == "00000000")

    my_bs.flip_site(2)
    my_bs.flip_site(7)
    my_bs.flip_site(0)
    assert(str(my_bs) == "10100001")
    assert(len(my_bs) == 8)
    

def test_BS_2():
    my_bs = montecarlo.BitString(13)
    my_bs.set_config([0,1,1,0,0,1,0,0,1,0,1,0,0])
    assert(my_bs.on() == 5)
    assert(my_bs.off() == 8)
    assert(my_bs.int() == 3220)

def test_BS_3():
    my_bs = montecarlo.BitString(20)
    my_bs.set_int_config(3221)

    # Let's make sure this worked:
    tmp = np.array([0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0,1,0,1])
    assert((my_bs.config == tmp).all())

    # We can provide an even stronger test here:
    for i in range(1000):
        my_bs.set_int_config(i) # Converts from integer to binary
    assert(my_bs.int() == i) # Converts back from binary to integer and tests

def test_BS_4():
    my_bs1 = montecarlo.BitString(13)
    my_bs1.set_config([0,1,1,0,0,1,0,1,1,0,1,0,0])

    my_bs2 = montecarlo.BitString(13)
    my_bs2.set_int_config(3252)


    assert(my_bs1 == my_bs2)

    my_bs2.flip_site(5)
    assert(my_bs1 != my_bs2)


def test_average_values():
    N=10
    T = 2.0
    # now test the general ising hamiltonian
    Jval = 1.0
    mu = [.1 for i in range(N)]
    J = []
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])
    ham2 = montecarlo.IsingHamiltonian(J=J, mu=mu)
    E, M, HC, MS = ham2.compute_average_values(T) 
 
    assert(np.isclose(E, -4.6378514858094695))
    assert(np.isclose(M, -0.1838233606011354 ))
    assert(np.isclose(HC, 1.9883833749653714 ))
    assert(np.isclose(MS, 1.8391722085614428))

def test_energy():
    x = [] # Store list of indices
    y = [] # Store list of energies
    xmin = None # configuration of minimum energy configuration
    emin = 0 # minimum of energy
    my_bs = montecarlo.BitString(10)
    my_bs.set_config([0,0,1,0,1,0,0,1,0,1])

    N=10
    T = 2.0
    # now test the general ising hamiltonian
    Jval = 1.0
    mu = [.1 for i in range(N)]
    J = []
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])
    ham2 = montecarlo.IsingHamiltonian(J=J, mu=mu)

    assert(np.isclose(ham2.energy(my_bs), -6.2))


    
if __name__== "__main__":
    test_montecarlo_imported()
    test_average_values()
    test_energy()
    test_BS_1()
    test_BS_2()
    test_BS_3()
    test_BS_4()