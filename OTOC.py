from qiskit import QuantumCircuit, assemble, Aer
import qiskit as qiskit
from qiskit import execute
from qiskit.visualization import plot_histogram, plot_bloch_vector
from math import sqrt, pi
import random
from multiprocessing import Pool, Queue
import multiprocessing
from multiprocessing import Pool
import csv
from multiprocessing import Process
import cmath
import time
import pandas as pd
import time
import numpy as np
import pickle
import os
import sys

# python OTOC.py Nq otoc_type sample_size output_path_folder_name

FINAL_RES = []
FINAL_HAAR = []
FINAL_PAULI = []
FINAL_CLIF = []

Nq = int(sys.argv[1])
otoc_type = str(sys.argv[2])
sample_size = int(sys.argv[3])
folder_name = str(sys.argv[4])

if (folder_name[-1] != "/"):
    folder_name = folder_name + "/" + "Qubit_"+ str(Nq) + "/"
else :
    folder_name = folder_name + "Qubit_"+ str(Nq) + "/"
print(folder_name)    
seeds = np.arange(0,10000,1000)

##############################################################################################################

def OTOC_4_1(ensemble, q):
    """Implementation of Equation 17 (Alves)
    Args:
        ensemble(list) : Contains ensemble of unitaries belonging to 
                        a given design
        Nq (int): no of qubits
       Returns:
       results(list) : returns a list (length : Nq x Nq) --> need to reshape to (Nq,Nq)
    """

    results = []
    for i in range(Nq):
        for j in range(Nq):
            seed = (random.randint(0,100000))
            sample_sum = 0
            for m in range(5):
                seed1 = (random.randint(0,len(ensemble)-1))
                U = ensemble[seed1] # replace count with seed and check results
                qc = qiskit.QuantumCircuit(Nq, Nq)
                vec = qiskit.quantum_info.random_statevector(2**Nq)
                bra = np.matrix(vec.conjugate().data)
                qc.initialize(vec)
                qc.append(U, range(Nq))
                qc.y(j)
                qc.append(U.inverse(), range(Nq))
                qc.x(i)
                qc.append(U, range(Nq))
                qc.y(j)
                qc.append(U.inverse(), range(Nq))
                qc.x(i)
                qc.measure_all()
                job = execute(qc, shots = 100, backend = Aer.get_backend("statevector_simulator"), seed_simulator = seed)
                result = job.result()
                output = result.get_statevector(qc, decimals=5)
                ket = np.matrix(output.data)
                value = np.inner(bra, ket)[0,0] 
                sample_sum = sample_sum + value 
            results.append(sample_sum)
    q.put(results)
#     return results


##############################################################################################################


def OTOC_4_2(ensemble, q):
    """Implementation of Equation 18 (Alves)
    Args:
        ensemble(list) : Contains ensemble of unitaries belonging to 
                        a given design
        Nq (int): no of qubits
        Ns (int) : no of shots
       Returns:
       results(list) : returns a list (length : Nq x Nq) --> need to reshape to (Nq,Nq)
    """
    results = []
    for i in range(0,Nq):
        for j in range(0, Nq):
            seed = (random.randint(0,100000))
            sample_sum = 0
            for m in range(5):
                seed1 = (random.randint(0,len(ensemble)-1))
                U = ensemble[seed1] # replace count with seed and check results
                qc = qiskit.QuantumCircuit( Nq, Nq)
                vec = qiskit.quantum_info.random_statevector(2**Nq)
                bra = np.matrix(vec.conjugate().data)
                
                qc.initialize(vec)
                qc.append(U, range(Nq))
                qc.y(j)
                qc.append(U.inverse(), range(Nq))
                qc.y(i)
                qc.append(U, range(Nq))
                qc.x(j)
                qc.append(U.inverse(), range(Nq))
                qc.x(i)
                qc.measure_all()
                job = execute(qc, shots = 100, backend = Aer.get_backend("statevector_simulator"), seed_simulator = seed)
                result = job.result()
                output = result.get_statevector(qc, decimals=5)
                ket = np.matrix(output.data)
                value = np.inner(bra, ket)[0,0] 
                sample_sum = sample_sum + value 
            results.append(sample_sum)
    q.put(results)
#     return results



################################################################################################################

def OTOC_2_1(ensemble, q):
    """Implementation of Equation 19 (Alves)
    Args:
        ensemble(list) : Contains ensemble of unitaries belonging to 
                        a given design
        Nq (int): no of qubits
        Ns (int) : no of shots
       Returns:
       results(list) : returns a list (length : Nq x Nq) --> need to reshape to (Nq,Nq)
    """
    results = []
    for i in range(0,Nq):
        for j in range(0, Nq):
            seed = (random.randint(0,100000))
            sample_sum = 0
            for m in range(5):
                seed1 = (random.randint(0,len(ensemble)-1))
                U = ensemble[seed1] # replace count with seed and check results
                qc = qiskit.QuantumCircuit(Nq, Nq)
                vec = qiskit.quantum_info.random_statevector(2**Nq)
                bra = np.matrix(vec.conjugate().data)
                qc.initialize(vec)
                qc.append(U, range(Nq))
                qc.y(j)
                qc.append(U.inverse(), range(Nq))
                qc.x(i)
                qc.measure_all()
                job = execute(qc, shots = 100, backend = Aer.get_backend("statevector_simulator"), seed_simulator = seed)
                result = job.result()
                output = result.get_statevector(qc, decimals=5)
                ket = np.matrix(output.data)
                value = np.inner(bra, ket)[0,0] 
                sample_sum = sample_sum + value 
            results.append(sample_sum)
   
    q.put(results)
#     return results



################################################################################################################

def OTOC_2_2(ensemble, q):
    """Implementation of Equation 19 (Alves)
    Args:
        ensemble(list) : Contains ensemble of unitaries belonging to 
                        a given design
        Nq (int): no of qubits
        Ns (int) : no of shots
       Returns:
       results(list) : returns a list (length : Nq x Nq) --> need to reshape to (Nq,Nq)
    """
    results = []
    for i in range(0,Nq):
        for j in range(0, Nq):
            seed = (random.randint(0,100000))
            sample_sum = 0
            for m in range(5):
                seed1 = (random.randint(0,len(ensemble)-1))
                U = ensemble[seed1] # replace count with seed and check results
                qc = qiskit.QuantumCircuit( Nq, Nq)
                vec = qiskit.quantum_info.random_statevector(2**Nq)
                bra = np.matrix(vec.conjugate().data)
                qc.initialize(vec)
                qc.append(U, range(Nq))
                qc.z(j)
                qc.append(U.inverse(), range(Nq))
                qc.z(i)
                qc.measure_all()
                job = execute(qc, shots = 100, backend = Aer.get_backend("statevector_simulator"), seed_simulator = seed)
                result = job.result()
                output = result.get_statevector(qc, decimals=5)
                ket = np.matrix(output.data)
                value = np.inner(bra, ket)[0,0] 
                sample_sum = sample_sum + value 
            results.append(sample_sum)

    q.put(results)
#     return results


#############################################################################################################
def load_data():
        """Load ensemble data 
    Args:
       Returns:
       results(list) : returns a list containing clifford, haar and pauli
    """

    
    clif1 = []
    haar1 = []
    pauli1 = []
    PIK = folder_name + "clifford_ens_"+ str(Nq) +".dat"
    with open(PIK, "rb") as fp:   # Unpickling
        clif1 = pickle.load(fp)
    
    PIK = folder_name + "haar_ens_"+ str(Nq) +".dat"
    with open(PIK, "rb") as fp:   # Unpickling
        haar1 = pickle.load(fp)
    
    PIK = folder_name + "pauli_ens_"+ str(Nq) +".dat"
    with open(PIK, "rb") as fp:   # Unpickling
        pauli1 = pickle.load(fp)
    
    return [clif1,haar1,pauli1]

#############################################################################################################


if __name__ == '__main__':
    
    HAAR = []
    PAULI = []
    CLIF = []
    otoc_dict = {   '1':OTOC_4_1,
                    '2':OTOC_4_2,
                    '3':OTOC_2_1,
                    '4':OTOC_2_2
                }
    
    CLIF, HAAR, PAULI = load_data()
    
    q1 = Queue()
    q2 = Queue()
    q3 = Queue()
    
    
    for i in range(sample_size):
        target = otoc_dict[otoc_type]
        p1 = Process(target=target, args=(HAAR, q1,))
        p1.start()
        p2 = Process(target=target, args=(PAULI, q2,))
        p2.start()
        p3 = Process(target=target, args=(CLIF, q3,))
        p3.start()

    for i in range(sample_size):
        FINAL_HAAR.append(q1.get())
        FINAL_PAULI.append(q2.get())
        FINAL_CLIFF.append(q3.get())
     
    PIK = folder_name + "haar_" + otoc_choice + ".dat"
    outfile = open(PIK,'wb')
    with open(PIK, "wb") as fp:
        pickle.dump(FINAL_HAAR, fp)
    
    PIK = folder_name + "pauli_" + otoc_choice + ".dat"
    outfile = open(PIK,'wb')
    with open(PIK, "wb") as fp:
        pickle.dump(FINAL_PAULI, fp)

    PIK = folder_name + "clif_" + otoc_choice + ".dat"
    outfile = open(PIK,'wb')
    with open(PIK, "wb") as fp:
        pickle.dump(FINAL_CLIF, fp)
