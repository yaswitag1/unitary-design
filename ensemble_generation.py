from qiskit import QuantumCircuit, assemble, Aer
import qiskit as qiskit
from qiskit import execute
from qiskit.visualization import plot_histogram, plot_bloch_vector
from math import sqrt, pi
import random
from multiprocessing import Pool
import csv
import cmath
import time
import pandas as pd
import time
import numpy as np
import pickle
import os
from multiprocessing import Process
import multiprocessing as mp
from pathlib import Path
from multiprocessing import Pool, Queue
import sys

# python ensemble_generation_v2.py Nq Nu output_path_folder_name

FINAL_RES = []
FINAL_HAAR = []
FINAL_PAULI = []
FINAL_CLIF = []

Nq = int(sys.argv[1])
Nu = int(sys.argv[2])
folder_name = str(sys.argv[3])

if (folder_name[-1] != "/"):
    folder_name = folder_name + "/" + "Qubit_"+ str(Nq) + "/"
else :
    folder_name = folder_name + "Qubit_"+ str(Nq) + "/"

Path(folder_name).mkdir(parents=True, exist_ok=True)
    
print(folder_name)    
print("Nq : " + str(Nq))
print("Nu : " + str(Nu))
print(folder_name)

def gen_clifford_ensemble(seed, Q1):
    ensemble = []
    for i in range(Nu): 
        ensemble.append(qiskit.quantum_info.random_clifford(Nq, seed = seed).to_instruction())
        seed = seed + 1
    Q1.put(ensemble)
#     return ensemble

def gen_haar_ensemble(seed, Q2):
    ensemble = []
    for i in range(Nu):
        ensemble.append(qiskit.quantum_info.random_unitary(2**Nq, seed = seed).to_instruction())
        seed = seed + 1
    Q2.put(ensemble)
#     return ensemble

def gen_pauli_ensemble(seed, Q3):
    ensemble = []
    for i in range(Nu):
        ensemble.append(qiskit.quantum_info.random_pauli(Nq, seed = seed).to_instruction())
        seed = seed + 1
    Q3.put(ensemble)
#     return ensemble

if __name__ == '__main__':

    Q1 = Queue()
    Q2 = Queue()
    Q3 = Queue()
    
    seed = (random.randint(0,100000))
    p1 = Process(target=gen_clifford_ensemble, args=(seed,Q1))
    p1.start()
    
    seed = (random.randint(0,100000))
    p2 = Process(target=gen_haar_ensemble, args=(seed,Q2))
    p2.start()

    seed = (random.randint(0,100000))
    p3 = Process(target=gen_haar_ensemble, args=(seed,Q3))
    p3.start()
    
    p1.join()
    p2.join()
    p3.join()

    CLIFF = Q1.get()
    HAAR = Q2.get()
    PAULI = Q3.get()
    
    PIK = folder_name+"clifford_ens_"+str(Nq)+".dat"
    with open(PIK, "wb") as fp:
        pickle.dump(CLIFF, fp)
        
    PIK = folder_name+"haar_ens_"+str(Nq)+".dat"
    with open(PIK, "wb") as fp:
        pickle.dump(HAAR, fp)
        
    PIK = folder_name+"pauli_ens_"+str(Nq)+".dat"
    with open(PIK, "wb") as fp:
        pickle.dump(PAULI, fp)

        
