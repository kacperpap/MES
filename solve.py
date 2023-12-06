import numpy as np
import scipy.sparse
import scipy.sparse.linalg

class Solve:
    def __init__(self, H, P, C, initialTemperature, simulationStepTime, simulationTime):
        
        self.simulationStepTime = simulationStepTime
        self.simulationTime = simulationTime
        
        #ustalona wymiana ciepla
        self.t_stationary = np.full((len(P), 1), 0)
        
        H_sparse = scipy.sparse.csr_matrix(H)
        C_sparse = scipy.sparse.csr_matrix(C)
        
        #Rozwiazanie [H] * {t} = -{P}
        self.t_stationary = scipy.sparse.linalg.spsolve(H_sparse, P)    
    
    
        #################################################
        #  nieustalona/ niestacjonarna wymiana ciepla   #
        #  ([H] + [C]/dT){t1} - ([C]/dT){t0} + {P} = 0  #
        #                                               #
        #   ([H] + [C]/dT){t1} = ([C]/dT){t0} - {P}     #
        #                                               #
        #################################################
        
        self.t_nonstationary_table = np.empty((int(simulationTime/simulationStepTime), len(P)), dtype=float)
        self.t_nonstationary_end = np.full((len(P)), 0)
        
        t0 = np.full((len(P)), initialTemperature)
        
        A = H_sparse + C_sparse/simulationStepTime
        
        for i in range(int(simulationTime/simulationStepTime)):
            B = (C_sparse/simulationStepTime) @ t0 - (-P)
            self.t_nonstationary_table[i] = scipy.sparse.linalg.spsolve(A,B)
            t0 = self.t_nonstationary_table[i].copy()
        
        self.t_nonstationary_end = self.t_nonstationary_table[int(simulationTime/simulationStepTime) - 1] 
            

        
        
    def print_tvector(self,t_vec, f = None):
        if f is None:
            print('[', end='')
            for i in t_vec:
                if i == 0:
                    print("{:6}".format("0"), end="  ")
                else:
                    print("{:6.3f}".format(i), end="  ")
            print(']')
            
        else:
            print('[', end='',file=f)
            for i in t_vec:
                if i == 0:
                    print("{:6}".format("0"), end="  ",file=f)
                else:
                    print("{:6.3f}".format(i), end="  ",file=f)
            print(']',file=f)
            
            
            
    def print_simulation(self, f=None):
        step_time = self.simulationTime / len(self.t_nonstationary_table)

        for i in range(len(self.t_nonstationary_table)):
            time = (i+1) * step_time

            max_temp = np.max(self.t_nonstationary_table[i])
            min_temp = np.min(self.t_nonstationary_table[i])

            if f is None:
                print(f"Czas: {time:.2f}, Max temp: {max_temp:.2f}, Min temp: {min_temp:.2f}")
            else:
                print(f"Czas: {time:.2f}, Max temp: {max_temp:.2f}, Min temp: {min_temp:.2f}", file=f)

        
        