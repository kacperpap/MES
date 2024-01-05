import numpy as np



class Agregation:
  def __init__(self, N):
    self.N = N
    self.HG = np.zeros((N,N), dtype=float)
    self.PG = np.zeros((N),dtype=float)
    self.CG = np.zeros((N,N), dtype=float)
   
    
  def aggregate(self, H, HBC, P, C, elementIDs):
    H_HBC = H + HBC
    
    n = len(elementIDs)
    self.IDsFromLocalToGLobalH_HBC = np.zeros((n,n,2),dtype=int)
    self.IDsFromLocalToGLobalP = np.zeros((n),dtype=int)
    
    #konwertujemy indeksy wezlow w elemencie na globalna tablice wszystkich wezlow
    for i in range(n):
      for j in range(n):
        self.IDsFromLocalToGLobalH_HBC[i][j][0] = elementIDs[i]
        self.IDsFromLocalToGLobalH_HBC[i][j][1] = elementIDs[j]
      
    for side in range(n):
      self.IDsFromLocalToGLobalP[side] = elementIDs[side]
        
    for i in range(n):
      for j in range(n):
        self.HG[self.IDsFromLocalToGLobalH_HBC[i][j][0], self.IDsFromLocalToGLobalH_HBC[i][j][1]] += H_HBC[i,j]
            
    for side in range(n):
      self.PG[self.IDsFromLocalToGLobalP[side]] += P[0][side]
      
    for i in range(n):
      for j in range(n):
        self.CG[self.IDsFromLocalToGLobalH_HBC[i][j][0], self.IDsFromLocalToGLobalH_HBC[i][j][1]] += C[i,j]
        
            
            
            
  def print_matrix(self, matrix, format_str, f=None):
    assert isinstance(matrix, np.ndarray), "Expected a numpy array"
    if matrix.ndim == 1:
        print('[', end='', file=f)
        for val in matrix:
            if val == 0:
                print("{:6}".format("0"), end="  ", file=f)
            else:
                print(format_str.format(val), end="  ", file=f)
        print(']', file=f)
    else:
        for row in matrix:
            print('[', end='', file=f)
            for val in row:
                if val == 0:
                    print("{:6}".format("0"), end="  ", file=f)
                else:
                    print(format_str.format(val), end="  ", file=f)
            print(']', file=f)

  def print_HG(self, f=None):
    self.print_matrix(self.HG, "{:6.1f}", f)

  def print_PG(self, f=None):
    self.print_matrix(self.PG, "{:6.3f}", f)

  def print_CG(self, f=None):
    self.print_matrix(self.CG, "{:6.1f}", f)
    
    