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
        
            
            
            
  def print_HG(self, f = None):
    if f is None:
      for row in self.HG:
          print('[', end='')
          for val in row:
              if val == 0:
                  print("{:7}".format("0"), end="  ")
              else:
                  print("{:7.3f}".format(val), end="  ")
          print(']')
    else:
      for row in self.HG:
          print('[', end='',file=f)
          for val in row:
              if val == 0:
                  print("{:7}".format("0"), end="  ",file=f)
              else:
                  print("{:7.3f}".format(val), end="  ",file=f)
          print(']',file=f)
          
          
          
  def print_PG(self, f = None):
    if f is None:
      print('[', end='')
      for i in self.PG:
        if i == 0:
          print("{:6}".format("0"), end="  ")
        else:
          print("{:6.3f}".format(i), end="  ")
      print(']')
    else:
      print('[', end='',file=f)
      for i in self.PG:
        if i == 0:
          print("{:6}".format("0"), end="  ",file=f)
        else:
          print("{:6.3f}".format(i), end="  ",file=f)
      print(']',file=f)
      
      
      
  def print_CG(self, f = None):
    if f is None:
      for row in self.CG:
          print('[', end='')
          for val in row:
              if val == 0:
                  print("{:6}".format("0"), end="  ")
              else:
                  print("{:6.1f}".format(val), end="  ")
          print(']')
    else:
      for row in self.CG:
          print('[', end='',file=f)
          for val in row:
              if val == 0:
                  print("{:6}".format("0"), end="  ",file=f)
              else:
                  print("{:6.1f}".format(val), end="  ",file=f)
          print(']',file=f)

    
    