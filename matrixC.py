import numpy as np
import real_element


class matrixC:
    def __init__(self,integrationPointsNumber,nodes ,specificHeat, density):
        
        self.realElement = real_element.RealElement2D(integrationPointsNumber,nodes)
        self.NTab = self.realElement.universalElement2D.shapeFunTab
        self.dV = self.realElement.detJ_iPC
        
        self.CMatrixesTable = np.zeros((integrationPointsNumber**2,4,4),dtype=float)
        
        for i in range(integrationPointsNumber**2):
            self.CMatrixesTable[i] = specificHeat * density * np.outer(self.NTab[i], self.NTab[i]) * self.dV[i,0]
            
            
        self.gaussWages = self.realElement.universalElement2D.schema2D.A
        
        self.C = np.zeros((4,4),dtype=float)
        
        for i in range(integrationPointsNumber**2):
            self.C += self.CMatrixesTable[i] * self.gaussWages[i]


