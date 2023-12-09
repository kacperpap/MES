import numpy as np


class matrixHBC:
    def __init__(self,integrationPointsNumber,universalElement,nodes,nodesBounderConditionFlags,alpha):
        #Obliczanie macierzy Hbc dla kazdej plaszczyzny (ściany ze wzoru: SUMA(0...n_pkt_na_ścianie)(alfa * (({N}{N}^T * w1) + ({N}{N}^T * w2))))
        W = universalElement.W
        BCis_shapeFuncsTab = universalElement.BCis_shapeFuncsTab
        

        self.HBCs = np.zeros((4,4,4), dtype=float)
        
        
        for side in range(4):
            if(nodesBounderConditionFlags[side][0] == 1 and nodesBounderConditionFlags[(side + 1) % 4][0] == 1):
                x1, y1 = nodes[side]
                x2, y2 = nodes[(side + 1) % 4]
                L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)                #jeśli krawędzie nie są poziome lub pionowe to ich długość
                detJ = L/2.0                                             #obliczamy jako sqrt((x2 - x1)**2 + (y2 - y1)**2)
                for borderPoint in range (integrationPointsNumber):
                    self.HBCs[side] += alpha * np.outer(BCis_shapeFuncsTab[side*integrationPointsNumber + borderPoint], BCis_shapeFuncsTab[side*integrationPointsNumber + borderPoint]) * W[side * integrationPointsNumber + borderPoint]
                self.HBCs[side] *= detJ
            
        self.HBC = np.zeros((4,4),dtype=float)
        for i in range(4):
            self.HBC += self.HBCs[i]
            
    def printHBC(self):
        np.set_printoptions(formatter={'float': '{: 0.5f}'.format})
        print(np.array2string(self.HBC,precision=5,suppress_small=True))
        