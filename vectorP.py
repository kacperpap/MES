import numpy as np


class vectorP:
    def __init__(self,integrationPointsNumber,universalElement ,nodes,nodesBounderConditionFlags,alpha, Tambient):
        # self.universalElement2D = ue.UniversalElement2D(integrationPointsNumber)
        W = universalElement.W
        BCis_shapeFuncsTab = universalElement.BCis_shapeFuncsTab
        
        self.vectorsP = np.zeros((4,4), dtype=float)

        for side in range(4):
            if(nodesBounderConditionFlags[side][0] == 1 and nodesBounderConditionFlags[(side + 1) % 4][0] == 1):
                x1, y1 = nodes[side]
                x2, y2 = nodes[(side + 1) % 4]
                L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)                #jeśli krawędzie nie są poziome lub pionowe to ich długość
                detJ = L/2.0                                             #obliczamy jako sqrt((x2 - x1)**2 + (y2 - y1)**2)
                for borderPoint in range (integrationPointsNumber):
                    self.vectorsP[side] += alpha * BCis_shapeFuncsTab[side*integrationPointsNumber + borderPoint] * W[side * integrationPointsNumber + borderPoint] * Tambient
                self.vectorsP[side] *= detJ
            
        self.vectorP = np.zeros((1,4), dtype=float)
        for i in range(4):
            self.vectorP += self.vectorsP[i]
            
    def printP(self):
        np.set_printoptions(formatter={'float': '{: 0.5f}'.format})
        print(np.array2string(self.vectorP,precision=5,suppress_small=True))