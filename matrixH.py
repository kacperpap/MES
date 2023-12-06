import numpy as np
import real_element


                

class matrixH:
    """
    
    Argumenty:
        integrationPointsNumber -> liczba punktow calkowania
        nodes_ -> liczba krotek reprezentujacych wspolrzedne pojednyczego elementu siatki globalnej,
              dla siatki "kwadratowej" (element siatki to kwadrat z 4 wierzcholkami), funkcja przyjmuje 
              np.array([4,2]) czyli [[x1,y1], [x2,y2], [x3,y3], [x4,y4]]
              UWAGA, nalezy poprawnie okreslic wierzcholki na podstawie danych zawartych w grid.Elements
        heatTransferCoeff -> wspolczynnik wymiany ciepla (konwekcyjno- radiacyjny)
    """
    def __init__(self,integrationPointsNumber,nodes, heatTransferCoeff = 30):
        self.realElement = real_element.RealElement2D(integrationPointsNumber,nodes)
        self.dNdXTab = self.realElement.dNdXTab
        self.dNdYTab = self.realElement.dNdYTab
        self.dV = self.realElement.detJ_iPC
        #tablica3D, dla kazdego punktu calkowania mamy macierz 4x4 
        self.HiMatrixesTable = np.zeros((integrationPointsNumber**2,4,4),dtype=float)
        
        
        for i in range(integrationPointsNumber**2):
            self.HiMatrixesTable[i] = heatTransferCoeff * (np.outer(self.dNdXTab[i],self.dNdXTab[i]) + np.outer(self.dNdYTab[i], self.dNdYTab[i])) * self.dV[i,0]
            
            
        self.gaussWages = self.realElement.universalElement2D.schema2D.A
        self.H = np.zeros((4,4),dtype=float)
        for i in range(integrationPointsNumber**2):
            self.H += self.HiMatrixesTable[i] * self.gaussWages[i]
            

    