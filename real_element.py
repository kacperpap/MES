import universal_element as ue
import numpy as np
from scipy.linalg import lu_factor, lu_solve


class RealElement2D:
    """
    Przeksztalca pochodne funkcji ksztaltu z ukladu lokalnego na uklad globalny
    
    Argumenty:
        integrationPointsNumber -> liczba punktow calkowania
        nodes_ -> liczba krotek reprezentujacych wspolrzedne pojednyczego elementu siatki globalnej,
              dla siatki "kwadratowej" (element siatki to kwadrat z 4 wierzcholkami), funkcja przyjmuje 
              np.array([4,2]) czyli [[x1,y1], [x2,y2], [x3,y3], [x4,y4]]
              UWAGA, nalezy poprawnie okreslic wierzcholki na podstawie danych zawartych w grid.Elements
    """
    def __init__(self,integrationPointsNumber,universalElement,nodes):
        self.universalElement2D = universalElement
        self.dNdKsiTab = universalElement.dNdKsiTab
        self.dNdEtaTab = universalElement.dNdEtaTab
        
        self.dNdXTab = np.zeros((integrationPointsNumber**2,4), dtype=float)
        self.dNdYTab = np.zeros((integrationPointsNumber**2,4), dtype=float)
        self.detJ_iPC = np.zeros((integrationPointsNumber**2,1), dtype=float)

        
        for i in range(integrationPointsNumber**2):                
            #obliczamy jakobian dla KAZDEGO i-tego PUNKTU CALKOWANIA OSOBNO
            self.jacobyObj = ue.JacobyMatrix(nodes,universalElement,i, integrationPointsNumber)
            self.jacobian = self.jacobyObj.jacobian
            
            #############################################
            #   Bardziej optymalna metoda obliczania    #
            #   det oraz inv opiera się o macierz LU    #
            #############################################
            
            lu, piv = lu_factor(self.jacobian)
            identity = np.eye(self.jacobian.shape[0])
            
            # self.detJ_iPC[i,0] = np.linalg.det(self.jacobian)
            # self.scaleMatrix = np.linalg.inv(self.jacobian)
            self.detJ_iPC[i,0] = np.prod(np.diagonal(lu))
            self.scaleMatrix = lu_solve((lu,piv), identity)
            
            for j in range(4):
                self.dNdXTab[i][j] = self.scaleMatrix[0][0]*self.dNdKsiTab[i][j] + self.scaleMatrix[0][1]*self.dNdEtaTab[i][j]
                self.dNdYTab[i][j] = self.scaleMatrix[1][0]*self.dNdKsiTab[i][j] + self.scaleMatrix[1][1]*self.dNdEtaTab[i][j]
                
                