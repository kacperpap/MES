import gauss_integration as gauss
import numpy as np

class ShapeFunctions2D:

    def __init__(self):
        self.dNdKsiFun = [
            self.dN1_dKsi,
            self.dN2_dKsi,
            self.dN3_dKsi,
            self.dN4_dKsi
        ]

        self.dNdEtaFun = [
            self.dN1_dEta,
            self.dN2_dEta,
            self.dN3_dEta,
            self.dN4_dEta
        ]
        
        self.NFun = [
            self.N1,
            self.N2,
            self.N3,
            self.N4
        ]


    def dN1_dKsi(self, eta):
        return (-1.0 / 4.0) * (1.0 - eta)
    def dN2_dKsi(self, eta):
        return (1.0 / 4.0) * (1.0 - eta)
    def dN3_dKsi(self, eta):
        return (1.0 / 4.0) * (1.0 + eta)
    def dN4_dKsi(self, eta):
        return (-1.0 / 4.0) * (1.0 + eta)
    def dN1_dEta(self, ksi):
        return (-1.0 / 4.0) * (1.0 - ksi)
    def dN2_dEta(self, ksi):
        return (-1.0 / 4.0) * (1.0 + ksi)
    def dN3_dEta(self, ksi):
        return (1.0 / 4.0) * (1.0 + ksi)
    def dN4_dEta(self, ksi):
        return (1.0 / 4.0) * (1.0 - ksi) 
    def N1(self,ksi,eta):
        return 1.0/4.0 * (1 - ksi) * (1 - eta)
    def N2(self,ksi,eta):
        return 1.0/4.0 * (1 + ksi) * (1 - eta)
    def N3(self,ksi,eta):
        return 1.0/4.0 * (1 + ksi) * (1 + eta)
    def N4(self,ksi,eta):
        return 1.0/4.0 * (1 - ksi) * (1 + eta)



class UniversalElement2D:
    """
    UniversalElement2D, to element w ukladzie lokalnym w przestrzeni 2D
    
    Tworzy dwie tablice:
        dNdKsiTab = np.zeros((integrationPointsNumber**2,4))
        dNdEtaTab = np.zeros((integrationPointsNumber**2,4))
        
    Zawieraja one pochodne funkcji ksztaltu po osi ksi oraz eta, co jest konieczne do obliczenia
    pochodnych funkcji ksztaltu po x i y (czyli w globalnym ukladzie wspolrzednych) ze wzoru:
    | dNi/dKsi |   | dx/dksi  dy/dksi |   | dNi/dx |
    |          | = |                  | * |        |
    | dNi/dEta |   | dx/deta  dy/deta |   | dNi/dy |
    """
    def __init__(self, integrationPointsNumber):
        self.integrationPointsNumber = integrationPointsNumber
        self.schema2D = gauss.gaussIntegrationSchemaiD(integrationPointsNumber,2)
        P = self.schema2D.P
        
        #zdefiniowanie punktow brzegowych (wspolrzednych dla konkretnych punktow brzegowych)
        #zakladamy ze kazdy punkt to krotka B[numer_ściany,punkt_na_ścianie,x,y]
        #ściany definiujemy jako 0...3 od dołu przeciwnie do wskazówek zegara
        #UWAGA dwie pierwsze wartosci to int dwie nastepne to float
        self.B = np.zeros([4*integrationPointsNumber,2])
        self.W = np.zeros([4*integrationPointsNumber,1])

        
        if integrationPointsNumber == 2:
            self.B[0] = [P[0,0],-1]; self.B[1] = [P[1,0],-1]
            self.B[2] = [1,P[1,1]]; self.B[3] = [1,P[3,1]]
            self.B[4] = [P[3,0],1]; self.B[5] = [P[2,0],1]
            self.B[6] = [-1,P[2,1]]; self.B[7] = [-1,P[0,1]]
            self.W[0] = self.W[1] = self.W[2] = self.W[3] = self.W[4] = self.W[5] = self.W[6] = self.W[7] = 1.0 
        elif integrationPointsNumber == 3:
            w1 = w3 = 5.0/9.0; w2 = 8.0/9.0
            self.B[0] = [P[0,0],-1]; self.B[1] = [P[1,0],-1]; self.B[2] = [P[2,0],-1]
            self.B[3] = [1,P[2,1]]; self.B[4] = [1,P[5,1]]; self.B[5] = [1,P[8,1]]
            self.B[6] = [P[8,0],1]; self.B[7] = [P[7,0],1]; self.B[8] = [P[6,0],1]
            self.B[9] = [-1,P[6,1]]; self.B[10] = [-1,P[3,1]]; self.B[11] = [-1,P[0,1]]
            self.W[0] = self.W[3] = self.W[8] = self.W[11] = w1
            self.W[1] = self.W[4] = self.W[7] = self.W[10] = w2
            self.W[2] = self.W[5] = self.W[6] = self.W[9] = w3
        elif integrationPointsNumber == 4:
            w1 = (18.0 + pow(30,1/2) ) / 36.0
            w2 = (18.0 - pow(30,1/2) ) / 36.0
            self.B[0] = [P[0,0],-1]; self.B[1] = [P[1,0],-1]; self.B[2] = [P[2,0],-1]; self.B[3] = [P[3,0],-1]
            self.B[4] = [1,P[3,1]]; self.B[5] = [1,P[7,1]]; self.B[6] = [1,P[11,1]]; self.B[7] = [1,P[15,1]]
            self.B[8] = [P[15,0],1]; self.B[9] = [P[14,0],1]; self.B[10] = [P[13,0],1]; self.B[11] = [P[12,0],1]
            self.B[12] = [-1,P[12,1]]; self.B[13] = [-1,P[8,1]]; self.B[14] = [-1,P[4,1]]; self.B[15] = [-1,P[0,1]]
            self.W[1] = self.W[2] = self.W[5] = self.W[6] = self.W[9] = self.W[10] = self.W[13] = self.W[14] = w1
            self.W[0] = self.W[3] = self.W[4] = self.W[7] = self.W[8] = self.W[11] = self.W[12] = self.W[15] = w2


        #pierwszy indesks np.zeros(x,y) to wiersze drugi to kolumny
        self.dNdKsiTab = np.zeros((integrationPointsNumber**2,4), dtype=float)
        self.dNdEtaTab = np.zeros((integrationPointsNumber**2,4), dtype=float)

        self.shapeFun = ShapeFunctions2D()

        #punkty zapisane jako krotki P[0][x,y] P[1][x,y] ... P[n-1][x,y]
        for i in range(integrationPointsNumber**2):
            for j in range(4):
                self.dNdKsiTab[i][j] = self.shapeFun.dNdKsiFun[j](self.schema2D.P[i,1])
                self.dNdEtaTab[i][j] = self.shapeFun.dNdEtaFun[j](self.schema2D.P[i,0])
                
        #obliczanie macierzy Hbc dla wszystkich punktow B[i,j]
        self.BCis_shapeFuncsTab = np.zeros((4*integrationPointsNumber,4),dtype=float)
        
        for i in range(4):
            for j in range(integrationPointsNumber):
                self.BCis_shapeFuncsTab[i*integrationPointsNumber+j, 0] = self.shapeFun.N1(self.B[i*integrationPointsNumber+j,0],self.B[i*integrationPointsNumber+j,1])
                self.BCis_shapeFuncsTab[i*integrationPointsNumber+j, 1] = self.shapeFun.N2(self.B[i*integrationPointsNumber+j,0],self.B[i*integrationPointsNumber+j,1])
                self.BCis_shapeFuncsTab[i*integrationPointsNumber+j, 2] = self.shapeFun.N3(self.B[i*integrationPointsNumber+j,0],self.B[i*integrationPointsNumber+j,1])
                self.BCis_shapeFuncsTab[i*integrationPointsNumber+j, 3] = self.shapeFun.N4(self.B[i*integrationPointsNumber+j,0],self.B[i*integrationPointsNumber+j,1])
           
            
        self.shapeFunTab = np.zeros((integrationPointsNumber**2,4),dtype=float)
        
        for i in range(integrationPointsNumber**2):
            for j in range(4):
                self.shapeFunTab[i][j] = self.shapeFun.NFun[j](self.schema2D.P[i,0], self.schema2D.P[i,1])



class JacobyMatrix:
    """
    Jakobian dla wymiaru 1D:
    |dx/dksi|
    
    Jakobian dla wymiaru 2D dla i-tego punktu calkowania:
    | dx/dksi  dy/dksi |
    | dx/deta  dy/deta |
    
    UWAGA: 
    dx/dksi = dN_1/dksi * x1 + ... + dN_4/dksi * x4, gdzie [x1,...,x4] to kolejne x-owe wartości wierzchołków naszego
    elementu całkowania (np. kwadratu), a dN_i/dksi oraz dN_i/deta są różne dla każdego punktu całkowania
    (tabela dNdKsiTab to tabela 4x4 poniewaz)
    
    Parametry:
    nodes_ -> liczba krotek reprezentujacych wspolrzedne pojednyczego elementu siatki globalnej,
              dla siatki "kwadratowej" (element siatki to kwadrat z 4 wierzcholkami), funkcja przyjmuje 
              np.array([4,2]) czyli [[x1,y1], [x2,y2], [x3,y3], [x4,y4]]
              UWAGA, nalezy poprawnie okreslic wierzcholki na podstawie danych zawartych w grid.Elements
    (integrationPointsNumber) -> liczba ta nie ma wplywu na jakobian, poniewaz jest on niezalezny
                                 od liczby punktow calkowania DLA TEGO SAMEGO ELEMENTU SIATKI
                                 (jakobian mowi o "skali" przejscia z elementu lokalnego do globalnego
                                 dlatego zalezy tylko od tego elementu)
    iPC -> punkt całkowania dla którego obliczamy jakobian 
    Zwraca:
    jacobian -> macierz o rozmiarach wymiar_przestrzeni x wymiar_przestrzeni (dla 2D to 2x2)
    detJ -> wyznacznik macierzy jacobian
    """
    def __init__(self, nodes_,iPC, integrationPointsNumber):
        if not (isinstance(nodes_, np.ndarray) and (nodes_.shape ==(4,2))):
            raise TypeError("\n\tJacobyMatrix: nodes_ must be type np.array((4,2))")
        
        self.nodes = nodes_
        self.shapeDerivatives = UniversalElement2D(integrationPointsNumber)
        self.jacobian = np.zeros((2,2), dtype=float)
        
        self.dXdKsi = 0; self.dYdKsi = 0; self.dXdEta = 0; self.dYdEta = 0
             
        #tylko dla punktu calkowania ipc                
        for i in range(4):
            self.dXdKsi += self.nodes[i,0] * self.shapeDerivatives.dNdKsiTab[iPC][i]
            self.dYdKsi += self.nodes[i,1] * self.shapeDerivatives.dNdKsiTab[iPC][i]
            self.dXdEta += self.nodes[i,0] * self.shapeDerivatives.dNdEtaTab[iPC][i]
            self.dYdEta += self.nodes[i,1] * self.shapeDerivatives.dNdEtaTab[iPC][i]
                

        self.jacobian = np.array([[self.dXdKsi, self.dYdKsi], [self.dXdEta, self.dYdEta]])
            
        
    
        
            

