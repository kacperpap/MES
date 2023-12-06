import numpy as np

class gaussSchema1D:
    def __init__(self, A, x):
        self.A = A
        self.x = x 

class gaussSchema2D:
    def __init__(self, A, P):
        self.A = A
        self.P = P


def gaussIntegrationSchemaiD(integrationPointsNumber, iDimensions):
    """
    (integrationPointsNumber = 2/3/4; iDimensions = 1/2)
    
    Zwraca dla:
        iDimension = 1:
            A - np tablica przechowujaca wagi punktow siatki o rozmiarze integrationPointsNumber
            x - np tablica przechowujaca wspolrzedne x-owe punktow o rozmiarze integrationPointsNumber
        iDimension = 2:
            A - np tablica przechowujaca wagi punktow o rozmiarze integrationPointsNumber**2 
            P - np tablica przechowujaca wspolrzedne x-owe i y-owe punktow o rozmiarze integrationPointsNumber**2, 
                w postaci krotek [x,y], (aby odwolac sie do wszystkich x-ow w tabeli przechodzimy w petli
                o odpowiednim rozmiarze P[i,0] w kazdej iteracji, dla punktow y-owych P[i,1]
    """

    if(iDimensions == 1):
        A = np.zeros(integrationPointsNumber, dtype=float)
        x = np.zeros(integrationPointsNumber, dtype=float)
    elif(iDimensions == 2):
        A = np.zeros(integrationPointsNumber**2, dtype=float)
        P = np.zeros([integrationPointsNumber**2,2], dtype=float)
    else:
        raise ValueError("\n\tgaussIntegrationSchemaiD: iDimension must be 1 or 2\n")
        
    
    if(iDimensions == 1):
        if integrationPointsNumber == 2:
            x[0] = -pow(1/3,1/2)
            x[1] = -x[0]
            A[0] = A[1] = 1.0
        elif integrationPointsNumber == 3:
            x[0] = 0.0
            x[1] = -pow(3/5,1/2)
            x[2] = -x[1]
            A[0] = 8.0 / 9.0
            A[1] = A[2] = 5.0 / 9.0
        elif integrationPointsNumber == 4:
            x[0] = -pow(3/7 - 2/7 * pow(6/5,1/2),1/2)
            x[1] = -x[0]
            x[2] = -pow(3/7 + 2/7 * pow(6/5,1/2),1/2)
            x[3] = -x[2]
            A[0] = A[1] = (18.0 + pow(30,1/2) ) / 36.0
            A[2] = A[3] = (18.0 - pow(30,1/2) ) / 36.0
        else:
            raise ValueError("\n\tgaussIntegrationSchemaiD: integrationPointsNumber must be in [2,3,4]\n")

        return gaussSchema1D(A,x)
    
    elif(iDimensions == 2):
        if integrationPointsNumber == 2:
            a = -pow(1/3,1/2)
            P[0] = [a, a]; P[1] = [-a, a]; P[2] = [a,-a]; P[3] = [-a,-a]
            A[0] = A[1] = A[2] = A[3] = 1.0
        elif integrationPointsNumber == 3:
            c = -pow(3/5,1/2); d = 8.0 / 9.0; e = 5.0 / 9.0
            P[0] = [c,c]; P[1] = [0,c]; P[2] = [-c,c]; P[3] = [c,0]; P[4] = [0,0]; P[5] = [-c,0]; P[6] = [c,-c]; P[7] = [0,-c]; P[8] = [-c,-c]
            A[0] = A[2] = A[6] = A[8] = e*e
            A[1] = A[3] = A[5] = A[7] = e*d
            A[4] = d*d
        elif integrationPointsNumber == 4:
            f = -pow(3/7 - 2/7 * pow(6/5,1/2),1/2)
            g = -pow(3/7 + 2/7 * pow(6/5,1/2),1/2)
            h = (18.0 + pow(30,1/2) ) / 36.0
            j = (18.0 - pow(30,1/2) ) / 36.0
            P[0] = [g,g]; P[1] = [f,g]; P[2] = [-f,g]; P[3] = [-g,g]; P[4] = [g,f]; P[5] = [f,f]; P[6] = [-f,f]; P[7] = [-g,f]; P[8]= [-g,-f]
            P[9] = [-f,-f]; P[10] = [f,-f]; P[11] = [g,-f]; P[12] = [g,-g]; P[13] = [f,-g]; P[14] = [-f,-g]; P[15] = [-g,-g]
            A[0] = A[3] = A[12] = A[15] = j*j
            A[1] = A[2] = A[4] = A[7] = A[8] = A[11] = A[13] = A[14] = h*j
            A[5] = A[6] = A[9] = A[10] = h*h
        else:
            raise ValueError("\n\tgaussIntegrationSchemaiD: integrationPointsNumber must be in [2,3,4]\n")

        return gaussSchema2D(A,P)

        


def gaussIntegrate(f, integrationPointsSchema, dimension):
    """
    Całkuje w przedziale [-1;1] za pomocą metody Gaussa.

    Argumenty: 
    f -- funkcja (dla 1D funkcja jednej zmiennej f(x), dla 2D funkcja f(x,y))
    integrationPointsSchema -- ilość punktów w schemacie całkowania iD,
        2 oznacza że wykres funkcji na przedziale [-1;1] dzielimy na 2 zakresy i dla
        każdego obliczamy wartość funkcji w punkcie,dla 2D liczba zakresów to integrationPointsSchema**2
    dimension -- określa wymiar 1D - oblicza pole pod wykresem, 2D 
        objętość pod płaszczyzną

    Zwraca wartosc typu float
    """
    result = 0

    schema = gaussIntegrationSchemaiD(integrationPointsSchema, dimension)

    if(dimension == 1):
        for i in range(integrationPointsSchema):
            result += schema.A[i] * f(schema.x[i])
    elif(dimension == 2):
        for i in range(integrationPointsSchema**2):
            result += schema.A[i] * f(schema.P[i,0], schema.P[i,1])

    return result





