import numpy as np


class GlobalData:
    """
    UWAGA:
    Do poprawnego zaladowania pliku tekstowego konieczna zmiana pliku
    Zmienne NodesNumber i ElementsNumber musza byc pisane bez spacji
    """
    def __init__(self, filename = None):
        self.SimulationTime: int = None
        self.SimulationStepTime: int = None
        self.Conductivity: int = None
        self.Alfa: int = None
        self.Tot: int = None
        self.InitialTemp: int = None 
        self.Density: int = None
        self.SpecificHeat: int = None
        self.NodesNumber: int = None
        self.ElementsNumber: int = None

        if filename:
            self.loadDataFromFile(filename)
        else:
            raise ValueError("\n\tGlobalData: filename must be specified correctly\n")
        
    def loadDataFromFile(self, filename: str):
        with open(filename, "r") as file:
            lines = file.readlines()

            startLine = 0
            endLine = None

            for line in lines:
                if line.strip() == "*Node":
                    endLine = lines.index(line)
                    break
                else:
                    endLine = len(lines)
            

            file.seek(0) #powrot wskaznika na poczatek pliku
            data = {}

            for line in lines[startLine:endLine]:
                parts = line.strip().split()
                if len(parts) >= 2:
                    key, value = parts[0], parts[1]
                    data[key] = value


        self.SimulationTime = int(data.get("SimulationTime", None))
        self.SimulationStepTime = int(data.get("SimulationStepTime", None))
        self.Conductivity = int(data.get("Conductivity", None))
        self.Alfa = int(data.get("Alfa", None))
        self.Tot = int(data.get("Tot", None))
        self.InitialTemp = int(data.get("InitialTemp", None))
        self.Density = int(data.get("Density", None))
        self.SpecificHeat = int(data.get("SpecificHeat", None))
        self.NodesNumber = int(data.get("NodesNumber", None))
        self.ElementsNumber = int(data.get("ElementsNumber", None))

    def printGlobalData(self,f):
        if f is not None:
            print("\n\tGLOBAL DATA:\n",file=f)
            print(f"SimulationTime: {self.SimulationTime}",file=f)
            print(f"SimulationStepTime: {self.SimulationStepTime}",file=f)
            print(f"Conductivity: {self.Conductivity}",file=f)
            print(f"Alfa: {self.Alfa}",file=f)
            print(f"Tot: {self.Tot}",file=f)
            print(f"InitialTemp: {self.InitialTemp}",file=f)
            print(f"Density: {self.Density}",file=f)
            print(f"SpecificHeat: {self.SpecificHeat}",file=f)
            print(f"NodesNumber: {self.NodesNumber}",file=f)
            print(f"ElementsNumber: {self.ElementsNumber}",file=f)
        else:
            print("\n\tGLOBAL DATA:\n")
            print(f"SimulationTime: {self.SimulationTime}")
            print(f"SimulationStepTime: {self.SimulationStepTime}")
            print(f"Conductivity: {self.Conductivity}")
            print(f"Alfa: {self.Alfa}")
            print(f"Tot: {self.Tot}")
            print(f"InitialTemp: {self.InitialTemp}")
            print(f"Density: {self.Density}")
            print(f"SpecificHeat: {self.SpecificHeat}")
            print(f"NodesNumber: {self.NodesNumber}")
            print(f"ElementsNumber: {self.ElementsNumber}") 



class Node:
    def __init__(self, nodes):
        self.XY = np.empty([nodes,2], dtype=float)
        self.nodeBC = np.zeros([nodes,1], dtype=int)
    
    def printNodes(self):
        for i in range(len(self.XY)):
            print('Node {:<4}: [{:<10} {:<10}], BC: {}'.format(i+1, round(self.XY[i][0], 6), round(self.XY[i][1], 6), self.nodeBC[i]))
              
    
class Element:
    def __init__(self, elements):
        self.ID = np.empty([elements,4], dtype=int)
    
    def printElements(self):
        print("\n\tELEMENTS:\n")
        for element in self.ID:
            print(element)
            

class BorderConditions:
    def __init__(self):
        self.BCs = np.array([],dtype=int)
        
    def printBCs(self):
        print("\n\tBORDER CONDITIONS:\n")
        print(', '.join(map(str,self.BCs)))

    
class Grid:
    """
    Umozliwia zaladowanie odpowiedniego pliku tekstowego z siatka
    Zawiera struktury GlobalData, Nodes, Elements
    Filename: obowiazkowy parametr
    
    Nodes: struktura np.array([nodesNumber,2]) -> krotki o strukturze (x,y)
    Elements: struktura np.array([elementsNumber,4]) -> krotki o strukturze (a,b,c,d)
    
    Aby odwolac sie do konkretnej krotki:
        Grid obj = Grid(filename)
        obj.Nodes.XY[i] lub obj.Elements.ID[i]
    
    Aby odowlac sie do wartosci krotki:
        Grid obj = Grid(filename)
        for i in range(obj.GlobalData.NodesNumber):
            obj.Nodes.XY[i,0]   //iterowanie po x-ach
            obj.Nodes.XY[i,1]   //iterowanie po y-ach
            
        for i in range(obj.GlobalData.ElementsNumber):
            obj.Elements.ID[i,0]   //iterowanie po indeksach a
            obj.Elements.ID[i,1]   //iterowanie po indeksach b
            obj.Elements.ID[i,2]   //iterowanie po indeksach c
            obj.Elements.ID[i,3]   //iterowanie po indeksach d
            
    Funkcja getElement(elementId):
        Argumenty:
            elementId - indeks krotki przechowujacej informacje o indeksach wierzcholkow
                        elementu siatki (kwadratu) w tablicy Nodes
        Zwraca:
            nodes = np.array([4,2]), czyli wierzcholki kwadratu o indeksie elementId

    """
    def __init__(self, filename = None):
        self.GlobalData = GlobalData(filename)
        self.Elements = Element(int(self.GlobalData.ElementsNumber))
        self.Nodes = Node(int(self.GlobalData.NodesNumber))
        self.BC = BorderConditions()
        
        if filename:
            self.loadDataFromFile(filename)
    
    def loadDataFromFile(self, filename: str):
        with open(filename, "r") as file:
            lines = file.readlines()

            startLineNodes = 0
            endLineNodes = len(lines)

            #okresl zakres nodes

            for line in lines:
                if line.strip() == "*Node":
                    startLineNodes = lines.index(line)
                    endLineNodes = startLineNodes + int(self.GlobalData.NodesNumber) + 1
                    break
            
            #okresl zakres elements

            file.seek(0)

            startLineElements = 0
            endLineElements = len(lines)

            for line in lines:
                parts = line.strip().split(",")
                if parts[0] == "*Element":
                    startLineElements = lines.index(line)
                    endLineElements = startLineElements + int(self.GlobalData.ElementsNumber) + 1
                    break
                
            #okresl zakres BC's
            
            file.seek(0)
            
            startLineBCs = endLineElements
            endLineBCs = len(lines)
            
            for line in lines:
                if line.strip() == "*BC":
                    startLineBCs = lines.index(line)
                    break


            #pobierz dane
            file.seek(0)

            for line in lines[startLineNodes:endLineNodes]:
                parts = line.strip().split(',')
                if len(parts) >= 3:
                    for part in parts:
                        epsilon = 10**-20
                        if float(part) < epsilon:
                            part = 0.0

                    key = int(parts[0]) - 1           #indeksujemy od zera (w pliku od 1)   
                    value1 = float(parts[1]) 
                    value2 = float(parts[2])

                    self.Nodes.XY[key] = [value1,value2]                
  
        
            file.seek(0)

            for line in lines[startLineElements:endLineElements]:
                parts = line.strip().split(',') 
                if len(parts) >= 5:
                    key = int(parts[0]) - 1               
                    value1 = int(parts[1])   
                    value2 = int(parts[2])    
                    value3 = int(parts[3])    
                    value4 = int(parts[4])    
                                                       
                    self.Elements.ID[key] = [value1, value2, value3, value4]
                        
            
            file.seek(0)
            
            for line in lines[startLineBCs:endLineBCs]:
                parts = line.strip().split(',')
                if len(parts) >= 1 and parts[0] != "*BC":
                    self.BC.BCs = np.append(self.BC.BCs, [int(part) for part in parts])
                    for part in parts:
                        self.Nodes.nodeBC[int(part)-1]=1
                            
    def getElement(self, elementId):
        if(elementId >= int(self.GlobalData.ElementsNumber)):
            raise IndexError("\n\tGrid: getElement: elementId is out of table") 
        
        nodes = np.empty([4,2], dtype=float)
        nodesBCs = np.empty([4,1], dtype=int)
        for i in range(4):
            # nodes[i] = self.Nodes.XY[self.Elements.ID[elementId,i] - 1]
            nodeId = self.Elements.ID[elementId,i] - 1
            nodes[i] = self.Nodes.XY[nodeId]
            nodesBCs[i] = self.Nodes.nodeBC[nodeId]

        return [nodes, nodesBCs]
    
    
    def getElementIDs(self, elementId):
        if(elementId >= int(self.GlobalData.ElementsNumber)):
            raise IndexError("\n\tGrid: getElement: elementId is out of table") 

        nodeIDs = np.empty([4], dtype=int)
        for i in range(4):
            nodeId = self.Elements.ID[elementId,i] - 1
            nodeIDs[i] = nodeId

        return nodeIDs


    

    