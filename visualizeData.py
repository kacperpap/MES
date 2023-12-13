import os


class ParaView:
    def __init__(self, selectedGrid, nodesTemperatureDataSimulationTab):
        self.grid = selectedGrid
        self.nodesNum = selectedGrid.GlobalData.NodesNumber 
        self.elementsNum = selectedGrid.GlobalData.ElementsNumber
        self.dataTab = nodesTemperatureDataSimulationTab
        
        #utworz pliki
        self.create_vtk_simulation_files()
    
    def create_vtk_file(self,filename, iteration):
        with open(filename, 'w') as f:
            f.write("# vtk DataFile Version 2.0\n")
            f.write("Unstructured Grid Example\n")
            f.write("ASCII\n")
            f.write("DATASET UNSTRUCTURED_GRID\n\n")
            f.write(f"POINTS {self.nodesNum} float\n")
            
            for point in self.grid.Nodes.XY:
                f.write(f"{point[0]} {point[1]} 0\n")   #trzeci wymiar z == 0
            
            f.write("\n")
            f.write(f"CELLS {self.elementsNum} {self.elementsNum * 5}\n")
            
            for cell in self.grid.Elements.ID:
                f.write(f"4 {cell[0] - 1} {cell[1] - 1} {cell[2] - 1} {cell[3] - 1}\n")
            
            f.write("\n")
            f.write("CELL_TYPES 9\n")
            
            for i in range(self.elementsNum):
                f.write("9\n")
            
            f.write("\n")
            f.write(f"POINT_DATA {self.nodesNum}\n")         
            f.write("SCALARS Temp float 1\n")
            f.write("LOOKUP_TABLE default\n")
            
            for point in self.dataTab[iteration]:
                f.write(f"{point}\n")


    def create_vtk_simulation_files(self):
        directory = "simulation"
        if not os.path.exists(directory):
            os.makedirs(directory)
        for i in range(len(self.dataTab)):
            filename = os.path.join(directory, f"step{i+1}.vtk")
            self.create_vtk_file(filename, i)