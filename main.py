import grid
import gauss_integration
import universal_element
import real_element
import matrixH
import matrixHBC
import numpy as np
import argparse
import vectorP
import aggregation
import solve
import matrixC
import cProfile
import pstats
import visualizeData

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", help="Wyświetl kroki rozwiązania", action="store_true")
parser.add_argument("-p", "--profile", help="Wyświetl statystyki efektywności kodu", action="store_true")


args = parser.parse_args()



# selectedGrid = grid.Grid("Test1_4_4.txt")
selectedGrid = grid.Grid("Test2_4_4_MixGrid.txt")
# selectedGrid = grid.Grid("Test3_31_31_kwadrat.txt")
# selectedGrid = grid.Grid("Test4_31_31_trapez.txt")

integrationPoints = 4


def main():
    agreg = aggregation.Agregation(selectedGrid.GlobalData.NodesNumber)
    ue = universal_element.UniversalElement2D(integrationPoints)
        
    for i in range(selectedGrid.GlobalData.ElementsNumber):
        nodesData = selectedGrid.getElement(i)
        nodesGlobalIDs = selectedGrid.getElementIDs(i)
        realElementInstance = real_element.RealElement2D(integrationPoints,ue,nodesData[0])
        Hinstance = matrixH.matrixH(integrationPoints,realElementInstance,selectedGrid.GlobalData.Conductivity)
        HBCinstance = matrixHBC.matrixHBC(integrationPoints,ue,nodesData[0],nodesData[1],selectedGrid.GlobalData.Alfa)
        Pinstance = vectorP.vectorP(integrationPoints,ue,nodesData[0],nodesData[1],selectedGrid.GlobalData.Alfa,selectedGrid.GlobalData.Tot)
        Cinstance = matrixC.matrixC(integrationPoints,realElementInstance,selectedGrid.GlobalData.SpecificHeat, selectedGrid.GlobalData.Density)

        agreg.aggregate(Hinstance.H, HBCinstance.HBC, Pinstance.vectorP, Cinstance.C, nodesGlobalIDs)
          
    s = solve.Solve(agreg.HG, agreg.PG, agreg.CG, selectedGrid.GlobalData.InitialTemp,selectedGrid.GlobalData.SimulationStepTime, selectedGrid.GlobalData.SimulationTime)
    s.print_simulation()
    paraView = visualizeData.ParaView(selectedGrid,s.t_nonstationary_table)


def debug_output():
    with open('output.txt', 'w') as f:
            print("Uruchomiono program z opcja -v\n", file=f)
            print("\nDane globalne\n",file=f)
            selectedGrid.GlobalData.printGlobalData(f)
            print("\nWezly elementow siatki:\n", file=f)
            print(selectedGrid.Nodes.XY, file=f)
            print("\nIndeksy elementow siatki:\n", file=f)
            print(selectedGrid.Elements.ID, file=f)
            print("\nIndeksy punktow warunkow brzegowych:\n", file=f)
            print(selectedGrid.BC.BCs, file=f)
            print("\nTablica NodesBC (wartosc zmiennej warunku brzegowego dla kazdego indeksa node'a):\n", file=f)
            print(selectedGrid.Nodes.nodeBC, file=f)
            pc = gauss_integration.gaussIntegrationSchemaiD(integrationPoints,2)
            print("\nKolejnosc punktow calkowania:", file=f)
            j = 0
            for i in pc.P:
                print("(" + str(j) + ")" +"ksi: " + str(i[0]) + "     eta: " + str(i[1]) + "\n", file=f)
                j+=1
            print("\nWagi kolejnych punktow calkowania:", file=f)
            j = 0
            for i in pc.A:
                print("(" + str(j) + ")" +"A: " + str(i) + "\n", file=f)
                j+=1
            ue = universal_element.UniversalElement2D(integrationPoints)
            print("\n\nTablica dN / dEta dla wszystkich punktow calkowania\n", file=f)
            print(ue.dNdEtaTab, file=f)
            print("\n\nTablica dN / dKsi dla wszystkich punktow calkowania\n", file=f)
            print(ue.dNdKsiTab, file=f)
            print("\nZbior B[x,y] dla wszystkich punktow warunkow brzegowych (kolejnosc od skrajnie dolny-lewy, przecienie do wskazowek zegara (universal element))\n", file=f)
            print(ue.B, file=f)
            print("\nMacierze funkcji ksztaltu dla wszystkich punktow brzegowych zawierajace [N1,N2,N3,N4] (universal element)\n", file=f)
            print(ue.BCis_shapeFuncsTab, file=f)
            print("\nMacierze funkcji ksztaltu [N1,N2,N3,N4] (universal element)\n", file=f)
            print(ue.shapeFunTab, file=f)
                
                
            agreg = aggregation.Agregation(selectedGrid.GlobalData.NodesNumber)

            for i in range(selectedGrid.GlobalData.ElementsNumber):    
                print("\n(" + str(i) + ")" + "Wezly elementow siatki", file=f)
                nodesData = selectedGrid.getElement(i)
                print(nodesData[0],file=f)
                print("\n(" + str(i) + ")" + "Jakobiany oraz wyznaczniki jakobianow dla podanych wezlow dla kolejnych punktow calkowania", file=f)
                for k in range(integrationPoints**2):
                    j = universal_element.JacobyMatrix(nodesData[0],ue,k, integrationPoints)
                    print("\nPC[" + str(k) + "]:", file=f)
                    print(j.jacobian, file=f)
                    print("det(" + "PC[" + str(k) + "]" + "):", file=f)
                    print(np.linalg.det(j.jacobian), file=f)
                r= real_element.RealElement2D(integrationPoints,ue,nodesData[0])
                print("\n(" + str(i) + ")" + "Macierz dN / dX dla wskazanych wezlow", file=f)
                print(r.dNdXTab, file=f)
                print("\n(" + str(i) + ")" + "Macierz dN / dY dla wskazanych wezlow", file=f)
                print(r.dNdYTab, file=f)
                h = matrixH.matrixH(integrationPoints,r,selectedGrid.GlobalData.Conductivity)
                print("\n(" + str(i) + ")" + "Macierz H wskazanych wezlow", file=f)
                print(h.H, file=f)
                hbc = matrixHBC.matrixHBC(integrationPoints,ue,nodesData[0],nodesData[1],selectedGrid.GlobalData.Alfa)
                print("\n(" + str(i) + ")" + "Macierze HBC dla poszczegolnych scian elementu", file=f)
                print(hbc.HBCs, file=f)
                print("\n(" + str(i) + ")" + "Macierz HBC dla elementu", file=f)
                print(hbc.HBC, file=f)
                print("\n(" + str(i) + ")" + "Wektory P dla poszczegolnych scian elementu", file=f)
                P = vectorP.vectorP(integrationPoints,ue,nodesData[0],nodesData[1],selectedGrid.GlobalData.Alfa,selectedGrid.GlobalData.Tot)
                print(P.vectorsP, file=f)
                print("\n(" + str(i) + ")" + "Wektor P dla elementu", file=f)
                print(P.vectorP, file=f)
                c = matrixC.matrixC(integrationPoints,r ,selectedGrid.GlobalData.SpecificHeat, selectedGrid.GlobalData.Density)
                print("\n(" + str(i) + ")" + "Macierz C dla elementu", file=f)
                print(c.C, file=f)
                nodesGlobalIDs = selectedGrid.getElementIDs(i)
                agreg.aggregate(h.H, hbc.HBC, P.vectorP, c.C, nodesGlobalIDs)
            print("\n\nMacierz zagregowana H i HBC", file=f)
            agreg.print_HG(f)
            print("\n\nWektor P zagregowany", file=f)
            agreg.print_PG(f)
            print("\n\nMacierz zagregowana C", file=f)
            agreg.print_CG(f)
            s = solve.Solve(agreg.HG, agreg.PG, agreg.CG, selectedGrid.GlobalData.InitialTemp,selectedGrid.GlobalData.SimulationStepTime, selectedGrid.GlobalData.SimulationTime)
            print("\n\nWektor temperatury po zakonczeniu symulacji we wszystkich liczonych wezlach", file=f)
            s.print_tvector(s.t_nonstationary_end,f)
            print("\n\nTemperatury max i min w kazdym kroku symulacji (krok -> [s])", file=f)
            s.print_simulation(f)
            print("\n\nTabela obliczen symulacji dla wszystkich krokow i wszystkich wezlow", file=f)
            print(s.t_nonstationary_table, file=f)

    f.close()




if __name__ == "__main__":
    try:
        if args.verbose:
            debug_output()
        elif args.profile:
            profiler = cProfile.Profile()
            profiler.enable()
            main()
            profiler.disable()
            with open('efficiency.txt', 'w') as f:
                ps = pstats.Stats(profiler, stream=f)
                ps.sort_stats('tottime')
                ps.print_stats()
        else:
            main()
                        
    except (TypeError, ValueError, IndexError) as e:
        print(str(e))






