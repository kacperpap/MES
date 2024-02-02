# Application of Finite Element Method for heat transfer simulation

## License
For details, see [LICENSE.md](https://github.com/kacperpap/MES/edit/master/LICENSE.md)

## Finite Element Method (FEM) Concept

The core idea of the Finite Element Method (FEM) is discretization. This process involves transforming any continuous value into a discrete model based on a finite number of nodes.

The entire algorithm consists of several key steps:

1. **Node Definition**: The first part is defining the nodes of the grid.

2. **Parameterization**: Next, we parameterize values, such as temperature, at each node and treat them as parameters to be determined.

3. **Subdivision into Finite Elements**: We then divide the entire grid into subareas, or finite elements, which share common nodes and approximate the shape of the area.

4. **Interpolation within Each Element**: Within each element, we interpolate values (e.g., temperature) using a polynomial based on its nodal values.

5. **Functional Minimization**: Finally, we minimize a functional by selecting appropriate nodal values, which are determined by solving a system of algebraic equations. For heat conduction, we minimize a functional that corresponds to Fourier's equation with imposed boundary conditions.

This method provides a powerful approach for numerically solving a variety of problems that can be represented by partial differential equations, such as heat conduction, fluid dynamics, and structural analysis.

## Detailed Description

A more comprehensive description is available in the report located in the repository. To familiarize yourself with a more detailed theoretical explanation of the method and a description of the code, please refer to this [report](https://github.com/kacperpap/MES/edit/master/MES_sprawozdanie_Papuga_Kacper.pdf)

## Usage

1. **Clone the Repository**

2. **Install Dependencies**: All dependencies are listed in the `requirements.txt` file. To install them, open the main folder in the command prompt and enter the following command: `python.exe -m pip install -r requirements.txt`.

3. **Set Parameters**: Parameters are manually set in the `main.py` file. You need to specify the `grid` and the number of `integrationPoints` you want to start the simulation with. Available grids are commented in the code. If you want to add your own mesh (the mesh must be composed of quadrangles), familiarize yourself with the structure of the added one. Implemented Gauss integration schemes include 2, 3, and 4 integration point schemes. To set these parameters, modify the code in the `main.py` file as shown below:
   ```python
      selectedGrid = grid.Grid("Test2_4_4_MixGrid.txt")
      integrationPoints = 2

   
4. **Run code**: Run your code by entering the command `python.exe ./main.py` or by double-clicking on the `main.py` file.

5. **Run with additional options**: You can run your code with additional flags: `python.exe ./main.py --verbose` / `python.exe ./main.py -v` or `python.exe ./main.py --profile` / `python.exe ./main.py -p`. The verbose option generares `output.txt` file which is useful for debugging purposes. It shows step-by-step calculations. Be aware that running that on a large mesh will take some time. The profile option uses `cProfile` and `pstats` to generate an ```efficiency.txt``` file, which shows the efficiency stats of computing the simulation.


## Example of terminal simulation results:
<p align="center">
  <img src="https://i.imgur.com/sAM7Luz.png" width="350" height="500" />
</p>


## ParaView simulation

The code generates `.vtk` files, which can be uploaded into the ParaView. Example of simulation in ParaView shown below:
![paraview_simulation](https://i.imgur.com/0TaG2A5.gifv)


## Sources

The primary substantive basis of the code presented in this repository is derived from exercises and lectures. Additionally, content is sourced from the following book:

Milenin, A. (2010). *Podstawy metody elementów skończonych: Zadania termomechaniczne*. AGH University of Science and Technology Press. Available at https://wydagh.agh.edu.pl/produkt/855-podstawy-metody-elementow-skonczonych



