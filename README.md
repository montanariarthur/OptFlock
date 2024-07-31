# Optimal Flocking Dynamics
Codes for the real-time optimization of flocking dynamics with heterogeneous or homogeneous agents.

See the references below for more details.

# Usage

- `main_CMflock.m` : Compares the real-time optimization of the "centroid model" of flocking (Eq. 1 of Ref. 1) for heterogeneous and homogeneous systems. In this model, the agents are tasked to converge to a pre-specified formation while tracking a virtual target moving across space. The communication network is defined by an all-to-all network, weighted according to the relative distance between agents. Data packages are exchanged periodically, resulting in a piecewise-constant model.

- `main_TDconsensus.m` : Compares the stability optimization of the time-delay consensus model (Eq. 8 of Ref. 1) for heterogeneous and homogeneous systems. In this model, the agents must achieve consensus in position and velocity, converging to the same state. If the system is unstable, the flock of agents fragments into isolated groups.

- `main_OSflock.m` : Compares the real-time optimization of the "Olfati-Saber model" of flocking (Eq. 9 of Ref. 1) for heterogeneous and homogeneous systems. In this model, the agents are tasked to form a cohesive lattice structure while tracking a virtual target. The communication range is limited, resulting in a sparse interaction network. The lattice structure is emergent since the final configuration depends on the agents' initial conditions. The simulation can be performed for agents moving in free space or maneuvering around obstacles.

# Dependences

- `EigOptimization` : This folder contains optimization routines (`beta_opt.m`) used to minimize the Lyapunov exponent (calculated in `opteigreal.m`) associated with the corresponding flocking models (`CM`, `TD`, `OS`) for heterogeneous (`het`) and homogeneous (`hom`) systems.
  
- `FlockODEs` : This folder contains the ODEs describing the CM (`CMflock_piecewise.m`) and OS (`OSflock.m`) flocking models. The functions `agent_coord.m` compute the performance metrics (tracking error, lattice deviation energy, relative connectivty) and the XY coordinates of the agents/target for each time instant in the corresponding models.

- The  time-delay consensus model (`main_TDconsensus.m`) requires installation of the DDE-BIFTOOL toolbox. These codes were tested on version 3.0. See Ref. 3 and their website (https://sourceforge.net/projects/ddebiftool/) for more details.

# License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".

# References
1.  AN Montanari, AED Barioni, C Duan, AE Motter. Optimal flock formation induced by heterogeneity. (2024)
2.  R Olfati-Saber. Flocking for multi-agent dynamic systems: Algorithms and theory. *IEEE Transactions on Automatic Control*, **51**:401-420 (2006).
3.  K Engelborghs, T Luzyanina, D Roose. Numerical bifurcation analysis of delay differential equations using DDE-BIFTOOL. *ACM Transactions on Mathematical Software*, **28**:1-21 (2002).

