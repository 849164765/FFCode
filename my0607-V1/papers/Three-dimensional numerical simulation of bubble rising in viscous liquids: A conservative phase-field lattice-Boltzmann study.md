# Three-dimensional numerical simulation of bubble rising in viscous liquids: A conservative phase-field lattice-Boltzmann study

Cite as: Phys. Fluids 31, 063106 (2019); https://doi.org/10.1063/1.5096390 

Submitted: 17 March 2019 . Accepted: 01 June 2019 . Published Online: 21 June 2019 

Ang Zhang , Zhipeng Guo, Qigui Wang, and Shoumei Xiong 

## COLLECTIONS

F This paper was selected as Featured 

![](images/57a0f3a8badbe8112dd4dfe5a6a7305a3bb9bbf224931368ec5b35c76f16ae05.jpg)



View Online


![](images/d276f244f732d3307b530f12c3ab48783d16bbe436853acaf8c7262413b5d190.jpg)


![](images/dd43ba05315a3528c412fb73923cec226632ff0a07928348f46b3d2cad6aff37.jpg)



CrossMark


## ARTICLES YOU MAY BE INTERESTED IN

Fluid dynamics of oscillatory flow in three-dimensional branching networks 

Physics of Fluids 31, 063601 (2019); https://doi.org/10.1063/1.5093724 

On the concept and theory of induced drag for viscous and incompressible steady flow 

Physics of Fluids 31, 065106 (2019); https://doi.org/10.1063/1.5090165 

Large-eddy simulation of an open-channel flow bounded by a semi-dense rigid filamentous canopy: Scaling and flow structure 

Physics of Fluids 31, 065108 (2019); https://doi.org/10.1063/1.5095770 

![](images/2503e2cb34d0cd1f848ca8566e1b7684ca418a9248946a2c954ae4c0335377a9.jpg)


# Three-dimensional numerical simulation of bubble rising in viscous liquids: A conservative phase-field lattice-Boltzmann study

Cite as: Phys. Fluids 31, 063106 (2019); doi: 10.1063/1.5096390 

Submitted: 17 March 2019 • Accepted: 1 June 2019 • 

Published Online: 21 June 2019 

![](images/ed66cb5abb4346d077eae7f1ad93cccfcbe05633ec92114ae4311d423d84e3b6.jpg)


![](images/e0460c8fbb142d0e009b2769ae58be400a14fc0ce37767f9219cda0e726beae9.jpg)


![](images/83e1934a3e86832052b41f3d4639d6ce1d94e0d5773e005f9afce9add79727e2.jpg)


Ang Zhang,1 Zhipeng Guo,1,a) Qigui Wang,2 and Shoumei Xiong1,3,a) 

## AFFILIATIONS

1 School of Materials Science and Engineering, Tsinghua University, Beijing 100084, China 

2Materials Technology, GM Global Propulsion Systems, Pontiac, Michigan 48340-2920, USA 

3Key Laboratory for Advanced Materials Processing Technology, Ministry of Education, Tsinghua University, Beijing 100084, China 

a)Authors to whom correspondence should be addressed: zhipeng_guo@mail.tsinghua.edu.cn and smxiong@tsinghua.edu.cn 

## ABSTRACT

Simulating bubble rising in viscous liquids is challenging because of the large liquid-to-gas density ratio and complex topological evolution of the gas-liquid interface. In this study, a conservative phase-field model is employed to accurately track the interface during bubble rising, and the lattice Boltzmann model is used to determine the flow field driven by the buoyancy force and the surface tension force. To facilitate large-scale three-dimensional simulations, a parallel-adaptive mesh refinement algorithm is developed to reduce the computing overhead. The simulated bubble shapes under different configurations are compared with the shape chart through experiments [D. Bhaga and M. E. Weber, “Bubbles in viscous liquids: shapes, wakes, and velocities,” J. Fluid Mech. 105, 61–85 (1981)]. The influence of the numerical parameters (including domain size, surface tension, liquid viscosity, gravity, and density ratio) on the bubble dynamics is investigated, which demonstrates the capability of the current numerical scheme in simulating multiphase flow. Furthermore, complex topology changes including the bubble coalescence, splitting, and interplay with obstacles (i.e., squeeze deformation and bubble splitting) are simulated and compared in different cases, i.e., with different Reynolds, Eötvös, and Morton numbers. The effect of the initial bubble spacing on the coalescence of the two bubbles and the influence of boundary conditions on multiple bubble dynamics are investigated. When the bubbles can be completely blocked by the obstacle is quantified in terms of the obstacle width. Numerical results validate the robustness of the present numerical scheme in simulating multiphase flow. 

Published under license by AIP Publishing. https://doi.org/10.1063/1.5096390 

## I. INTRODUCTION

As a typical case of multiphase flows, bubble rising in quiescent viscous liquids has aroused increasing interest in engineering applications, geological sciences, and fundamental physics, such as sheet cavitation on propellers,1,2 porosity defect in casting,3 cavitation bubble during solidification,4 explosive volcanic eruptions,5 and vaporization and boiling.6 The complex interaction among buoyancy, surface tension, geometric evolution, and momentum inertia causes difficulties in characterizing the nonlinear behavior of rising bubbles. Further research on bubble dynamics can contribute to a better understanding of the multiphase flow mechanism, which helps designing two-phase flow devices in multiphase reactor systems.7 

Extensive work has been performed experimentally,8–11 theoretically,12–14 and computationally.7,15–17 correlated the bubble shape and terminal rise velocity in viscous liquids. Maxworthy et al.9 uncovered scaling regimes of rising bubbles and showed the importance of different dynamical interactions. Wang et al.10 investigated the interaction mechanism of multiple synchronized bubbles, especially on the coalescence of bubbles and the motion of the joined bubbles. Sharaf et al.11 studied the shapes and paths of rising bubbles in quiescent liquids and obtained a phase plot showing the distinct behavior of rising bubbles in a large range of Galilei and Eötvös numbers. Many experimental studies, however, are limited by the precise control of the flow field and accurate monitoring of the evolving interface.7 It is difficult to obtain the detailed velocity distribution surrounding and within the bubble in the experiment, and generalizing the gas-liquid regimes to other two-phase systems might be hazardous or imprecise.18 On the other hand, recent advances in computing power and numerical methods facilitate the investigation of the bubble dynamics, the core of which can be summarized as three parts: interface tracking, flow solver, and discrete mesh.7 

There are two basic categories to track the interface including the sharp-interface model and the diffuse-interface model. Different from the infinitesimally thin boundary with jump conditions in the sharp-interface model, the diffuse-interface model assumes a finitewidth interface, across which the fluid properties vary smoothly based on thermodynamics.19–22 Numerically speaking, the diffuseinterface model is preferred over the other because no extra treatment is required to track the interface. A typical technique of the diffuse-interface model is the phase-field method (PFM), in which the coexisting phases are distinguished by an introduced order parameter, and irregular topological changes can be easily handled by solving phase-field equations.23–25 The interface tracking in numerous multiphase-flow scenarios, such as bubble rising, milk crown (i.e., droplet falling and splashing),26 and Rayleigh-Taylor instability,27 has been achieved using the phase-field concept by Fakhari et al.,28,29 Chiu and Lin,26 Zu and He,27 and Shi et al.30 

The PFM can be categorized into nonconservative and conservative models. Sun and Beckermann31 presented a general interface tracking method based on the phase-field theory. Chiu and Lin26 extended this method to a conservative form, which could preserve the total mass like the Cahn-Hilliard equation.32 Takada et al. proposed a revised Allen-Cahn advection equation in a conservative form based on the work performed by Chiu and Lin.26 Geier et al.34 demonstrated the mass-conserving characteristic of the phase-field version reformulated by Chiu and Lin2 based on the flux density theory. The present conservative PFM is actually a modified Allen-Cahn approach with a second-order conservation property.34,35 Achieving mass conservation by avoiding numerical dissipation is the key to precisely modeling multiphase flow, which has also been explored in other numerical methods such as the modified level-set method proposed by Shu et al.36 and the mass-conserving multiphase lattice Boltzmann model developed by Li et al.1 

The flow solver can be divided into a classical Navier-Stokes (NS) solver and a kinetic-based lattice Boltzmann model (LBM). The LBM can be restored to the former by the Chapman-Enksog expansion37 but has more attractive advantages for the case with high solid fraction and geometry complexity.20 In addition, the solution of the Poisson equation for the pressure in the NS solver is avoided in the LBM, which simplifies the algorithm structure and becomes easier for parallel computing.38 Unlike the momentum conservation on the continuum level in the NS solver, the LBM solves the macroscopic flow by repeated collision-propagation dynamics of assumed discrete particles. This makes the LBM a mesoscopic approach, and thus, it can simulate problems with more sophisticated physics due to essentially mesoscopic interphase boundaries. 7,29 

Simulation of the multiphase flow is a challenging task because the diffuse interface must be characterized using sufficient fine grids to obtain high accuracy. For the bubble rising, the simulated bubble size must be much larger than the interface width, and the minimum grid number across the diameter needs to increase with the deformation degree.7 Moreover, the domain size also needs to be significantly larger than the bubble size to avoid the wall effects. Accordingly, using uniform grids to discrete the multiphase systems is unnecessary and time-consuming due to huge computing consumption and memory usage. It is thus highly anticipated to develop new robust algorithms which can meet the requirement in both mesh resolution and computation cost. 

The combination of parallel computing and adaptive mesh refinement (AMR) makes efficiently solving the equations tractable. Grids with high resolution cover the phase interface, while coarser grids are employed elsewhere. The whole grid structure is adjusted dynamically as the interface evolves. The AMR strategy has been integrated into the multiphase-flow LBMs by Töelke et al.,39 Yu and Fan,7 Janssen and Krafczyk,40 Fakhari et al.,29 and Wang et al.41 The development of the massively parallel computing algorithm enlarges the simulation scale, but further employment of AMR to simulate a large-scale multiphase flow is very limited. 

As a result, most simulations can only handle small bubbles with spherical or ellipsoidal shape, and the number of bubbles is limited to two to three when exploring the interbubble interaction.42 Zhang et al.43 modeled the coalescence of two bubbles using the boundary integral method. Chen et al.44 studied the formation and interaction of two bubbles growing from the orifices using the freeenergy-based LBM. To maintain numerical stability, current studies are still struggling because of the large liquid-to-gas density ratio and complex topological changes during bubble coalescence and breakup. Besides, because of the sharply evolving topology, investigation on the dynamical behavior between bubbles and obstacles is very limited except for those reported by Alizadeh et al.18 Recently, Qian et al.45 investigated the interaction between the rising bubble and a solid wall using an arbitrary Lagrangian-Eulerian method, but the bubble deformation is limited compared with those where the bubble breaks up and/or wraps the obstacle. 

In this work, a conservative PFM, which is derived from Fick’s second law and Cahn-Hilliard theory,46 is employed to simulate the evolution of the phase interface. The conservative interfacetracking equation is essentially equivalent to those proposed by Chiu and Lin26 and Geier et al.,34 and the mass-conserving feature has been thoroughly validated by reproducing highly complex topology changes under different steady flows in our recent work.46 The velocity field determined using the LBM is coupled into the PFM to move and deform the phase interface, in which the driving force comprises the buoyancy force induced by the density difference and the surface tension force. To achieve large-scale three-dimensional (3D) simulations with multiple bubbles, a parallel-adaptive mesh refinement algorithm (Para-AMR), which was developed by the current authors to efficiently simulate dendritic and eutectic evolution,47–50 is combined with the present phase-field lattice-Boltzmann approach. The objective of this work is to track complex topology changes during bubble rising, which include bubble coalescence, splitting, and interplay with obstacles. By combining the conservative phase-field lattice-Boltzmann approach and parallel-adaptive mesh refinement algorithm, large-scale 3D simulations become tractable, which lays the foundation for further study on the multiphase flow. 

This work is organized as follows. In Sec. II, the mathematical methods, including the conservative PFM, the LBM, and the parallel-adaptive mesh refinement algorithm, are introduced, followed by the numerical validation and discussion on parameter influence in Sec. III. Typical numerical results are presented in Sec. IV to illustrate multiple bubble dynamics and the interplay with the obstacle in 3D cases. Finally, the concluding remarks are summarized in Sec. V. 

## II. MATHEMATICAL METHODS

## A. Conservative phase-field model

A conservative phase-field model (PFM), which is derived from Fick’s second law and Cahn-Hilliard theory (see the $\mathrm { A p p e n d i x } ) , ^ { 4 6 }$ is employed to characterize the shape evolution during bubble rising. The phase field variable φ changes smoothly from 0 in the gas to 1 in the liquid, which acts as a counterpart of the level-set function in the level-set method51 or the volume fraction in the volume of fluid method,52 

$$
\partial_ {t} \phi + \nabla \cdot (\boldsymbol {v} \phi) = \nabla \cdot \left[ M \left(\nabla \phi - \frac {\phi - \phi^ {2}}{\delta} \frac {\nabla \phi}{| \nabla \phi |}\right) \right], \tag {1}
$$

where t stands for time, M is the mobility, δ measures the width of the diffuse interface, and v is the macroscopic velocity vector. The interface location can be extracted by setting $\phi = 0 . 5 .$ . As a phase-field version of the modified Allen-Cahn equation, Eq. (1) is essentially an extension of the general interface-tracking technique proposed by Sun and Beckerman,31 which is conservative and driven by curvature. At equilibrium, the change in φ normal to the interface can be described using a hyperbolic tangent profile,53 

$$
\phi^ {e q} = \frac {1}{2} + \frac {1}{2} \tanh \left(\frac {\varepsilon}{2 \delta}\right), \tag {2}
$$

where ε is the coordinate normal to the interface. 

Equation (1) is discretized using the central difference method in space and the fourth-order Runge-Kutta method in time. For the divergence operator, a net flux control volume method is employed;47 that is, the phase field variable is fixed at the cell center, and the divergence is calculated by summing the fluxes at all walls of the cell. 

Assuming that the gas density variation caused by the liquid hydraulic pressure is neglected under isothermal conditions, the bubble rising in quiescent liquids can be regarded as a multifluid system with two incompressible and immiscible Newtonian fluids.15 But in the PFM, only one single set of governing equations is solved for the two different fluids, and the fluid properties vary across the diffuse interface, i.e., as a function of the phase field variable φ. Taking the macroscopic density ρ for instance, it can be calculated using a linear interpolation between the bulk phases, i.e., $\rho = \rho _ { l } \phi + \rho _ { g }$ (1 $- \ \phi )$ , where $\rho _ { l }$ and $\rho _ { g }$ are the density of the liquid and gas phases, respectively. With this “one-fluid approach,” extra treatment on the jump condition and discontinuity across the interface is avoided, which decreases computational complexity and can handle complex topological changes with high stability. 

## B. Lattice Boltzmann model

Bubble rising is an intricate behavior that involves interface separation and coalescence. Compared with the conventional Navier-Stokes solver, the kinetic-based LBM is a mesoscopic approach and has appealing advantages in simulating multiphase flow because the phase interface is mesoscopic in nature.29 The consistency between the LBM and the Navier-Stokes solver can be revealed by the Chapman-Enskog analysis.37 The distribution function with the discrete force term $\dot { F } _ { i } ( r , t )$ is expressed $\mathsf { a s } ^ { 7 , 5 4 }$ 7,54 

$$
f _ {i} (\boldsymbol {r} + \delta \boldsymbol {r}, t + \delta t) = f _ {i} (\boldsymbol {r}, t) - \left(f _ {i} (\boldsymbol {r}, t) - f _ {i} ^ {e q} (\boldsymbol {r}, t)\right) / \tau + F _ {i} (\boldsymbol {r}, t) \delta t, \tag {3a}
$$

with 

$$
F _ {i} = \left(1 - \frac {1}{2 \tau}\right) w _ {i} \left(3 \frac {\boldsymbol {c} _ {i} - \boldsymbol {v}}{c ^ {2}} + \frac {9 (\boldsymbol {c} _ {i} \cdot \boldsymbol {v}) \boldsymbol {c} _ {i}}{c ^ {4}}\right) \cdot \boldsymbol {F}, \tag {3b}
$$

where $f _ { i } ( \pmb { r } , t )$ is the particle distribution function and represents the number of particles along the ith direction at the position r and time $t , \tau = 0 . 5 + 3 \nu / ( c ^ { 2 } \delta t )$ is the relaxation factor, ν is the liquid kinetic viscosity, F is the force vector, $c = 1$ is the lattice velocity, and ci and wi are the discrete velocity and corresponding weight coefficient along the ith direction, respectively, both of which are dependent on the employed discrete velocity model. The 2D nine-velocity (D2Q9) and 3D nineteen-velocity (D3Q19) models37 are employed in this work, and taking the D3Q19 model for instance, $w _ { 0 } = 1 / 3$ , w1–6 = 1/18, and $w _ { 7 - 1 8 } = 1 / 3 6 ,$ and three components of the discrete velocity $\begin{array} { r } { \pmb { c } _ { i } = ( u , \nu , } \end{array}$ w) are listed in Table I. 

The equilibrium distribution function $f _ { i } ^ { e q } ( \pmb { r } , t )$ is given by 

$$
f _ {i} ^ {e q} = \rho w _ {i} \left(1 + \frac {3 \boldsymbol {c} _ {i} \cdot \boldsymbol {v}}{c ^ {2}} + \frac {9 (\boldsymbol {c} _ {i} \cdot \boldsymbol {v}) ^ {2}}{2 c ^ {4}} - \frac {3 \boldsymbol {v} \cdot \boldsymbol {v}}{2 c ^ {2}}\right), \tag {4}
$$

where $\begin{array} { r } { \pmb { \nu } = \sum _ { i } f _ { i } \pmb { c } _ { i } / \rho + \delta t \pmb { F } / ( 2 \rho ) } \end{array}$ is the flow velocity in the presence of the force, which guarantees a second-order space-time accuracy.37 

To characterize the force interaction between the rising bubble and the viscous liquid, the force term is expressed as $\dot { F } = F _ { b }$ $+ \ F _ { s } ,$ , where $F _ { b } = - \ ( \rho - \rho _ { l } ) \pmb { g }$ is the buoyancy force and $F _ { s } = \mu _ { \phi } \nabla \phi$ is the surface tension $\mathrm { f o r c e } , ^ { 2 9 }$ in which g is the gravity acceleration and $\mu _ { \phi } = 4 \beta \phi ( \phi - 1 ) ( \phi - 1 / 2 ) - \kappa \nabla ^ { 2 } \phi$ is the chemical potential, where $\beta$ and κ are the constants determined by the surface tension σ and interface width δ, i.e., $\beta = 3 \sigma / \delta$ and $\kappa = 6 \sigma \delta$ . The derivative of the chemical potential $\mu _ { \phi }$ is canceled in calculating the surface tension force, which diminishes the numerical stiffness and enhances the locality of the collision process.29 It is noted that the buoyancy force induced by the density difference works as the driving force for the bubble rising, and the surface tension force which is dependent on the phase gradient accounts for the interface deformation. 


TABLE I. Three velocity components of ci = (u, v, w) in the D3Q19 model.


<table><tr><td>i</td><td>0</td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td><td>13</td><td>14</td><td>15</td><td>16</td><td>17</td><td>18</td></tr><tr><td>u</td><td>0</td><td>1</td><td>-1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>-1</td><td>-1</td><td>1</td><td>1</td><td>1</td><td>-1</td><td>-1</td><td>0</td><td>0</td><td>0</td><td>0</td></tr><tr><td>v</td><td>0</td><td>0</td><td>0</td><td>1</td><td>-1</td><td>0</td><td>0</td><td>1</td><td>1</td><td>-1</td><td>-1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>1</td><td>-1</td><td>-1</td></tr><tr><td>w</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>-1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>-1</td><td>-1</td><td>1</td><td>1</td><td>-1</td><td>-1</td><td>1</td></tr></table>

To maintain high stability, extra attention should be paid to the discretization scheme of the divergence and Laplace operators in the force formulas.55 Taking 3D for instance, the two derivatives are evaluated using an eighteen-point difference scheme, i.e., 

$$
\nabla \phi = \frac {3}{2 c \delta t} \left\{ \begin{array}{l} \boldsymbol {i} \sum_ {\alpha = 1} ^ {1 8} w _ {\alpha} \boldsymbol {c} _ {\alpha} \cdot \boldsymbol {i} [ \phi (\boldsymbol {r} + \boldsymbol {c} _ {\alpha} \delta t) - \phi (\boldsymbol {r} - \boldsymbol {c} _ {\alpha} \delta t) ] \\ \boldsymbol {j} \sum_ {\alpha = 1} ^ {1 8} w _ {\alpha} \boldsymbol {c} _ {\alpha} \cdot \boldsymbol {j} [ \phi (\boldsymbol {r} + \boldsymbol {c} _ {\alpha} \delta t) - \phi (\boldsymbol {r} - \boldsymbol {c} _ {\alpha} \delta t) ], \\ \boldsymbol {k} \sum_ {\alpha = 1} ^ {1 8} w _ {\alpha} \boldsymbol {c} _ {\alpha} \cdot \boldsymbol {k} [ \phi (\boldsymbol {r} + \boldsymbol {c} _ {\alpha} \delta t) - \phi (\boldsymbol {r} - \boldsymbol {c} _ {\alpha} \delta t) ] \end{array} \right. \tag {5a}
$$

$$
\nabla^ {2} \phi = \frac {1}{6 (c \delta t) ^ {2}} \left\{\sum_ {\alpha = 1} ^ {1 8} \phi (\boldsymbol {r} + \boldsymbol {c} _ {\alpha} \delta t) + \sum_ {\alpha = 1} ^ {6} \phi (\boldsymbol {r} + \boldsymbol {c} _ {\alpha} \delta t) - 2 4 \phi (\boldsymbol {r}) \right\}, \tag {5b}
$$

where $i , j ,$ and k are the unit vectors along the x, y, and z axis, respectively. This scheme computes 18 nearest neighbors on the D3Q19 lattice to avoid sharp gradient that might cause instability and inaccuracy, which improves the numerical robustness in simulating bubble flow with large liquid-to-gas density ratio.55 

It is noted that both the spatial and time steps (δx and δt) in the LBM are rescaled to 1, and a dimensional transformation is required from the LBM to the PFM.25,38 In particular, the velocity calculated using the LBM needs to multiply dxmax/dt before being incorporated into Eq. (1), where $d x _ { m a x }$ and dt are the maximum mesh size (see the multilevel data structure in Sec. II C) and the time step during simulation. 

## C. Para-AMR algorithm

To characterize the diffuse interface, sufficiently small meshes must be generated to maintain a high resolution. Using AMR can not only satisfy the high accuracy but also reduce the computing overhead. The fine mesh resolution is adopted near the bubble interface, and the mesh size increases with the distance away from the interface. The mesh structure updates the layout dynamically according to the interface motion. 

A predefined gradient criterion is employed to trigger refinement, i.e., 

$$
\left| \nabla \phi \right| \geq \xi , \tag {6}
$$

where ξ is a threshold determined using numerical tests. Once Eq. (6) is reached, the local meshes are tagged and then subdivided by bisecting the parent cells. The refinement process is continued until the required maximum refinement level is reached. The tagged cells are clustered into so-called patch boxes which comprise a collection of meshes with the same size (depicted using the same color in Fig. 1). The data between neighboring levels are communicated when performing interpolation (from coarse to fine) and restriction (from fine to coarse) operations.47,50,56 Figure 1 shows the hierarchical grid structure during bubble rising, in which four grid levels (from level I to level IV) are used to decompose the whole domain. The gas-liquid interface is perfectly characterized by the minimum sized meshes, while the bulk phases are covered by coarser meshes. 

![](images/308a921b6f4f13cd36444b522fa3ee7dcf023941afcb73d20172e842901e217b.jpg)


![](images/189d77c4f36eba9499c777ad170f2e118e7574e091b0f89e68a4aef5ef99f4b1.jpg)


![](images/c5ee4fa3aeaedf941c2b784121761480a3ea251d35bcb39f4a8ecad6cc717033.jpg)


![](images/e3bc20c1cc1c43cedbffac03fb37589c2543a5a893deb2f2d3b4a9c1b56eac34.jpg)


![](images/0b3f1a6c9f217f07f601669fcbb349b3b1cc8266a3949c5f529928fddf64e272.jpg)



FIG. 1. (a) Grid layout for simulating bubble rising. $( \mathsf { b } _ { 1 } ) \mathsf { - } ( \mathsf { b } _ { 4 } )$ Each grid level in the hierarchical four-level grid structure. The mesh size increases with the distance away from the phase interface. Each color denotes a collection of meshes with the same size in the certain grid level, and a darker color corresponds to a finer mesh.


After constructing the multilevel data structure, all computing data are dispatched to multiple processes to achieve parallel computing based on Message Passing Interface (MPI).47 Then, each process has its own but different arrays of patch boxes. To update the data simultaneously, a layer of ghost cells is added at the external boundaries of each patch box to collect data. Based on the parallel and adaptive mesh refinement algorithm (Para-AMR), the 3D numerical simulations become tractable to implement. 

Table II lists the elapsed time during simulation of 3D rising bubbles for different conditions. The simulation parameters are set to be the same as case B in Table III, and all simulations last for a total of 10 000 time steps. The domain size is 160 × 384 × 160 in lattice units, i.e., nearly ten million meshes if a uniform grid is employed. When 24 parallel processes are dispatched, the elapsed time is shortened by about 16 times, i.e., from 46 h to less than 3 h. The employment of the AMR further speeds up the computation by reducing the computing data and improves the efficiency by 3 times. 


TABLE II. Elapsed time of the 3D rising bubble under different conditions.


<table><tr><td>Conditiona</td><td><eq>N_{p} = 1, N_{amr} = 1</eq></td><td><eq>N_{p} = 24, N_{amr} = 1</eq></td><td><eq>N_{p} = 24, N_{amr} = 4</eq></td></tr><tr><td>Elapsed time (s)</td><td>164 484</td><td>10 203</td><td>3472</td></tr></table>


$^ { \mathrm { a } } N _ { p }$ and $N _ { a m r }$ are the number of parallel processes and the number of grid levels. 



TABLE III. Simulation conditions of single bubble flow.


<table><tr><td>Case</td><td>Variable</td><td>Re</td><td>Eo</td><td>Mo</td><td>Shape<eq>^{8}</eq></td></tr><tr><td>A</td><td></td><td></td><td>0.327 35</td><td><eq>7.8047 \times 10^{-8}</eq></td><td>Spherical (s)</td></tr><tr><td>B</td><td></td><td></td><td>3.273 5</td><td><eq>7.8047 \times 10^{-5}</eq></td><td>Oblate ellipsoid (oe)</td></tr><tr><td>C</td><td>σ</td><td>25.892</td><td>65.47</td><td>0.624 38</td><td>Oblate ellipsoidal cap (oec)</td></tr><tr><td>D</td><td></td><td></td><td>654.7</td><td>624.38</td><td>Skirted with smooth skirt (sks)</td></tr><tr><td>E</td><td></td><td></td><td>935.29</td><td>1820.3</td><td>Skirted with smooth skirt (sks)</td></tr><tr><td>F</td><td></td><td>0.129 46</td><td></td><td><eq>1.2488 \times 10^{11}</eq></td><td>Spherical (s)</td></tr><tr><td>G</td><td></td><td>2.589 2</td><td></td><td><eq>7.8047 \times 10^{5}</eq></td><td>Oblate ellipsoidal cap (oec)</td></tr><tr><td>H</td><td>ν</td><td>25.892</td><td>327.35</td><td>78.047</td><td>Oblate ellipsoidal cap (oec)</td></tr><tr><td>I</td><td></td><td>89.284</td><td></td><td>0.552 01</td><td>Spherical cap with closed wake (scc)</td></tr><tr><td>J</td><td></td><td>647.31</td><td></td><td><eq>1.998 \times 10^{-4}</eq></td><td>Spherical cap with open wake (sco)</td></tr><tr><td>K</td><td></td><td>8.187 9</td><td>0.327 35</td><td><eq>7.8047 \times 10^{-6}</eq></td><td>Spherical (s)</td></tr><tr><td>L</td><td></td><td>36.617</td><td>6.547</td><td><eq>1.5609 \times 10^{-4}</eq></td><td>Oblate ellipsoid (oe)</td></tr><tr><td>M</td><td>g</td><td>81.879</td><td>32.735</td><td><eq>7.8047 \times 10^{-4}</eq></td><td>Oblate ellipsoidal disk (oed)</td></tr><tr><td>N</td><td></td><td>258.92</td><td>327.35</td><td><eq>7.8047 \times 10^{-3}</eq></td><td>Spherical cap with open wake (sco)</td></tr><tr><td>O</td><td></td><td>366.17</td><td>654.7</td><td><eq>1.5609 \times 10^{-2}</eq></td><td>Spherical cap with open wake (sco)</td></tr></table>

## III. MODEL EVALUATION

The conservative PFM was validated in our recent work by simulating interface evolution in different benchmark problems,46 in which the predefined velocity only depended on the physical location. For the bubble rising, the two-phase flow is induced through the difference of buoyancy force and surface tension force, both of which evolve over time and space. To validate the numerical accuracy, a series of numerical tests in both 2D and 3D cases are performed to evaluate the dependence of single bubble rising on the domain size and force parameters. Comparison with experiments is presented to further demonstrate the robustness of the numerical scheme. Unless stated otherwise, the default values are set in the lattice unit as $\rho _ { l } = 1 . 0 , \rho _ { g } = 1 . 0 \times 1 0 ^ { - 3 } , \nu = 1 . 8 7 5 \times 1 0 ^ { - 2 } , \sigma = 2 . 2 5 \times 1 0 ^ { - 3 }$ , $d _ { 0 } = 3 2$ , and δ = 2 in 2D and $\rho _ { l } = 1 . 0 , \rho _ { g } = 1 . 0 \times 1 0 ^ { - 3 } , \nu = 1 . 2 5 \times 1 0 ^ { - 2 } ,$ $\sigma = 1 \times 1 0 ^ { - 3 } , d _ { 0 } = 3 2$ , and δ = 2 in 3D, where d0 is the diameter of the initial circle/spherical bubble. A zero-Neumann boundary condition is set at all sides for the phase field, and a no-slip boundary condition is applied at all sides for the velocity. The bubble contour is extracted at $\phi = 0 . 5 ,$ , and the minimum mesh size $d x _ { m i n } \mathrm { i } s 0 . 8$ . 

The bubble rising in viscous liquids is a typical two-phase problem, and the complex nonlinear relation can be correlated by introducing three dimensionless numbers including the gravity Reynolds number (Re), Eötvös (Bond) number (Eo or Bo), and Morton number $\left( M o \right) , ^ { 2 9 }$ 

$$
R e = \frac {\sqrt {g \rho_ {l} \left(\rho_ {l} - \rho_ {g}\right) d _ {0} {} ^ {3}}}{\nu}, \quad E o = \frac {g \left(\rho_ {l} - \rho_ {g}\right) d _ {0} {} ^ {2}}{\sigma}, \quad M o = \frac {g \left(\rho_ {l} - \rho_ {g}\right) v ^ {4}}{\sigma^ {3} \rho_ {l} ^ {2}}. \tag {7}
$$

Accordingly, $R e \ : = \ : 2 5 . 8 9 ,$ Eo = 3.27, and $M o \ = 7 . 8 0 \times 1 0 ^ { - 5 }$ , and the bubble shape is expected to be oblate ellipsoid,8 the shape of which can be easily characterized by the length ratio of semimajor to semiminor axes (i.e., a/b, see the definition in the inset in Fig. 2). 

## A. Domain size dependence

During bubble rising in a confined domain, the distribution of flow field is affected by domain size, or the “wall effect” proposed by Krishna et $a l . ^ { 5 7 }$ Choosing a proper domain size, especially the transverse domain size, can both reduce the boundary effect and moderate the computing cost. Figure 2 shows the terminal velocity U and the length ratio of semimajor to semiminor axes vs the transverse domain size in both 2D and 3D cases, in which the insets show the bubble shape under different configurations. S is the domain size along both the x and z directions (note that the y direction is the rise direction of the bubble). As S/d0 increases, the terminal velocity first increases and then becomes stable, so does the length ratio, indicating the bubble shape becomes independent of the domain size for larger S/d0. 

![](images/1ed6a0be8ad1b0375f34ca8210c664414905e3a65dc340b9b6f25a2a348ff211.jpg)


![](images/0602d134f08a79772c096bb1e21f1c952296eeb7ae83eda2d61b720bad92be17.jpg)



FIG. 2. Terminal velocity and length ratio of semimajor to semiminor axes vs the transverse domain size in (a) 2D and (b) 3D cases, where S is the domain size along both x and z directions, $d _ { 0 }$ is the diameter of the initial circle/spherical bubble, and the insets show the bubble shape and the definition of semimajor and semiminor axes.


The viscous effect induced by the side walls constrains the transverse deformation and slows down the rise velocity of bubbles, which agrees with the experimental investigations reported by Krishna et al.57 In their experiments, the bubble shape and rise velocity in columns with different inside diameters were measured to investigate the wall effects. As the domain size increases, the viscous effect induced by the domain wall becomes weak and even completely disappear, which is essential to study the bubble dynamics.58 Figure 2 shows that both the bubble deformation and the terminal velocity become less affected by the side walls if $S / d _ { 0 } \ge 8 $ i n 2D or $S / d _ { 0 } \geq 5$ in 3D. When the domain size is less than 7d0 in 2D or 4d0 in 3D, the terminal velocity decreases, which reduces both transverse deformation and vertical movement. Accordingly, the domain size along both the x and z directions is fixed at 8d0 in 2D and 5d0 in 3D in the following simulations to optimize computing cost and numerical accuracy. In this case, the bubble dynamics only exhibits a deviation of ∼3.3% compared with that in the larger domain. For the vertical size along the rise direction (i.e., the y direction), the domain height is fixed at 12d0,15 which is large enough to observe the steady rise because the effect of the domain height can be negligible if the bubble does not approach the top wall. 

## B. Force dependence

Based on the experiment, Bhaga and Weber8 found that the bubble shape rising in viscous liquids can be divided into eight categories including spherical (s), oblate ellipsoid (oe), oblate ellipsoidal disk (oed), oblate ellipsoidal cap (oec), spherical cap with open or closed wake (sco or scc), and skirted bubbles with smooth or wavy skirts (sks or skw). Accordingly, different combinations of the Re, Eo, and Mo, which can be classified into three groups according to σ, ν, and g, are simulated, details of which are listed in Table III. The interested variable changes while the other two remain unchanged. Figure 3 shows the bubble shapes under different configurations in both 2D and 3D cases, in which the arrows designated using I, II, and III in Fig. 3(a) correspond to the simulation cases numbered A-E, F-J, and K-O in Table III, respectively. 

## 1. Surface tension

According to Eq. (7), changing surface tension σ alters both Eo and Mo but Re. A large surface tension tends to maintain the original shape of the bubble. As numbered A-E in Table III, the inertia force becomes predominant with decreasing surface tension (i.e., Eo changes from 0.327 35, 3.2735, 65.47, 654.7, to 935.29 while Re is always 25.892). As shown by arrow I in Figs. 3(b) and 3(c), the bubble shapes in 2D/3D keep the initial circle/sphere at $E o \ = \ 0 . 3 2 7 3 5 ,$ while they experience ellipse/oblate ellipsoid, ellipse/ellipsoidal cap, and skirted shapes at Eo = 3.2735, 65.47, and 654.7, respectively. A lower surface tension favors larger deformation, and the bubble deforms mainly at the bottom. 

## 2. Liquid viscosity

Similar to the shape selection along arrow I in Figs. 3(b) and 3(c), changing the liquid viscosity also causes different bubble shapes. As the liquid viscosity decreases, i.e., from F to J along arrow II (Re increases from 0.129 46, 2.5892, 25.892, 89.284, to 647.31 while Eo remains constant), the convection effect increases, and the bubble experiences more significant deformation. Under condition F, the bubble is circle/spherical, and no significant dent occurs at the bubble base. As Re increases, the indentation becomes more significant, accompanied by a thinner bubble tail and a more stretched bubble rear, and the bubble shape finally evolves into a skirted regime. 

The weak convection under large liquid viscosity makes the interface become more immobile or behave more like a rigid body, which reduces the rising velocity and makes bubbles exhibit little deformation. As the viscosity decreases, Re increases, and the effect of convection on bubble shape is strengthened, leading to large deformation and toroidal bubbles [see case H in Fig. 7(c)]. 

![](images/7e3f8971a65e0173de398995f06802ec20e853986c78c8aafaf1d6b98cd25233.jpg)


![](images/b8b96b338a0dd0c7b5c2bcf9948755208ac94931ce51ad461d8af7bad03b3174.jpg)



FIG. 3. (a) Locations of current simulation cases in the regime diagram proposed by Bhaga and Weber.8 Bubble shapes under different configurations in (b) 2D and (c) 3D cases. The arrows designated using I, II, and III in (a) correspond to the simulation cases numbered A-E, F-J, and K-O in Table respectively.


![](images/ea7d822e688e2757a4d3888f5c0ba3ec948f0a183a2d5931e9cd63d624ce1dcc.jpg)


![](images/9dbfd40e69129693ccfe788e69215ff1c31140e714aa6845c69effc4e74520db.jpg)



FIG. 4. Terminal velocity and length ratio of semimajor to semiminor axes vs density ratio in (a) 2D and (b) 3D cases, in which the insets show the bubble shapes by extracting $\phi = 0 . 5$ under different configurations.


## 3. Gravity acceleration

The bubble movement is driven by the volumetric buoyancy force. To investigate the effect of buoyancy force, the gravity acceleration is varied as numbered K-O in Table III [see arrow III in Fig. 3(a)]. The $R e ,$ Eo, and Mo numbers all increase with the gravity acceleration, and the bubbles exhibit significantly different shapes under different configurations, which agrees well with those predicted in Fig. 3(a). During the bubble rising, the flow forms a liquid jet beneath the bubble, causing a caplike or rimlike bubble shape. The jet effect increases with the Re number, as shown in $\mathrm { F i g s . } 3 ( \mathrm { b } )$ and 3(c), which agrees with those reported by Hua and Lou.59 

## 4. Density ratio

Besides the gravity acceleration, the buoyancy force is also dependent on the density difference. To study the effect of density difference, the density ratio $\rho _ { l } / \rho _ { g }$ is varied as 2, 5, 10, 20, 100, 200, and 1000 while other parameters maintain constant. Figure 4 shows the terminal velocity and length ratio of semimajor to semiminor axes vs the density ratio, in which the insets show the bubble shapes by extracting $\phi = 0 . 5$ under different configurations. The terminal velocity increases with the density ratio until reaching a steady state. When $\rho _ { l } / \rho _ { g } \geq 1 0 0$ , the bubble contours overlap with one another regardless of the domain dimension, indicating the density ratio has little effect on the terminal velocity and final shape, which agrees with those reported by Hua and Lou.59 

The buoyancy is proportional to the density difference $\Delta \rho = \rho _ { l }$ $- \ \rho _ { g }$ , and a larger density ratio corresponds to a larger buoyancy force and faster terminal velocity. Assuming $\rho _ { l } = 1 , \Delta \rho$ for each case is 0.5, 0.9, 0.95, 0.99, 0.995, and 0.999, respectively, and the points determined by both Re and Eo are all in the region designated using $^ { \omega } \mathrm { o e } ^ { \gamma }$ in Fig. 3(a). It is noted that the point with $\Delta \rho = 0 . 5$ locates near the boundary between the regions $\kappa _ { \boldsymbol { \mathsf { S } } } ^ { \bullet }$ and $^ \alpha { _ { 0 } e , ^ { \scriptscriptstyle 3 } }$ and the bubble has the smallest length ratio, i.e., approaching circular/spherical shape. $\mathrm { A s } \Delta \rho ( \mathrm { o r } \rho _ { l } / \rho _ { g } )$ increases, both the rising velocity and the length ratio of semimajor to semiminor axes increase due to increasing buoyancy effect, but such an effect is only significant when density ratio is low, i.e., $\rho _ { l } / \rho _ { g } < 1 0 0$ . 

## C. Bubble drag

In the experiments reported by Bhaga and $\mathrm { w e b e r , ^ { 8 } }$ the drag coefficient during bubble rising is a function of Re. For the fluid with high Morton number $( M o > \breve { 4 } \times 1 0 ^ { - 3 } ) ,$ the dependence of the drag coefficient $C _ { D }$ on Re can be described as8 

$$
C _ {D} = \left[ (2. 6 7) ^ {0. 9} + (1 6 / R e) ^ {0. 9} \right] ^ {1 / 0. 9}. \tag {8}
$$

Joseph60 proposed a hyperbolic drag law based on the potential flow, i.e., 

$$
C _ {D} = 0. 4 4 5 (6 + 3 2 / R e). \tag {9}
$$

Figure 5 shows the computed relationship between the drag coefficient and the Reynolds number compared with those proposed by Bhaga and Weber8 and Joseph.60 The drag coefficient is computed from a force balance, i.e., $\dot { C } _ { D } ~ = ~ { \pi } \Delta \rho g \dot { d _ { 0 } } / \big ( 2 \rho _ { l } U ^ { 2 } \big )$ and $C _ { D } \ = \ 4 \Delta \rho g d _ { 0 } / \big ( 3 \rho _ { l } U ^ { 2 } \big )$ in 2D and 3D, respectively.42,60 A reasonable agreement with the empirical correlation8 and potential flow solution60 is achieved for the current simulations. The consistency depends on both the Reynolds number and the domain dimension. It can be observed that the 3D bubble dynamics matches well with the correlations compared with the 2D simulations, especially at larger Re (approaching 100). 

![](images/10ee3cb8f174ff22eb45458c86b12c63533f08ccfeba38f9e5541553b3602090.jpg)



FIG. 5. Drag coefficient $C _ { D }$ vs the Reynolds number Re compared with the empirical correlation proposed by Bhaga and Weber8 and potential flow solution derived by Joseph.


![](images/89a2410f92efb423d000729da209e6dfe991fa5952784ea1cfb9c55f17d2de00.jpg)



FIG. 6. Evolution of the rising bubble in a 2D case. The parameter sets in (a), (b), (c), and (d) correspond to cases K, B, H, and N in Fig. 3 (left to right), respectively. For cases K, B, and H in (a), (b), and (c), the snapshots are extracted by setting φ = 0.5 after 0, 10 000, 20 000, 30 000, and 40 000 steps from the bottom up, respectively, while for case N in (d), the time is after 0, 500, 800, 1000, and 1500 steps from the bottom up, respectively.


## IV. RESULTS AND DISCUSSION

Four different cases, including single bubble, two bubbles, multiple bubbles, and the interaction between the bubble and the obstacle, are considered to validate the current numerical scheme and explore the interaction between gas and liquid. The parameter setting is consistent with those in Sec. III unless mentioned otherwise. 

## A. Single bubble

Four typical cases (i.e., K, B, H, and N in Fig. 3) are used to illustrate the bubble dynamics under different configurations. Figures 6 and 7 show the evolution of the rising bubble in 2D and 3D cases, respectively, and the first row in Fig. 7 shows the views along the [001], [010], and [111] directions. For case K, the surface tension force is predominant compared with the inertial force, and the bubble exhibits circular/spherical shape throughout the rectilinear trajectory. When the surface tension force becomes small [see case B in Figs. 6(b) and 7(b)], the bubble shape experiences a transition from the initial circular/spherical to the steady elliptical/oblate ellipsoidal shape. The length along the rising direction (i.e., y axis) is compressed, while those along the other two axes are stretched. In cases H and N, the shape evolution is more affected by the inertial force, and a dented region exhibits at the bottom. 

![](images/dabc533436c9946833f25eb551f57e4e0593f377d60adff25f1724f1548a2d00.jpg)



FIG. 7. Evolution of the rising bubble in a 3D case. The first row shows the three views along the [111], [010], and [001] directions. The parameter sets in (a), (b), (c), and (d) correspond to cases K, B, H, and N in Fig. 3 (left to right), respectively. For cases K, B, and H in (a), (b), and (c), the snapshots are extracted by setting φ = 0.5 after 0, 10 000, 20 000, 30 000, and 40 000 steps from the bottom up, respectively, while for case N in (d), the time is after 0, 500, 800, 1000, and 1500 steps from the bottom up, respectively.


![](images/c74c9d2019850e2d6ce114b01e37f9660507bc507577482454f550987fcb9d73.jpg)



FIG. 8. Distribution of the velocity field corresponding to case H in (a1) and (a2) 2D and (b1) and (b2) 3D, in which the bubble contours are depicted using black lines and the arrows denote the velocity vector.


![](images/8a70c17c883926fce3ea3cb19ab401425c870fef5fb2dcf07ba306ae995c788f.jpg)



FIG. 9. Rising velocity of the bubble vs time for case B, in which dt is the time step during simulation.


For a large Eo number $[ E o = 3 2 7 . 3 5 ,$ see Figs. 6(c) and 6(d)], the bubble undergoes larger deformation, exhibiting skirted shapes in the late stage. An increasingly large dented area exhibits as the bubble rises, accompanied by a more stretched bubble rear and a thinner tail. At the later stage, the 2D bubble breaks up, exhibiting a symmetric thin film of gas trailing behind the bubble rim. Thin satellite bubbles known as skirts trail the large caps, which has been observed in the experiments performed by Bhaga and Weber8 and simulations reported by Yu and Fan7 and Hua and Lou.59 

In the 3D case shown in Figs. 7(c) and 7(d), the bottom of the spherical bubble evolves into a dimpled form, indicating the rising velocity of the lower surface is larger than that of the upper surface, 

![](images/c308d7e73663cda31ce2adff7a2548e880b76f9c292ba54532dfc97e4283ed36.jpg)



FIG. 10. Evolution of the two bubbles under different configurations. (a), (b), and (c) correspond to cases K, B, and H, respectively. The initial state is shown in the insets in the first column. For case K in (a1), (a2), (a3), and (a4), the snapshots are extracted by setting $\dot { \phi } = 0 . \dot { 5 }$ after 20 000, 40 000, 50 000, and 80 000 steps from left to right, respectively, while for case B in (b1), (b2), (b3), and (b4) and case H in (c1), (c2), (c3), and (c4), the time is after 10 000, 20 000, 30 000, and 40 000 steps from left to right, respectively.


which agrees well with those reported by Amaya-Bower and Lee.58 The lower surface involutes, generating a high-speed liquid jet which directs along the reverse gravity (see Fig. 8). As the bubble rises, the jet moves toward the opposite surface until forming a toroidal bubble which occurs upon the collapse of transient cavities.61 After jet penetration, the bubble continues to collapse, and the cross-sectional area of the top penetration region increases to the diameter of the liquid jet, whilst the rest of the bubble keeps shrinking; that is, the bubble undergoes a topological transition from the initial sphere to torus. The bubble with a vortex ring has been validated experimentally and theoretically by observing the bubble approaching a rigid boundary.6 2 

Figure 8 shows the distribution of the velocity field corresponding to case H, in which the bubble contours are depicted using black 

![](images/c012fffbfb3f6b5d556dad2a3cf2e67e9c5d530c963d51f0c463e7580f505bff.jpg)



FIG. 11. Evolution of the flow field at the central slice according to Fig. 10. (a), (b), and (c) correspond to cases K, B, and H, respectively. (a1)–(a4) correspond to (a1)–(a4) in Fig. 10. (b1)–(b4) correspond to (b1)–(b4) in Fig. 10. (c1)–(c4) correspond to (c1)–(c4) in Fig. 10.


![](images/6d99c8a02b32b6b38ca00de62be9f1a74b44ecbccf5f75f375ed8eb4d5118d33.jpg)


![](images/e3eb0672650d1ab92017691799ed7247f2e2f7480504a8f5251ee68a12d1a197.jpg)


![](images/bb591d3a6fc67d6e6ba1cb74fbe66028901932e7639a0b8c97df41644bac05a2.jpg)



FIG. 12. Position of the bubble top and bottom vs time, in which the coalescence of the two bubbles is marked by the dashed lines. (a), (b), and (c) correspond to cases K, B, and $\mathsf { H } ,$ respectively. The insets show the bubble shape at different time.


lines. In the 2D case shown in Figs. $8 ( \mathrm { a } _ { 1 } )$ and $8 ( \mathrm { a } _ { 2 } )$ , both the flow field and the bubble shape exhibit axisymmetric distribution, while those in 3D [see $\mathrm { F i g s . ~ 8 ( b _ { 1 } ) }$ and $8 ( \mathsf { b } _ { 2 } ) ]$ are in radial symmetry. The difference of the flow field induced by the domain dimension causes different bubble shapes. The liquid in 2D reaches the downstream by necessarily crossing the two horizontal tips of the bubble. However, the liquid in 3D can bypass the bubble, which increases the flow flexibility and changes the convection effect, i.e., equivalently weakening the effect of convection on the side of the bubble. Besides, different from the 2D case that the velocity can obtain two extremes along the symmetry axis [see the two red regions in $\mathrm { F i g } . \ 8 ( { \mathrm { a } } _ { 2 } ) ]$ , a local high-velocity region concentrates as a whole in the middle in 3D [see Fig. $8 ( \mathrm { b } _ { 2 } ) ]$ . The resultant liquid jet forms the vortex round the bubble ring, which wraps the bubble to facilitate the bubble rising. 

The influence of the domain dimension is examined by comparing the rising velocity between 2D and 3D cases. Figure 9 shows the rising velocity vs time for case B. The two curves (for 2D and 3D) change in a similar manner, i.e., increasing rapidly until approaching a constant and then decreasing after approaching the domain top. The middle state is not totally flat, which can be considered to be due to the bubble deformation during the rising process. When the bubble reaches the domain top, slight velocity fluctuation occurs due to the bouncing of the bubble from the top wall. 

## B. Two bubbles

Two equal sized bubbles are initially placed at the domain with a separation distance of $1 . 2 d _ { 0 }$ between their centers. Three typical parameter sets including cases K, B, and H are simulated to study the bubble coalescence. Figure 10 shows the evolution of the two bubbles under different configurations, in which the initial state is shown in the insets in the first column. The two bubbles rise with different velocities and are always aligned vertically along the symmetry axis without lateral offset. The trailing bubble undergoes more significant deformation because it is more affected by the hydrodynamic interaction between the two bubbles. The thin liquid film between bubbles is squeezed during collision and then ruptured due to the coalescence. After that, the merged bubble returns to the situation of a single bubble rising. 

![](images/2f06f6957886dbc20d120b4354a32e8ee8cee268846215a7e7f22c44f981020d.jpg)



FIG. 13. Effect of bubble spacing on the coalescence of the two bubbles, in which “coal.” and $\mathsf { s e p a } . \mathsf { i n }$ the legend denote the “coalescence” and “separate” and dt and $d _ { 0 }$ are the time step and the bubble diameter, respectively.


In case K, the surface tension dominates during bubble rising, and the two bubbles exhibit little deformation except the peanutlike shape at the merging stage [see Fig. 10(a2)]. As the surface tension decreases, the initial spherical shape cannot be maintained, and the shape of the two bubbles becomes disklike [see Fig. 10(b4)] or caplike [see Fig. 10(c3)]. The trailing bubble has larger deformation and moves toward the wake zone of the leading bubble. For case B, a local high-curvature region develops at the top of the trailing bubble, while the leading one exhibits an oblate ellipsoidal cap shape [see Fig. 10(b2)]. The bottom surface of the leading bubble becomes flat in the presence of the trailing bubble. Once the peak of the taperlike trailing bubble encounters the leading one, coalescence occurs and the radius of the contact line increases until approaching the horizontal radius of the leading one. Finally, a new bigger bubble rises in a disklike shape due to buoyancy [see Fig. 10(b4)]. 

In case H, the bubbles rise with significant deformation due to smaller surface tension (larger Eo). The trailing bubble is immersed in the wake of the leading one and interacts with the wake flow, which causes elongated deformation along the vertical symmetry axis. The trailing bubble behaves like a jet into the dimple of the leading one. A similar phenomenon has been observed in experiments reported by Brereton and Korotney63 and simulated using the LBM by Takada et al.64 and Cheng et al.55 Then, the two bubbles merge into a caplike whole and evolve into the toroidal shape at last. Before the bubbles collapse, less fluid is drawn into the jet region due to the interaction between the two bubbles, which results in a weaker jet compared with that in the one-bubble case. 

Figure 11 shows the evolution of the flow field at the central slice corresponding to Fig. 10. The flow velocity evolves with the rising bubble and obtains local extreme near the bubbles. A pronounced feature is that the flow velocity near the trailing bubble is always larger than that near the leading one regardless of the flow pattern. Influenced by the wake flow of the leading one, the trailing bubble is attracted toward the leading one until coalescing into a larger one. Besides, the flow pattern in case K changes less significantly compared with those in cases B and H, which is consistent with the shape difference among different configurations in Fig. 10. 

![](images/80a8b644f0b2d5f98ce0c7bf0f065ae0492015a9e0ac8b3d7c715ede8db61f44.jpg)



FIG. 14. Evolution of the 20 bubbles under different configurations. (a), (b), and (c) correspond to cases K, B, and H, respectively. The initial state is shown in the insets in $( \mathsf { b } _ { 4 } ) .$ . For case K in (a1), (a2), (a3), and (a4) and case B in (b1), $( \mathsf { b } _ { 2 } ) , ( \mathsf { b } _ { 3 } ) ,$ and $( \mathsf { b } _ { 4 } ) _ { : }$ , the snapshots are extracted by setting $\phi = 0 . 5$ 5 after 10 000, 50 000, 100 000, and 150 000 steps from left to right, respectively. For case H in (c1), (c2), $( \mathsf { c } _ { 3 } ) ,$ and $\left( \mathsf { c } _ { 4 } \right)$ , the time is after 10 000, 20 000, 40 000, and 60 000 steps from left to right, respectively.


Figure 12 shows the position of the bubble top and bottom vs time, in which the insets show the bubble shape at different time. The coalescence of the two bubbles is marked by the dashed lines. Before coalescence, the trailing bubble rises faster with larger deformation, and the average rising velocity is 0.065, 0.196, and 0.293 for cases K, B, and H respectively, which is larger than that of the leading one (i.e., 0.051, 0.166, and 0.224). The bubbles experience morphological transition during coalescence and then behave like a single bubble. For case K in Fig. 12(a), the spacing between the bubble top and bottom remains unchanged, while for case B in Fig. 12(b), the spacing first decreases and then increases. For case H in Fig. 12(c), the coalesced bubbles evolve into toroidal shape, which is similar to that in Fig. 7(c). 

Besides being influenced by the parameter setting (i.e., Re, Eo, and Mo numbers), the interaction and coalescence of two bubbles are also dependent on the initial bubble spacing. The trailing bubble is fixed at $y = 0 . 2 Y ,$ and the position of the leading one is varied to adjust the bubble spacing, where $Y = 1 . 2 d _ { 0 }$ is the domain height. Figure 13 shows the effect of bubble spacing on the coalescence. Figure 13 can be divided into two parts by the dashed lines with respect to the coalescence. On the left, the time at the y axis denotes the starting time when the coalescence occurs, and it increases with the bubble spacing. Once the initial bubble spacing is larger than the critical value indicated by the dashed lines (i.e., 1.32d0, 1.8d0, and 3.12d0 for cases K, B, and H respectively), no coalescence can be observed during the rising process. On the right, the time denotes the moment when the leading bubble meets the top wall, and it decreases with the initial bubble spacing due to less rising distance for the leading one. As the initial bubble spacing increases, the effect of the wake flow on the trailing bubble is weakened, and the two bubbles tend to rise in a similar trend. In addition, the difference among the three cases further indicates the dependence of bubble dynamics on three dimensionless numbers. 

## C. Multiple bubbles

Twenty equal sized bubbles are randomly initialized in a larger domain of $1 0 d _ { 0 } \times 2 0 d _ { 0 } \times 1 0 d _ { 0 }$ . Different from the two-bubble case in Sec. IV B, the interaction between multiple bubbles comprises both the coaxial coalescence and oblique coalescence. Different configurations under the three typical parameter sets, i.e., numbered K, B, and H, are compared to study the interaction mechanism. 

![](images/c55bd07c9bbc5d958f4363541256df133249640a5590d95e40244c523a4e0f3c.jpg)



FIG. 15. Evolution of the flow field at the central slice according to Fig. 14. (a), (b), and (c) correspond to cases K, B, and H, respectively. $( \mathsf { a } _ { 1 } ) - ( \mathsf { a } _ { 4 } )$ correspond to $( { \sf a } _ { 1 } ) -$ $( \mathsf { a } _ { 4 } ) \mathsf { i n F i g . } 1 4 . ( \mathsf { b } _ { 1 } ) – ( \mathsf { b } _ { 4 } )$ correspond to $( { \mathsf { b } } _ { 1 } ) – ( { \mathsf { b } } _ { 4 } ) \mathrm { i n } \ F i g . \ 1 4 .$ $( \mathsf { c } _ { 1 } ) \mathsf { - } ( \mathsf { c } _ { 4 } )$ correspond to $( \mathsf { c } _ { 1 } ) - ( \mathsf { c } _ { 4 } ) \mathrm { i n } \mathsf { F i g . } \mathsf { \dot { 1 } } 4 .$


Figure 14 shows the evolution of the 20 bubbles under different configurations, in which the initial state is shown in the insets in Fig. 14(b4). All bubbles rise due to buoyancy force but with different velocities and tend to aggregate to form clusters or chimneys.65 The wake flow of each bubble influences the motion of other neighboring bubbles, and such an effect varies with both the distance between two bubbles and the parameter setting. The coalescence occurs between the two bubbles which are not far from each other; otherwise, the bubbles rise like an isolated bubble but with different levels of deformation, especially for those leading bubbles in the upper part of the domain. 

In case K, the bubbles remain spherical, and the coalescence occurs when approaching the top boundary. In case B, the bubble shape becomes disklike or caplike, which depends on the physical position. In case H, the bubbles evolve into vortex chains, which has been discussed by Brücker66 and is simulated for the first time to our best knowledge. The vortex chains have complicated geometric topologies, and the deformation becomes larger over time. Besides, regardless of the parameter setting, the alignment effect induced by the interaction among the bubbles deforms the trailing bubbles considerably and inclines those with respect to the vertical direction. 

![](images/27fbf71e7064e4ebc1e8445bd5f6a1bb11b06ddfb99690ce6a7bad8bacb27654.jpg)



FIG. 16. Evolution of the bubble shape and flow field under another two boundary conditions in case B. (a1)–(a4) and (c1)–(c4) correspond to the same time in Figs. 14(b1)– 14(b4). (b1)–(b4) and (d1)–(d4) correspond to the same time in Figs. 15(b1)–15(b4). In (a) and (b), a periodic boundary condition is set at four side walls for both phase field and velocity, while the boundary conditions at both top and bottom boundaries are kept the same as those in Fig. 14. In (c) and (d), a periodic boundary condition is set at all walls. The parameter configuration is set to be the same as case B, and the snapshots are extracted by setting φ = 0.5 after 10 000, 50 000, 100 000, and 150 000 steps from left to right, respectively.


Figure 15 shows the evolution of the flow field at the central slice corresponding to Fig. 14. The distribution of the flow field is extremely nonuniform, and the number of bubbles at the slice varies with time, indicating that the rising path of multiple bubbles is a space curve. The flow velocity always obtains the local extremes near the bubbles, and the flow tends to be gentle as the bubbles meet the top side. 

As discussed in Sec. III A, the wall effect, or equivalently the boundary effect, has remarkable influence on the bubble dynamics. Figure 16 shows the bubble evolution under another two boundary conditions, and the parameter configuration is set to be the same as case B. In Figs. 16(a1)–16(b4), a periodic boundary condition is set at four side walls for both phase field and velocity, while the boundary conditions at both top and bottom boundaries are kept the same as those in Fig. 14. In Figs. 16(c1)–16(d4), a periodic boundary condition is set at all walls to simulate an infinite reservoir. A significant difference can be observed in the bubble shape and the distribution of the flow field. In particular, the bubbles tend to become more spherical under the periodic boundary condition due to less influence from the domain boundary. But under the no-slip condition, the velocity at the wall is forced to be zero, which equivalently imposes an extra force to deform the bubbles, especially for those near the domain boundary [see Figs. 14 and 16(a1)–16(a4)]. 

It is noted that physical coalescence occurs within the range of the van der Waals force, which is about 10 nm and orders of magnitude smaller than the interface thickness in the diffuse interface methods.37 The interface with some mesh width causes a so-called premature coalescence, which has been examined in simulating droplet collisions and is one of the major difficulties in simulating coalescence behavior.67 In the present work, the bubble coalescence criterion is simplified by only considering the geometric constraint, which is adopted in other numerical schemes such as the front tacking method.1 

![](images/59563264de895b8367f7da90632eb87fa2333ab7adfbc86e8bc257be3ddf0b39.jpg)



FIG. 17. Evolution of the rising bubble in the presence of the cuboid static obstacle under different configurations, and four snapshots of the bubble deformation at different instants exhibit on the right of each subfigure. For case K in (a), (b), and (c), the snapshots are extracted by setting φ = 0.5 after 0, 60 000, 100 000, 120 000, 150 000, 160 000, 180 000, and 200 000 steps from the bottom up, respectively. The last six snapshots overlap with one another, and the four snapshots exhibited on the right of each subfigure are extracted after 100 000, 120 000, 150 000, and 160 000 steps, respectively. For case B in (d), (e), and (f), the time is after 0, 20 000, 30 000, 40 000, 50 000, 60 000, 70 000, and 80 000 steps from the bottom up, respectively. The four snapshots exhibited on the right of each subfigure are extracted after 30 000, 40 000, 50 000, and 60 000 steps, respectively. For case H in (g), (h), and (i), the time is after 0, 10 000, 20 000, 30 000, 40 000, 50 000, 60 000, and 70 000 steps from the bottom up, respectively. The four snapshots exhibited on the right of each subfigure are extracted after 20 000, 30 000, 40 000, and 50 000 steps, respectively.


## D. Interaction between bubble and obstacle

The impingement between the bubble and the obstacle is a common phenomenon in the gas-liquid system, which involves buoyancy-driven rising motion and splitting of bubbles. To investigate the impingement mechanism, both the size of the obstacle and the parameter setting during bubble rising are changed to observe the shape evolution. 

Figure 17 shows the evolution of the rising bubble in the presence of the obstacle under different configurations, and four snapshots of bubble deformation at different instants exhibit on the right of each subfigure. A cuboid obstacle with the size of a × a × b is placed at the center of the confined domain. $a = a = 5 0$ and b is varied as 5, 20, and 40 in the lattice unit. Each column corresponds to the same obstacle size, while each row is with the same parameter set, i.e., cases K, B, and H from top to bottom, respectively. Both the initialization of the bubble and the boundary conditions are the same as those in Fig. 7, and a bounce back scheme1 is employed at the surface of the obstacle. 

The presence of obstacles changes the distribution of the flow field and thus affects the bubble evolution. At the beginning, the bubble rises due to buoyancy and then touches and interacts with the static obstacle, causing different evolution processes of the bubble shape. When the size of the obstacle changes, whether the obstacle can block the rising bubble is dependent on the parameter setting. Taking Figs. 17(d)–17(f) for instance, the bubble splits into two symmetric separate bubbles when touching the smaller obstacle due to the vertical shearing effect [see Figs. 17(d) and 17(e)] but maintains as a whole under larger obstacles [see Fig. 17(f)]. In the first row, the bubble shape is spherical before contacting the obstacle. After it collides with the bottom surface of the obstacle, the bubble is squeezed due to being completely blocked by the obstacle. Such a blocking effect looks more significant under smaller obstacles, where the obstacle seems to be inlaid in the bubble [see Fig. 17(a)]. In the bottom row, Re equals that in the second row, but Eo is a 100 times larger, i.e., a much weaker surface tension effect. Regardless of the obstacle size, the bubble in the bottom row breaks into separate parts after impinging the obstacle, the number of which doubles as the size of the obstacle increases. The separate smaller bubbles exhibit caplike shape and rise along the side walls of the obstacle. 

![](images/bc930b1c0e36064485ea0e6cf20aef6bf517eba96e5e9cf2adb0bf4da8b2f8e4.jpg)



FIG. 18. Evolution of the flow field at the central slice according to the four snapshots selected using dashed lines in Figs. 17(b), 17(e), and 17(h), respectively. (a), (b), and (c) correspond to cases K, B, and H, respectively. (a1)–(a4) correspond to the four snapshots in Fig. 17(b). (b1)–(b4) correspond to the four snapshots in Fig. 17(e). (c1)–(c4) correspond to the four snapshots in Fig. 17(h).


Figures 18(a)–18(c) show the evolution of the flow field at the central slice corresponding to the four snapshots selected using dashed lines in Figs. 17(b), 17(e), and 17(h), respectively. The deformation of the bubbles under flow field can be clearly distinguished. Blocked by the static obstacle where the flow velocity is zero, the bubbles exhibit a symmetric structure regardless of the occurrence of the bubble splitting. 

Figure 19 shows the rising velocity vs time corresponding to cases K, B, and H in Figs. 17(b), 17(e), and 17(h), respectively. For case K, the bubble is totally blocked by the obstacle, and the rising velocity is reduced to zero. For cases B and H, the bubble breaks up after impacting on the obstacle, and the rising velocity decreases. The surface tension greatly alters the morphology of the split children bubbles. The shape of the two new children bubbles in case B is similar to that of the parent, and the rising velocity increases until approaching the top wall. But for case H, the surface tension is too small to integrate the children bubbles, and more children bubbles are generated during the rising process, which is accompanied by the decrease of the volume of the main body. 

The width of the obstacle b is changed to characterize the blocking effect of the cuboid obstacle. Figure 20 illustrates the blocking effect as a function of the width b in the three cases K, B, and H, in which such an effect manifests the breakup in the upper half and blocking in the lower half. The bubble in case H breaks up when impacting on the obstacle regardless of the obstacle width, while the bubble in case K is totally blocked. In case B, the interaction between the bubble and the obstacle changes from breakup to being blocked as b increases. Different interactions in the three typical cases further indicate the importance of the buoyancy force and surface tension. When the surface tension dominates, i.e., case K (a smaller Eo), the bubble exhibits little deformation. As surface tension decreases or the buoyancy force increases, the bubble deformation increases, and the bubble is sensible to the presence of obstacles; for example, the bubble in case H splits even if the width of the obstacle is only 0.08d0. 

![](images/a3e069b4187d51902560967fddd8efd3ac06363460a30e5f6746ee1570119d09.jpg)



FIG. 19. Rising velocity of the bubble vs time for cases K, B, and H in Figs. 17(b), 17(e), and 17(h), respectively, in which dt is the time step during simulation.


![](images/a8eb9f09f43d2c8d2c4caae4e7f302f6183febdf87325f5c5eca463f74c45d96.jpg)



FIG. 20. Blocking effect as a function of the obstacle width for the three cases K, B, and H. The blocking effect comprises the breakup in the upper half and being blocked in the lower half.


Both the obstacle size and the parameter setting affect the bubble deformation. When the rising bubble touches the obstacle, it collides, splits, or adheres to the bottom surface of the obstacle. Because the main purpose of this work is to investigate the deformation mechanism after interacting with the obstacle, the shape of the obstacle is simplified and other factors, such as the interfacial tension and contact angles between fluids and obstacle surface, are neglected, which will be further studied in a future work. 

## V. CONCLUSIONS

A conservative phase-field model is employed to track the interface during bubble rising, and the lattice Boltzmann model is used to solve the two-phase flow induced by the buoyancy force and the surface tension force. To facilitate 3D simulations, a parallel and adaptive mesh refinement algorithm is developed to reduce the computing overhead. The simulated terminal velocity and the corresponding bubble shape are compared with the shape chart based on experimental observations, which validates the robustness of the present numerical scheme in simulating multiphase flow. Both 2D and 3D analyses are performed to establish the dependence of the numerical parameters on the bubble rising including domain size, surface tension, liquid viscosity, gravity, and density ratio. It is found that the simulated results are in good agreement with those published in the literature. 

The simulation of multiple bubbles validates the capability of the numerical scheme in modeling the complicated topological change including the coaxial and oblique bubble coalescence. 

Bubbles rise due to buoyancy but with different velocities. The leading bubbles behave like isolated ones before coalescence, while the trailing ones undergo more significant deformation, especially when being trapped in the wake of the leading bubble. The effect of the initial bubble spacing on the coalescence of the two bubbles is investigated. The interaction between the rising bubble and the static obstacle is evaluated by changing both obstacle size and numerical parameters. When the bubbles can be completely blocked by the obstacle is quantified in terms of the obstacle width. Because of the blocking effect induced by the predefined static obstacle, squeeze deformation and bubble splitting occur. 

It is worth stressing that the bubble evolution under any configuration is largely dependent on the parameter setting or the three dimensionless numbers including the Reynolds number, Eötvös number, and Morton number. Different combinations of these numbers correspond to different bubble dynamics and thus lead to various flow features. To approach real situations and actual engineering practice, it is better that more aspects including bubble nucleation and growth are considered, which will be performed in a forthcoming investigation. 



28A. Fakhari and M. H. Rahimian, “Phase-field modeling by the method of lattice Boltzmann equations,” Phys. Rev. E 81, 036707 (2010). 



## ACKNOWLEDGMENTS



29A. Fakhari, M. Geier, and T. Lee, “A mass-conserving lattice Boltzmann method with dynamic grid refinement for immiscible two-phase flows,” J. Comput. Phys. 315, 434 (2016). 



This work was financially supported by the Joint Funds of the National Natural Science Foundation of China (Grant No. U1537202), the Tsinghua-General Motors International Collaboration Project (Grant No. 20153000354), the Tsinghua University Initiative Scientific Research Program (Grant No. 20151080370), and the Tsinghua Qingfeng Scholarship (Grant No. THQF2018-15). The authors would also like to thank the National Laboratory for Information Science and Technology in Tsinghua University for access to supercomputing facilities. 



30H. Liang, B. C. Shi, and Z. H. Chai, “An efficient phase-field-based multiplerelaxation-time lattice Boltzmann model for three-dimensional multiphase flows,” Comput. Math. Appl. 73, 1524 (2017). 



## APPENDIX: DERIVATION OF THE CONSERVATIVE PHASE-FIELD MODEL



31Y. Sun and C. Beckermann, “Sharp interface tracking using the phase-field equation,” J. Comput. Phys. 220, 626 (2007). 



The phase-field equation according to Fick’s second law68 is expressed as 



32J. W. Cahn and J. E. Hilliard, “Free energy of a nonuniform system. I. Interfacial free energy,” J. Chem. Phys. 28, 258 (1958). 



$$
\partial_ {t} \phi = - \nabla \cdot (- M \nabla \phi). \tag {A1}
$$



33N. Takada, J. Matsumoto, and S. Matsumoto, “Phase-field model-based simulation of motions of a two-phase fluid on solid surface,” J. Comput. Sci. Technol. 7, 322 (2013). 



The phase field variable at equilibrium changes in a hyperbolic tangent manner [see Eq. (2)] and is independent of the flux type.34,53 A phase-separation flux js is incorporated into Eq. (A1), which counteracts the diffusive flux at equilibrium, 



34M. Geier, A. Fakhari, and T. Lee, “Conservative phase-field lattice Boltzmann model for interface tracking equation,” Phys. Rev. E 91, 063309 (2015). 



$$
\boldsymbol {j} _ {s} = M \nabla \phi^ {e q} = \frac {M}{4 \delta} \left[ 1 - \tanh ^ {2} \left(\frac {\varepsilon}{2 \delta}\right) \right] \frac {\nabla \phi}{| \nabla \phi |} = \frac {M}{\delta} (\phi - \phi^ {2}) \frac {\nabla \phi}{| \nabla \phi |}. \tag {A2}
$$



35S. M. Allen and J. W. Cahn, “Mechanisms of phase-transformations within miscibility gap of Fe-rich Fe-Al alloys,” Acta Metall. 24, 425 (1976). 



The velocity v is coupled into Eq. (A1) to account for the convective flux.5 3 



36H. Z. Yuan, C. Shu, Y. Wang, and S. Shu, “A simple mass-conserved level set method for simulation of multiphase flows,” Phys. Fluids 30, 040908 (2018). 



## REFERENCES



37T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, and E. M. Viggen, The Lattice Boltzmann Method Principles and Practice (Springer, Cham, Switzerland, 2017). 





38A. Zhang, J. Du, Z. Guo, and S. Xiong, “Lamellar eutectic growth under forced convection: A phase-field lattice-Boltzmann study based on a modified Jackson-Hunt theory,” Phys. Rev. E 98, 043301 (2018). 





1G. Kaehler, F. Bonelli, G. Gonnella, and A. Lamura, “Cavitation inception of a van der Waals fluid at a sack-wall obstacle,” Phys. Fluids 27, 123307 (2015). 





39J. Töelke, S. Freudiger, and M. Krafczyk, “An adaptive scheme using hierarchical grids for lattice Boltzmann multi-phase flow simulations,” Comput Fluids 35, 820 (2006). 





2V. H. Arakeri, “Viscous effects on the position of cavitation separation from smooth bodies,” J. Fluid Mech. 68, 779 (1975). 





40C. Janssen and M. Krafczyk, “A lattice Boltzmann approach for free-surfaceflow simulations on non-uniform block-structured grids,” Comput. Math. Appl. 59, 2215 (2010). 





3A. Zhang, S. Liang, Z. Guo, and S. Xiong, “Determination of the interfacial heat transfer coefficient at the metal-sand mold interface in low pressure sand casting,” Exp. Therm. Fluid Sci. 88, 472 (2017). 





41H. Z. Yuan, Y. Wang, and C. Shu, “An adaptive mesh refinement-multiphase lattice Boltzmann flux solver for simulation of complex binary fluid flows,” Phys. Fluids 29, 123604 (2017). 





42A. Gupta and R. Kumar, “Lattice Boltzmann simulation to study multiple bubble dynamics,” Int. J. Heat Mass Transfer 51, 5192 (2008). 





4S. Wang, Z. P. Guo, X. P. Zhang, A. Zhang, and J. W. Kang, “On the mechanism of dendritic fragmentation by ultrasound induced cavitation,” Ultrason. Sonochem. 51, 160 (2019). 





43R. Han, S. Li, A. M. Zhang, and Q. X. Wang, “Modelling for three dimensional coalescence of two bubbles,” Phys. Fluids 28, 062104 (2016). 





5A. W. Woods, “The dynamics of explosive volcanic-eruptions,” Rev. Geophys. 33, 495, https://doi.org/10.1029/95rg02096 (1995). 





44G. Chen, X. Huang, A. Zhang, and S. Wang, “Simulation of three-dimensional bubble formation and interaction using the high-density-ratio lattice Boltzmann method,” Phys. Fluids 31, 027102 (2019). 





6M. Blander, “Bubble nucleation in liquids,” Adv. Colloid Interface Sci. 10, 1 (1979). 





45C. Zhang, J. Li, L. Luo, and T. Qian, “Numerical simulation for a rising bubble interacting with a solid wall: Impact, bounce, and thin film dynamics,” Phys. Fluids 30, 112106 (2018). 





7Z. Yu and L. Fan, “An interaction potential based lattice Boltzmann method with adaptive mesh refinement (AMR) for two-phase flow simulation,” J. Comput. Phys. 228, 6456 (2009). 





46A. Zhang, J. Du, Z. Guo, Q. Wang, and S. Xiong, “Conservative phase-field method with a parallel and adaptive-mesh-refinement technique for interface tracking,” Phys. Rev. E (submitted). 





8D. Bhaga and M. E. Weber, “Bubbles in viscous liquids: shapes, wakes, and velocities,” J. Fluid Mech. 105, 61–85 (1981). 





47Z. Guo and S. M. Xiong, “On solving the 3-D phase field equations by employing a parallel-adaptive mesh refinement (Para-AMR) algorithm,” Comput. Phys. Commun. 190, 89 (2015). 





9T. Maxworthy, C. Gnann, M. Kurten, and F. Durst, “Experiments on the rise of air bubbles in clean viscous liquids,” J. Fluid Mech. 321, 421 (1996). 





10P. Cui, Q. X. Wang, S. P. Wang, and A. M. Zhang, “Experimental study on interaction and coalescence of synchronized multiple bubbles,” Phys. Fluids 28, 012103 (2016). 





48J. Du, A. Zhang, Z. Guo, M. Yang, M. Li, and S. Xiong, “Atomistic determination of anisotropic surface energy-associated growth patterns of magnesium alloy dendrites,” ACS Omega 2, 8803 (2017). 





11D. M. Sharaf, A. R. Premlata, M. K. Tripathi, B. Karri, and K. C. Sahu, “Shapes and paths of an air bubble rising in quiescent liquids,” Phys. Fluids 29, 122104 (2017). 





49A. Zhang, Z. Guo, and S. Xiong, “Phase-field-lattice Boltzmann study for lamellar eutectic growth in a natural convection melt,” China Foundry 14, 373 (2017). 





12D. W. Moore, “The rise of a gas bubble in a viscous liquid,” J. Fluid Mech. 6, 113 (1959). 





50A. Zhang, Z. Guo, and S. M. Xiong, “Eutectic pattern transition under different temperature gradients: A phase field study coupled with the parallel adaptivemesh-refinement algorithm,” J. Appl. Phys. 121, 125101 (2017). 





13T. Bonometti and J. Magnaudet, “Transition from spherical cap to toroidal bubbles,” Phys. Fluids 18, 052102 (2006). 





51E. Olsson, G. Kreiss, and S. Zahedi, “A conservative level set method for two phase flow II,” J. Comput. Phys. 225, 785 (2007). 





14B. A. Puthenveettil, A. Saha, S. Krishnan, and E. J. Hopfinger, “Shape parameters of a floating bubble,” Phys. Fluids 30, 112105 (2018). 





52N. Van-Tu and W. Park, “A volume-of-fluid (VOF) interface-sharpening method for two-phase incompressible flows,” Comput Fluids 152, 104 (2017). 





15J. Hua, J. F. Stene, and P. Lin, “Numerical simulation of 3D bubbles rising in viscous liquids using a front tracking method,” J. Comput. Phys. 227, 3358 (2008). 





53C. Beckermann, H. J. Diepers, I. Steinbach, A. Karma, and X. Tong, “Modeling melt convection in phase-field simulations of solidification,” J. Comput. Phys. 154, 468 (1999). 





16A. Zhang, P. Sun, and F. Ming, “An SPH modeling of bubble rising and coalescing in three dimensions,” Comp. Methods Appl. Mech. Eng. 294, 189 (2015). 





54Z. Guo, C. Zheng, and B. Shi, “Discrete lattice effects on the forcing term in the lattice Boltzmann method,” Phys. Rev. E 65, 046308 (2002). 





17X. Niu, Y. Li, Y. Ma, M. Chen, X. Li, and Q. Li, “A mass-conserving multiphase lattice Boltzmann model for simulation of multiphase flows,” Phys. Fluids 30, 013302 (2018). 





55M. Cheng, J. Hua, and J. Lou, “Simulation of bubble-bubble interaction using a lattice Boltzmann method,” Comput. Fluids 39, 260 (2010). 





18M. Alizadeh, S. M. Seyyedi, M. T. Rahni, and D. D. Ganji, “Three-dimensional numerical simulation of rising bubbles in the presence of cylindrical obstacles, using lattice Boltzmann method,” J. Mol. Liq. 236, 151 (2017). 





56J. Du, A. Zhang, Z. Guo, M. Yang, M. Li, and S. Xiong, “Mechanism of the growth pattern formation and three-dimensional morphological transition of hcp magnesium alloy dendrite,” Phys. Rev. Mater. 2, 083402 (2018). 





19A. Zhang, Z. Guo, and S. M. Xiong, “Quantitative phase-field lattice-Boltzmann study of lamellar eutectic growth under natural convection,” Phys. Rev. E 97, 053302 (2018). 





57R. Krishna, M. I. Urseanu, J. M. van Baten, and J. Ellenberger, “Wall effects on the rise of single gas bubbles in liquids,” Int. Commun. Heat Mass Transfer 26, 781 (1999). 





20A. Zhang, J. Du, Z. Guo, Q. Wang, and S. Xiong, “A phase-field lattice-Boltzmann study on dendritic growth of Al-Cu alloy under convection,” Metall. Mater. Trans. B 49, 3603 (2018). 





58L. Amaya-Bower and T. Lee, “Single bubble rising dynamics for moderate Reynolds number using Lattice Boltzmann Method,” Comput. Fluids 39, 1191 (2010). 





21A. Zhang, J. Du, Z. Guo, Q. Wang, and S. Xiong, “Dependence of lamellar eutectic growth with convection on boundary conditions and geometric confinement: A phase-field lattice-Boltzmann study,” Metall. Mater. Trans. B 50, 517 (2019). 





59J. Hua and J. Lou, “Numerical simulation of bubble rising in viscous liquid,” J. Comput. Phys. 222, 769 (2007). 





22A. Zhang, J. Du, Z. Guo, Q. Wang, and S. Xiong, “Abnormal solute distribution near the eutectic triple point,” Scr. Mater. 165, 64 (2019). 





60D. D. Joseph, “Rise velocity of a spherical cap bubble,” J. Fluid Mech. 488, 213 (2003). 





23J. Du, A. Zhang, Z. Guo, M. Yang, M. Li, F. Liu, and S. Xiong, “Atomistic underpinnings for growth direction and pattern formation of hcp magnesium alloy dendrite,” Acta Mater. 161, 35 (2018). 





61J. P. Best, “The formation of toroidal bubbles upon the collapse of transient cavities,” J. Fluid Mech. 251, 79 (1993). 





24X. Fan, A. Zhang, Z. Guo, X. Wang, J. Yang, and J. Zou, “Growth behavior of γ′ phase in a powder metallurgy nickel-based superalloy under interrupted cooling process,” J. Mater. Sci. 54, 2680 (2019). 





62E. A. Brujan, G. S. Keen, A. Vogel, and J. R. Blake, “The final stage of the collapse of a cavitation bubble close to a rigid boundary,” Phys. Fluids 14, 85 (2002). 





25A. Zhang, S. Meng, Z. Guo, J. Du, Q. Wang, and S. Xiong, “Dendritic growth under natural and forced convection in Al-Cu alloys: From equiaxed to columnar dendrites and from 2D to 3D phase-field simulations,” Metall. Mater. Trans. B 50, 1514 (2019). 





63G. Brereton and D. Korotney, Coaxial and oblique coalescence of two rising bubbles. Dynamics of bubbles and vortices near a free surface (ASME, New York, 1991), Vol. AMD 119, p. 50. 





26P. Chiu and Y. Lin, “A conservative phase field method for solving incompressible two-phase flows,” J. Comput. Phys. 230, 185 (2011). 





64N. Takada, M. Misawa, A. Tomiyama, and S. Hosokawa, “Simulation of bubble motion under gravity by lattice Boltzmann method,” J. Nucl. Sci. Technol. 38, 330 (2001). 





27Y. Q. Zu and S. He, “Phase-field-based lattice Boltzmann model for incompressible binary fluid systems with density and viscosity contrasts,” Phys. Rev. E 87, 043301 (2013). 





65C. W. Stewart, “Bubble interaction in low-viscosity liquids,” Int. J. Multiphase Flow 21, 1037 (1995). 





66C. Brücker, “Structure and dynamics of the wake of bubbles and its relevance for bubble interaction,” Phys. Fluids 11, 1781 (1999). 





67O. Shardt, S. K. Mitra, and J. J. Derksen, “The critical conditions for coalescence in phase field simulations of colliding droplets in shear,” Langmuir 30, 14416 (2014). 





68W. F. Smith and J. Hashemi, Foundations of Materials Science and Engineering (McGraw-Hill Education-Europe, USA, 2009). 

