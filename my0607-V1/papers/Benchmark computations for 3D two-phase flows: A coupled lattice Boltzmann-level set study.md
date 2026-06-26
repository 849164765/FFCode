# Benchmark computations for 3D two-phase flows: A coupled lattice Boltzmann-level set study

Mohammad Amin Safia,b*, Nikolaos Prasianakisc, Stefan Turekb 

*Corresponding author; seyed.safi@psi.ch 

aPaul Scherrer Institute, Energy and Environment Division, Combustion Research Laboratory, CH-5232 Villigen-PSI, Switzerland 

bInstitute of Applied Mathematics (LSIII), TU Dortmund, Vogelpothsweg 87, D-44221 Dortmund, Germany 

cPaul Scherrer Institute, Nuclear Energy and Safety Division, Waste Management Laboratory, CH-5232 Villigen-PSI, Switzerland 

## Abstract

Following our previous work on the application of the diffuse interface coupled lattice Boltzmann-level set (LB-LS) approach to benchmark computations for 2D rising bubble simulations, this paper investigates the performance of the coupled scheme in 3D two-phase flows. In particular, the use of different lattice stencils, e. g., D3Q15, D3Q19 and D3Q27 is studied and the results for 3D rising bubble simulations are compared with regards to isotropy and accuracy against those obtained by finite element and finite difference solutions of the Navier-Stokes equations. It is shown that the method can eventually recover the benchmark solutions, provided that the interface region is aptly refined by the underlying lattice. Following the benchmark simulations, the application of the method in solving other numerically subtle problems, e. g., binary droplet collision and droplet splashing on wet surface under high Re and We numbers is presented. Moreover, implementations on general purpose GPUs are pursued, where the computations are adaptively refined around the critical parts of the flow. 

Keywords 3D rising bubble; Lattice Boltzmann method; Level set method; Droplet splashing; Binary droplet collision; GPGPU implementation 

## 1 Introduction

The application of the coupled lattice Boltzmann-level set scheme for two-phase flows to 2D benchmark simulations was discussed previously in the works of Safi and Turek in [36], coupling a one-fluid lattice Bolztmann equation (LBE) with signed distance level set equation (LSE) as a sharp interface method. This methodology was shown to be prone to unphysical velocities caused by weak pressure approximation in the one-fluid LBE. The alternative pressure evolution LBE is more robust at large pressure gradients but breaks down in calculating non-ideal part of the pressure if a sharp interface method is used. In the later work of Safi and Turek [35] the pressure evolution LBE was thus coupled with a phase field mass conserving LSE. It was shown in [35] that by applying averaged directional differencing to discretize the pressure forcing terms in LBE, one can use the continuum surface force (CSE) form of the surface tension forces and obtain stable and accurate solutions at high density and viscosity ratios as examined for 2D rising bubble benchmarks. In fact, such a methodology also offers an alternative to the two-LBE schemes as in the pressure evolution scheme of Lee et. al [24] or the free-energy implementation of Banari et. al. [6] which require to save and process data for a second group of distribution functions which are used to solve the LBE for the order parameter. Considering the fact that diffuse interface LB-based solutions often require high lattice resolutions and 3D lattice stencils need 15-27 populations, pursuing a two-LBE strategy for 3D flows ends up in extremely large amounts of computational time and memory. The memory issue is even more demanding for parallel implementations on general purpose GPUs (GPGPU) where memory resources per GPGPU are scarce. This in turn puts severe restrictions on the resolution as well as the size of the problem and could potentially limit the applicability of such schemes. 

On the other hand, to the knowledge of the authors, 3D simulations of two-phase flows have always been validated against experimental results in the so-called picture norms as in [3, 8, 25], where correspondence between the input setting of the numerical simulation to the setup data in the laboratory is often very hard, if not impossible. Consequently, from the users point of view, it is yet unclear how close the obtained numerical results could get to the direct numerical simulation of the two-phase Navier-Stokes equations (NSE), especially in situations where the exact extent and rate of interface deformation are of paramount importance. 

In this paper, we extend the methodology of diffuse coupled LB-LS to 3D flows, where the recent results obtained by finite element solutions of NSE for 3D rising bubbles [43] will be used as reference to examine the accuracy of the coupled scheme for different test cases. The benchmarking also enables us to investigate the effect of using various discrete velocity stencils in LB method. Having established the accuracy in 3D, the numerical robustness of the coupled LB-LS scheme will be examined for non-benchmark flows involving strong deformations, e. g., droplet splashing on liquid film and binary droplet collision under high Re and We numbers. On the computational side, some adaptive computation techniques for parallel GPGPU implementations are discussed so as to ease computational hurdles of two-phase 3D simulations. 

The paper is organized as follows. The coupled LB-LS scheme will be shortly explained in section 2. Benchmark test cases for 3D rising bubbles will be introduced in section 3, including a brief review of the numerical procedure used in the reference computations. The section eventually presents the numerical results for the rising bubble problems. Application of the scheme to other numerically and physically interesting two-phase flows will be the subject of section 5. Section 6 briefly describes the GPGPU implementation and the obtained computational performances. The paper will be closed with conclusions and discussions in section 7. 

## 2 Coupled LB-LS scheme

The LBE for the evolution of density distribution functions $f _ { \alpha }$ along discrete velocity directions $\alpha = 0 , . . . , n$ with BGK approximation for collision reads [15, 17] 

$$
\frac {\partial f}{\partial t} + c _ {\alpha} \cdot \nabla f = - \boldsymbol {\Lambda} (f _ {\alpha} (x, t) - f _ {\alpha} ^ {e q} (x, t)) - \frac {(c _ {\alpha} - \mathbf {u}) \cdot \mathbf {F}}{c _ {s} ^ {2}} f ^ {e q} (x, t) \tag {1}
$$

where F is the forcing term and $f _ { \alpha } ^ { e q }$ is the equilibrium distribution 

$$
f _ {\alpha} ^ {e q} (\rho , \mathbf {u}) = w _ {\alpha} \rho \left[ 1 + \frac {c _ {\alpha} \cdot \mathbf {u}}{c _ {s} ^ {2}} + \frac {[ c _ {\alpha} \cdot \mathbf {u} ] ^ {2}}{2 c _ {s} ^ {4}} + \frac {[ \mathbf {u} \cdot \mathbf {u} ]}{2 c _ {s} ^ {2}} \right] (2)
$$

with $w _ { \alpha }$ being the weight factors for each discrete direction. The general relaxation matrix Λ could be replaced by $1 / \tau$ for a single relaxation time (SRT) collision, while in a multiple relaxation time (MRT) collision it takes the following form [22] 

$$
\boldsymbol {\Lambda} = \mathbf {M} ^ {- 1} \hat {\boldsymbol {\Lambda}} \mathbf {M} - 2 \mathbf {I} \tag {3}
$$

where M constructs n hydrodynamic moments from $f _ { \alpha } \mathrm { s }$ and $\hat { \bf \cal N } = D i a g \{ s _ { 1 } , . . . , s _ { n } \}$ is the diagonal relaxation matrix so as to let different hydrodynamic moments relax to their corresponding equilibrium states via individual relaxation rates. The force term F arising from two-phase effects in the LB framework is 

$$
\mathbf {F} = \nabla (\rho c _ {s} ^ {2} - p) - \mathbf {F} _ {s} + \mathbf {G}. \tag {4}
$$

The first term in equation 4 is the gradient of the non-ideal part of the pressure √ $p ,$ where $\rho$ is the density and $c _ { s } = 1 / \sqrt { 3 }$ is the lattice speed of sound. G is a volume force, e. g. gravity. The second term, is the diffuse CSF form of surface tension force [32] 

$$
\mathbf {F} _ {s} = \sigma \kappa \mathbf {n} \delta_ {\varepsilon} \tag {5}
$$

where $\sigma$ is the surface tension coefficient, n is the normal to the interface, κ is the interface curvature. The delta function $\delta _ { \varepsilon }$ is used to apply the force in a diffuse manner over a thickness ε around the interface. In addition, one may define the continuous phase-field scalar ψ to denote the interface location as [31, 32] 

$$
\psi = \left\{ \begin{array}{l l} 1 & i f \quad \mathbf {x} \in \Omega_ {1} \\ 0. 5 & i f \quad \mathbf {x} \in \Gamma \\ 0 & i f \quad \mathbf {x} \in \Omega_ {2} \end{array} \right. \tag {6}
$$

where $\varOmega _ { 1 }$ and $\varOmega _ { 2 }$ point to phase 1 and 2, respectively and Γ is the location of the interface. Note that ψ decays continuously from 1 to 0 while crossing the interface. Moreover, based on the phase-field definition of ψ one may deduce $\nabla \rho = \nabla \psi ( \rho _ { 1 } - \rho _ { 2 } )$ . The definition also helps constructing interface properties as functions of $\psi$ 

$$
\mathbf {n} (\psi) = \frac {\nabla \psi}{| \nabla \psi |} \qquad , \qquad \kappa (\psi) = \nabla \cdot \mathbf {n} = \nabla \cdot \left(\frac {\nabla \psi}{| \nabla \psi |}\right) \qquad , \qquad \delta_ {\varepsilon} (\psi) = | \nabla \psi |. (7)
$$

In order to rectify the destructive numerical effect of the $\nabla ( p - \rho c _ { s } ^ { 2 } )$ at high density ratios, He et. al [17] proposed to use pressure distribution functions $g _ { \alpha }$ 

$$
g _ {\alpha} = f _ {\alpha} c _ {s} ^ {2} + (p - \rho c _ {s} ^ {2}) w _ {\alpha}, \tag {8}
$$

and thus solve the pressure evolution LBE by taking the total derivative of $g _ { \alpha }$ as 

$$
\frac {D g _ {\alpha}}{D t} = \frac {\partial g _ {\alpha}}{\partial t} + c _ {\alpha}. \nabla g _ {\alpha} = \boldsymbol {\Lambda} (g _ {\alpha} - g _ {\alpha} ^ {e q}) + (c _ {\alpha} - \mathbf {u}) \cdot [ \nabla \varphi (\Gamma_ {\alpha} (\mathbf {u}) - w _ {\alpha}) + (\mathbf {F} _ {s} + \mathbf {G}) \Gamma_ {\alpha} (\mathbf {u}) ]. \tag {9}
$$

Details for time integration and spatial discretization of equation (9) could be found in [24, 35]. In essence, by applying trapezoidal rule for time integration and using average and central directional differencings, one obtains 

$$
\bar {g} _ {k} (x + c _ {\alpha} \Delta t, t + \Delta t) = \bar {g} _ {k} (x, t) - (\Lambda + 2 \mathbf {I}) (\bar {g} _ {\alpha} (x, t) - \bar {g} _ {\alpha} ^ {e q} (x, t)) +
$$

$$
\left. \left(c _ {\alpha} - \mathbf {u}\right) \cdot \left[ \left(\sigma \nabla^ {c} \cdot \left(\frac {\nabla^ {c} \psi}{| \nabla^ {c} \psi |}\right) \nabla^ {c} \psi + \mathbf {G}\right) \Gamma_ {\alpha} (\mathbf {u}) \right] \right| _ {(x, t)} + \tag {10}
$$

$$
\left. \left(c _ {\alpha} - \mathbf {u}\right) \cdot \left[ \Delta \rho c _ {s} ^ {2} \left(\Gamma_ {\alpha} (\mathbf {u}) - w _ {\alpha}\right) \nabla^ {a v e} \psi \right] \right| _ {(x, t)}
$$

where the transformed distributions $\bar { g } _ { \alpha }$ and $\bar { g } _ { \alpha } ^ { e q }$ are defined as 

$$
\bar {g} _ {\alpha} = g _ {\alpha} + \frac {\Lambda}{2} \left(\bar {g} _ {\alpha} - \bar {g} _ {\alpha} ^ {e q}\right) - \frac {1}{2} \left(c _ {\alpha} - \mathbf {u}\right) \cdot \left[ \nabla \varphi \left(\Gamma_ {\alpha} (\mathbf {u}) - w _ {\alpha}\right) + \left(\mathbf {F} _ {s} + \mathbf {G}\right) \Gamma_ {\alpha} (\mathbf {u}) \right] \tag {11}
$$

$$
\bar {g} _ {\alpha} ^ {e q} = g _ {\alpha} ^ {e q} - \frac {1}{2} (c _ {\alpha} - \mathbf {u}) \cdot [ \nabla \varphi (\Gamma_ {\alpha} (\mathbf {u}) - w _ {\alpha}) + (\mathbf {F} _ {s} + \mathbf {G}) \Gamma_ {\alpha} (\mathbf {u}) ].
$$

Finally, the pressure and velocity of the flow are recovered as moments of new distribution functions $\bar { g } _ { \alpha }$ 

$$
p = \sum_ {\alpha} \bar {g} _ {\alpha} + \frac {c _ {s} ^ {2} \Delta \rho}{2} \mathbf {u} \cdot \nabla^ {c} \psi \tag {12}
$$

$$
\rho \mathbf {u} = \frac {1}{c _ {s} ^ {2}} \sum_ {\alpha} \bar {g} _ {\alpha} + \frac {1}{2} (\sigma \nabla^ {c} \cdot (\frac {\nabla^ {c} \psi}{| \nabla^ {c} \psi |}) \nabla^ {c} \psi + \mathbf {G}). \tag {13}
$$

In order to capture the interface, the LBE in equation (10) is coupled with a mass conserving LSE proposed by Olsson and Kreiss [30] so as to solve the evolution of the scalar field ψ as it moves with the flow velocity u 

$$
\partial_ {t} \psi + \mathbf {u} \cdot \nabla \psi = - \nabla \cdot (\psi (1 - \psi) \mathbf {n}) + \eta \nabla^ {2} \psi . \tag {14}
$$

Note that the initial profile of ψ is obtained by iteratively solving equation (14) to reach steady state. The parameter η controls the thickness and could be chosen as 

$$
\eta = \frac {h ^ {\beta}}{2} \tag {15}
$$

where h is the spatial resolution and the exponent β is chosen to be close to 1 so as to control the interface thickness and has a subsequent effect on the quality of the interface capturing as well as the overall mass conservation. We use second order Runge-Kutta for time integration of equation (14) and the fifth order WENO for its convective term. More elaborate discussions on the use of lower order schemes, e. g, ENO are given [26] and [31]. For the present simulation it is observed that using ENO leads to slightly poorer results for high Eo number flows where strong deformations create kinks and singularities which are hard to be resolved using lower order schemes. It must be noted that the computational overhead associated with using the WENO scheme would be only moderately more than that of ENO when implemented in an adaptive way on GPGPUs. This is discussed more in section 6. Finally, the terms on the right hand side are then treated explicitly and discretized in space using central differencing. Since the time scales enforced by the solution of the LBE are quite small, such a discretization for these non-linear terms will not introduce any noticeable numerical deficiency. 

## 2.1 MRT collision for 3D LBM

In order to reinforce the stability of the LBE solver at low viscosities and high Re numbers, the present implementation employs a MRT collision for 3D implementations through replacing Λ + 2I with $\mathbf { M } ^ { - 1 } \hat { \mathbf { A } } \mathbf { M }$ in equation (10). An extensive description of the principal moments for D3Q15 and D3Q19 stencils is done in [11, 28]. For the D3Q27 stencil, the MRT methodology given in [39] for turbulent flows in adopted here. The general choice of the relaxation times for all stencils would be analogous to those outlined in our previous 2D implementations [35, 36] which requires under-relaxing the energy moments so as to eliminate the initial fluctuations in the velocity of a rising bubble. A detailed discussion is beyond the scope of this paper and could be found in [34]. 

## 3 3D rising bubble benchmarks

Inspired by the 2D work in [18], the rising bubble benchmarks were extended for 3D effects in the works of Turek et. al [43] and Adelsberger et al. [2]. The problem consists of an initially stagnant bubble of radius $r _ { 0 } = 0 . 2 5$ and located at $( x , y , z ) = ( 0 . 5 , 0 . 5 , 0 . 5 )$ which starts to rise due to buoyancy through the liquid in a 3D rectangular domain of size $1 \times 1 \times 2$ as depicted in figure 1. The buoyancy force is 

$$
\mathbf {G} = \left(\rho_ {1} - \rho_ {2}\right) \mathbf {g} \tag {16}
$$

where $\mathbf { g } = ( 0 , - g )$ is the gravity and $\rho _ { 1 }$ and $\rho _ { 2 }$ are the liquid and gas densities, respectively. Note that G is only exerted as a net force to the gas inside the bubble. No-slip boundary conditions are applied to all boundaries including the side walls where the second order half-way bounceback rule is imposed. Moreover, due to the symmetry of the problem in the horizontal plane, only $1 / 4$ of the problem is solved to save computational cost and memory, where the center of the bubble will be located at $( x , y , z ) = ( 0 , 0 , 0 . 5 )$ and the symmetry boundary condition is applied to $x = 0$ and $y = 0$ planes. 

Two sets of input data are then used to define two test cases as summarized in table 1. The hydrodynamics in both cases are governed by the non-dimensional Re and Eo numbers 

$$
R e = \frac {\rho_ {1} \sqrt {g} (2 r _ {0}) ^ {3 / 2}}{\mu_ {1}} \tag {17}
$$

$$
E o = \frac {4 \rho_ {1} g (r _ {0}) ^ {2}}{\sigma}. \tag {18}
$$

In principle, the Eo number describes the ratio of momentum forces over the surface tension forces in a nondimensional fashion. The LB code uses the same values of density and viscosity as in table 1, while $g _ { l b }$ and $\sigma _ { l b }$ are obtained based on equations 17 and 18 using LB units for $r _ { 0 } .$ . In order to preserve consistency with the macroscopic time measurement in [43] and [2], the following equation is used to calculate the macroscopic time based on LB iterations $t _ { l b }$ 

$$
T = t _ {l b} \sqrt {\frac {g _ {l b}}{g L _ {0}}}
$$

where $L _ { 0 } = 1 / h$ is the characteristic length in lattice units, with h being the physical resolution. For both test cases, simulations are carried until $T = 3$ . Analogous to the 2D problems in [35], the first test case, referred to as TC1, is a low density ratio system where the the relatively small value of Eo number tends to keep the bubble in a round shape as it rises. In the second test case, however, the density and viscosity ratios are increased to 1000 and 100, respectively, and the Eo number is increased to 125. The latter particularly allows for large deformation rates implied by weak surface tension forces. It must be noted that in the preset simulations the maximum Mach number $M a = m a x ( \mid u \mid ) / c _ { s }$ , does not exceed 0.05 and corresponds to the velocity inside the gas phase at the coarsest lattice of $1 / h = 8 0$ . This falls well inside the incompressibility limit of $M a < 0 . 1 5$ in LB method [16]. 

![](images/b997b5d4bbced1da0067a43c7493605d5ca81ea5bd57785643214b20342a7957.jpg)



Fig. 1: Illustration of the rising bubble problems in the 3D problem.



Table 1: Physical parameters and dimensionless numbers for TC1 and TC2 rising bubble problems.


<table><tr><td>Test case</td><td><eq>\rho_1</eq></td><td><eq>\rho_2</eq></td><td><eq>\mu_1</eq></td><td><eq>\mu_2</eq></td><td>g</td><td>σ</td><td>Re</td><td>Eo</td><td><eq>\rho_1/\rho_2</eq></td><td><eq>\mu_1/\mu_2</eq></td></tr><tr><td>TC1</td><td>1000</td><td>100</td><td>10</td><td>1</td><td>0.98</td><td>24.5</td><td>35</td><td>10</td><td>10</td><td>10</td></tr><tr><td>TC2</td><td>1000</td><td>1</td><td>10</td><td>0.1</td><td>0.98</td><td>1.96</td><td>35</td><td>125</td><td>1000</td><td>100</td></tr></table>

In order to quantify the dynamics of the bubble during its course of rise, a number of benchmark quantities are defined and used throughout this paper. To track the bubble position during the rise process, the bubble centroid coordinate $( x _ { c } , y _ { c } )$ is computed as 

$$
X _ {c} = (x _ {c}, y _ {c}) = \frac {\int_ {\Omega_ {2}} \mathbf {x} d x}{\int_ {\Omega_ {2}} 1 d x}
$$

where $\varOmega _ { 2 }$ encompasses all the lattice points inside the bubble. The degree of roundness of a 3D bubble can be measured based on the sphericity parameter defined in [44] as 

$$
\Psi = \frac {A _ {a}}{A _ {b}} = \frac {\mathrm{areaofvolume-equivalentsphere}}{\mathrm{areaofbubble}} = \frac {\pi^ {1 / 3} (6 V _ {b}) ^ {2 / 3}}{A _ {b}}.
$$

The mean vertical velocity with which the bubble is rising is defined as 

$$
U _ {c} = \frac {\int_ {\Omega_ {2}} u _ {z} d x}{\int_ {\Omega_ {2}} 1 d x}.
$$

where $u _ { z }$ is the velocity component in the z direction. The last quantity is the bubble size which describes the maximum extension (diameter) of the bubble in the main coordinate directions which together with the bubble circularity and sphericity further specifies the deformation state of the bubble. It is calculated as 

$$
d _ {i} = m a x _ {p, q \in \Omega_ {2}} | p _ {i} - q _ {i} |, \qquad i = x, y, z.
$$

where $p$ and $q$ are any two arbitrary points which belong to the bubble. 

## 3.1 Review of the numerical tools

The numerical simulations performed by the current coupled LB-LS scheme, will be compared against two sets of solutions of the two-phase flows. The first group, also considered here as the reference solution, is the one obtained by the finite element FeatFlow package [1] designed for solving incompressible NSE. The two-phase solver is based on coupling the NSE with the signed distance-based LSE. Space discretization in the 3D code makes use of $Q _ { 2 }$ elements for velocity and discontinuous $P _ { 1 }$ elements for pressure. Besides $Q _ { 2 }$ elements are used for solving the LS equation. Time integration is carried out via a Crank-Nicolson scheme which is used in the 3D implementations. A full description of the solvers could be found in [43]. 

Another two-phase flow solver for 3D simulations is the finite difference-based NaSt3D code [2]. It solves the Navier-Stokes equation on an equidistant grid. Chorin’s projection method is used to decouple velocity and pressure fields, where a second order explicit Adams-Bashforth scheme is employed to solve the velocity while the pressure is recovered via solving the Poisson equation. The interface is captured through a signed distance LS function and the bubble mass is corrected using the local correction scheme of Sussman and Fatemi [40]. 

## 4 Numerical Results

The D3Q19 stencil is adopted as the default discrete velocity stencil which is known to posses a higher isotropy than the D3Q15 stencil and imposes only a moderate increase in the computational workload. The choice will be verified later on through detailed comparisons of the results obtained by D3Q15, D3Q19 and D3Q27 velocity models against the benchmark data. Yet, regardless of the type of the lattice stencil in use and considering the 5.7 GB memory available on Kepler 20x GPGPU used for the present simulations, the maximum lattice resolution is selected to be $1 / h = 2 5 6$ corresponding to a grid of $1 2 8 \times 1 2 8 \times 5 1 2$ which occupies 3.9 GB of memory for $1 / 4$ of the problem using the D3Q19 stencil. 

## 4.1 Test case 1

The evolution of the 3D bubble interface in TC1 up to $T = 3$ is presented in figure 2 on a lattice of $1 / h = 2 5 6$ and the reference shapes obtained by FeatFlow on a grid of $1 / h = 1 2 8$ . The bubble develops into an ellipsoidal shape as predicted in the bubble shape regime map of Clift et. al [9]. However, compared to the 2D cases in [35], it preserves its initial sphericity to a higher degree such that almost no deformation occurs after $T = 2 . 0$ . 

The convergence of the interface shape is also studied in figure 3, followed by more detailed comparisons between the bubble quantities obtained by LBM NaSt3D, and the reference values of FeatFlow in figures 4 to 7. For the low deformation rate of TC1, the present maximum lattice resolution seems to be sufficient to repeat the reference data by FeatFlow. For the rise velocity, the discrepancies grow around the peak values although the difference between the LBM data on $1 / h = 2 5 6$ and that of FeatFlow amounts to less than 0.5%. A similar behavior is observed for the interface-related sphericity and diameter values around $T = 1 . 5$ where the overall convergence trend is towards the reference curves. 

![](images/0c58cb4f40672530a8bc6e8fbe11138b0a4393d03a43e3b184696e423d928758.jpg)



(b)



Fig. 2: Time evolution of the bubble in 3D problem of TC1 using (a) coupled LB-levelset, (b) FeatFlow.


![](images/f4d3ff13366e00b7c4bd5265684fa06f2a0b5b5f614412c4742b2259a65b569e.jpg)



Fig. 3: Terminal shape of bubble at $T = 3$ for 3D problem of TC1 with enlarged view (right), on different lattice resolutions; $1 / h = 9 6$ (black), $1 / h = 1 9 2$ (blue), $1 / h = 2 5 6$ (green) and FeatFlow on $1 / h = 1 2 8 { \mathrm { ~ ( r e d ) } }$ .


![](images/3fb3ab7475cc52c4891e1357313977573682538302d45839bbdf5e1de9673701.jpg)


![](images/ce19b319e992254797b58eef554fd64a1a92ffdad474c8e74f7eeee8d1d6727f.jpg)



Fig. 4: Time evolution of center of mass for 3D rising bubble TC1 (left), and enlarged view (right).


![](images/e528c28384bcf80f0a098e7af531a31e66de3b9c30d14e42873ceedb20398434.jpg)


![](images/77660527f44e1bacd27764bb95dbde00a3794e89616768e904799840485924d4.jpg)



Fig. 5: Time evolution of rise velocity for 3D rising bubble TC1 (left), and enlarged view (right).


![](images/004a799087031d6534a9cdda07f3db685ef4f44c8af5c62fc4e878178189d106.jpg)


![](images/54c55858aa60432d9194c6a24cb0757aaacd5c79be313d79afce13a3bc2d852f.jpg)



Fig. 6: Time evolution of circularity for 3D rising bubble TC1 (left), and enlarged view (right).


![](images/1655a84b0300052477c340955b5faebef54af1d117cbb47c3409fde305f3006b.jpg)


![](images/8dab6115b12a026119001b79802cad8245f0dcf54cd9011698e614f4b9fbad9f.jpg)



Fig. 7: Time evolution of diameters in horizontal and vertical directions for 3D rising bubble TC1 (left), and enlarged view (right).


## 4.2 Test case 2

The 3D bubble shown in figure 8(a) experiences strong deformations provoked by the weaker surface tension effects due to a relatively large Eo number. The bubble first deforms into a dimpled cap and then proceeds to extend the edges so that the eventual cusp-shaped bubble is formed as predicted in [9]. Yet again, the diffuse nature of the present coupled scheme prevents the bubble to develop sharp edges, e. g., those predicted by Featflow in figure 8(b). Nevertheless, the picture norm admits the close relevance of the overall evolution trend in both schemes. Examining the convergence trend in the interface shape in figure 9 shows that LBM may eventually converge to the sharp interface solution of the two-phase Navier-Stokes equations, provided that the lattice would be sufficiently refined, as in the 2D counterpart of TC2 in [35] where a resolution of $1 / h = 6 4 0$ was used to approach the sharp interface results. 

The convergence of the bubble quantities together with the NaSt3D data and the reference data of FeatFlow are collected in figures 10 to 13. The agreement between the rise velocities keeps favorably up to $T = 1 . 5$ which confirms a decent performance of the D3Q19 discrete model in recovering the correct pressure and velocity fields under high density and viscosity differences. For $T > 1 . 5 ,$ as the bubble starts to develop the cusp shape, the phase-field LSM falls behind the sharp LSM, and tries to keep the edges prolonged and smooth. In terms of bubble quantities, this results in growing deviations from the sharp interface data, which is more clearly seen in the bubble sphericity, diameter and centroid position. Similar to the 2D test case, the sharp interface methods tend to create kinks or teeth-like irregularities on the skirt of the bubble towards T = 2 which then makes the solution strongly dependent on the specific parametrization used in different codes. 

![](images/58b5beaee2e2994caf0073cb714495846da0b9b218a8585d56169dbf4662ffce.jpg)



(b)



Fig. 8: Time evolution of the bubble in 3D problem of TC2 using (a) coupled LB-levelset (b) FeatFlow. The dark blue regions reflect the interior surface.


![](images/da5bb4e402605ac391497265556b0e2c7cc38ff12d26a26ae87d9cd95a72e2fa.jpg)



Fig. 9: Terminal shape of bubble at $T = 3$ for 3D problem of TC2 with enlarged view (right), on different lattice resolutions; $1 / h = 9 6 \ \mathrm { ( b l a c k ) }$ , $1 / h = 1 9 2$ (blue), $1 / h = 2 5 6$ (green) and FeatFlow on $1 / h = 1 2 8 { \mathrm { ~ ( r e d ) } }$ .


![](images/2f031c0458019d48da3000cd44c78f2bd874ba4c479cd0d9d66c9abde9dcbc13.jpg)


![](images/1ff776745b95bc9d9448a7efb0e98c689e3ff19ba8578813f6cfe18a6faa3b59.jpg)



Fig. 10: Time evolution of center of mass for the 3D rising bubble TC2 (left), and enlarged view (right).


![](images/d3468ed86eb30bd6fb50a9a0e6aa3e59e08b422b82f0da1d11201b6895976e71.jpg)


![](images/fa3d7ebd903947c160af9e7b886a782fb3729e9653773bc5fe23960c3adfe6b6.jpg)



Fig. 11: Time evolution of rise velocity for the 3D rising bubble TC2 (left), and enlarged views (right).


![](images/d96e1ef60d8946f21a8eb9bdbb67d98959c87e89627112fdcc504322d445ba90.jpg)


![](images/e2475a83248f270949796f5db6de48ea9753a70bcbee4530957041fab138f562.jpg)



Fig. 12: Time evolution of circularity for the 3D rising bubble TC2 (left), and enlarged view (right).


![](images/9c26939d258b68477718a891b3b9ecc0e824e891435d5f50821a001564cadd38.jpg)


![](images/cb8cb2c82e229ee1bcd215e77fbbbfc0d857bcf4ad2d5200b19c85043fa14cda.jpg)



Fig. 13: Time evolution of diameters in horizontal and vertical directions for 3D rising bubble TC2 (left), and enlarged view (right).


## 4.3 Isotropy and the choice of lattice stencil

The most prominent distinction between discrete velocity stencils, is their associated degrees of isotropy. This property could be best studied through examining the performance of different stencils for TC2, where the presence of four orthogonal walls and the large deformations of the bubble creates a natural anisotropy in the terminal shape of the bubble as shown in figure 15 for the front view and the view from the middle section. The interface lines for the plane perpendicular to the wall $( \theta = 0 )$ and diagonal plane (θ = 45) are also depicted in figure 16. In the first look, the comparisons unveil the overestimated deformation produced by the D3Q15 model on the diagonal plane. The deviation is caused by the inherent lack of isotropy of the D3Q15 stencil in the diagonal directions on $x - y , x - z , y - z$ planes (see figure 14). Although one expects to gain the best relevance to the FeatFlow results using the D3Q27 stencil, it is the D3Q19 stencil which exhibits a closer similarity to the FeatFlow shape. 

Unlike the observations in the picture norm, the bubble quantities well coincide over time and one could hardly notice significant differences as seen in figures 17 and 18 for bubble velocity and diameter, respectively, on a $1 / h = 2 5 6$ lattice. For the bubble diameter which is more sensitive to minimal changes in the interface, the D3Q15 stencil shows a more steep return in the vertical diameter up to nearly $d _ { z } = 0 . 5$ which is also seen by the interface lines in figure 16 where the strong anisotropy yields a longer skirt of the bubble in the diagonal plane. The rather obviously poor performance of D3Q15 again could largely be related to its lower degree of isotropy compared to D3Q19 and D2Q27 stencils. Similar effects are seen even for D3Q19 stencil at very high Re numbers as in the work of Geller et al. [13]. Here, however, the slightly superior accuracy and isotropy of the D3Q19 stencil over the D3Q27 discrete model is questionable. A likely reason could lie in the choice of the relaxation times in the MRT scheme for different stencils. Exhaustive analysis regarding the optimum choice for the values of the relaxation times is presented in [34]. Therein, the effect of underrelaxation of specific moments has been highlighted. Extending the stencil from 19 to 27 discrete velocities inevitably extends also the range of the hydrodynamic moments that are under consideration. In fact this adds 8 more moments and their respective relaxation times that have to be chosen carefully. The relaxation times of the D3Q27 stencil have been regulated based on the analysis of the D3Q19 stencil and it cannot be excluded that a better set of relaxation parameters could smear out the differences between the two stencils. On the contrary, such an uncertainty does not exist in the case of the D3Q15 stencil, since it is a subset of the D3Q19, and the optimum relaxation times remain exactly the same (for the same moments). Nevertheless, the above argument could not be firmly verified as almost all dynamic two-phase LB simulations in the literature are limited to D3Q15 or D3Q19 stencils for computational reasons [7, 41, 46, 47]. Lee and Liu [24] used the D3Q27 discrete model to simulate droplet impact on dry surfaces. Yet they brought no justification for this choice and provided no comparisons against other discrete models. 

![](images/d5fbe761c5468b2c4814bb1a06486b321d0c58ad1742f2596d916b4b6eceb1d5.jpg)



(a) D3Q15


![](images/204dea499554fc68f742b9cf431ad171de0b97d43f6cb8f7e9434b58e28ccc3d.jpg)



(b) D3Q19


![](images/80b913472059b8e4dc813d79c2422d1f8eb8de74c227d8c7c008c8a2e655b1b5.jpg)



(c) D3Q27



Fig. 14: Different lattice stencils for three dimensional LBM



D3Q15


![](images/39d12560224d8ea39faa542d6f040b7bd4ca7d08579c6c0145411ab48015dd1e.jpg)



D3Q19


![](images/de4bd631488192b6673a8fe7d644306cd8c07714e1170aa5f30d266d9cbae7f1.jpg)



D3Q27


![](images/60c5760f50d5b848166592672f4b73c374f231f6f342b8514a3a63e304a6dc08.jpg)



FeatFlow


![](images/a5b7b29c53597a911b1adf0d3eafb18353ed52fc55a92b079573026e5923446a.jpg)



(a)


![](images/55ccf8296c1232567532d027f018f20e27a25817c7eb950f96e4e76350d22717.jpg)


![](images/6b9ddb64f68d979ff45a1d5f61b6e915edc3849d9c24fabc2055ec901de1667c.jpg)


![](images/c03ec032cc3319f4284ef927b80492e09336a05143aadba7b83efc2271607f97.jpg)


![](images/18a9b51763a2b0ddc1a460d84836806284de2df19014378ddfdf1e43f4d2cd7f.jpg)



(b)



Fig. 15: Bubble interface T = 3.0 in 3D problem of TC2, (a) exterior view and (b) middle section view.


![](images/02a6d3a448a553814bcc62d74d311fd0698235bda68842a13ee8a913e9daa494.jpg)



D3Q15


![](images/c3d124313a7ada4abd7bab2a88291e2b4c8d40d4f304a742e44d741530e1d220.jpg)



D3Q19


![](images/a6c13e1351679e61a51186da6ecc5bd83b338c33d1f23df0edcea93aa3835a75.jpg)



D3Q27


![](images/82281de72dec56171d96227f5eab587ba33917313cbc4e61c5ce2a3e669058cd.jpg)



FeatFlow



Fig. 16: Bubble interface at T = 3 in 3D problem of TC2 for θ = 0 (red) and θ = 45 (blue) planes.


![](images/53ced69aed5de6c005d789740801d345e1d3d12ad2b28e08852bf88e4f64f9eb.jpg)


![](images/63f973765394a904f2eeeabb84a7d12433e9c031c4e3783f449e8cabcdbdf3cd.jpg)



Fig. 17: Time evolution of diameters in horizontal and vertical directions for 3D rising bubble TC2 (left), and enlarged view (right).


![](images/2832f7be693cc6032cbce37bda296122eb249fc07e9d6be95e77c2e8fba92232.jpg)


![](images/7711fe78b85f1817b64505e0fa6c0dd73208a24079ba2b0547d6f7fb84bf677f.jpg)



Fig. 18: Time evolution of rise velocity for 3D rising bubble TC2 (left), and enlarged views (right).


## 5 Numerical performance and stability in strong deformation singularities

After evaluating the accuracy of the 3D coupled scheme through benchmarking for rising bubbles, we may now proceed with simulations involving strong deformation singularities at high speeds so as to push the scheme towards its limits focusing more in the numerical stability and qualitative physical description. This includes the two problems of droplet splashing on a thin liquid film and the binary droplet collision. For this class of problems, a detailed quantitative comparison with other numerical and experimental data is beyond the scope of the current paper and only limited studies are presented to show the overall accuracy of the solutions. 

## 5.1 Droplet Splashing on thin liquid film

Splashing of a droplet over a thin liquid film is an attractive problem for both physicists and mathematicians []. In fact, the resulting splashing patterns could vary significantly with the falling velocity and the density of the droplet. The different splashing regimes are classified using the Reynolds Re and Weber W e numbers as 

$$
R e = \frac {\rho_ {1} U _ {0} (2 r _ {0})}{\mu_ {1}} \qquad , \qquad W e = \frac {\rho_ {1} U _ {0} (2 r _ {0}) ^ {2}}{\sigma}
$$

where $U _ { 0 }$ is the impact velocity of the droplet. Both 2D and 3D simulations are performed using D2Q9 and D3Q19 stencils, respectively. The 2D configuration as seen in figure 19 consists of a domain of size $2 \times 1$ where symmetry and periodic boundary conditions are used at $x = 0$ and $x = 2$ planes, respectively. The lower boundary is set to no-slip wall and a second order extrapolation is used at the top boundary. Density and viscosity ratios are set to 1000 and 40, respectively and the surface tension forces are realized via $W e = 8 0 0 0$ for all cases. Figure 20 shows the splashing patterns on a lattice resolution of $2 5 6 \times 5 1 2$ , where the time is measured based on $T = U _ { 0 } t _ { l b } / ( 2 r _ { 0 } )$ , with $r _ { 0 }$ in lattice units. For the low speed case of $R e = 2 0$ , the droplet will be slowly deformed and flattened upon the film. By increasing the Re number to 100 and 500, the splashing becomes fairly pronounced with a large rim and elevated fingers in the form of thin filaments. A widespread criteria to quantify the splashing effect is the spreading factor $R _ { s p } / 2 r _ { 0 }$ , where $R _ { s p }$ is the spreading radius at which the fluid velocity is maximum [20]. The numerical values of the spreading factor are compared in figure 21 against those in the work of Lee and Lin [23]. The straight line represents the power law trend $R _ { s p } = \sqrt { 2 r _ { 0 } U _ { 0 } t _ { l b } }$ proposed by Josserand and Zaleski [20]. In principle, agreement with the other numerical data is very close. Although the general trend of the power law is captured, one must note that the power law is only accurate at the very early times after the impact [23]. 

![](images/efbcbc1ae2b6f3c5eebcecbf95a620c07accd7db08919f8d1fb70e540b587635.jpg)



Fig. 19: Problem definition for the 2D droplet splashing on a thin liquid film.


Fro the 3D simulations, the solution covers only $1 / 4$ of the problem with symmetry boundary conditions on $x = 0$ and $y = 0$ planes, where the bubble hits the film along the z axis. The lattice resolution is $1 / h = 1 2 8$ , giving a domain size of $1 2 8 \times 2 5 6 \times 2 5 6$ . Figures 22 and 23 show the 3D frames extended in the x direction. It could be seen that the edges of the ring become more straight while the ring itself attains a smaller radial extension as compared to the 2D case. Moreover, the spin-off of the satellite bubbles from the thinning rim for $R e = 5 0 0$ is more clearly captured and visible in 3D. 

![](images/4ef99ceac801c308775d1624806b1670b5f0be94364c760ee3f72fdc8d0e52d4.jpg)



Re = 20


![](images/8ff5f0f1a336e31204d8fb9277a88ccceb067023c68671a3112d2c7784eba984.jpg)



Re = 100


![](images/080599fdb1b2fa62a00a85cc0bf52c04dc56250ea18b6c06550c3cf098602a80.jpg)



Re = 500



$T = 0 . 2$



T = 0.4



T = 0.8



T = 1.2



Fig. 20: Temporal evolution of the 2D droplet splashing over a thin film under different values of Re number.


![](images/057f3bf61a59c5b8d2bff80f8a9a711fca8aacbfe78c380be1e0e850eaf7844a.jpg)



Fig. 21: Plot of the spread factor ${ R _ { s p } } / { 2 r _ { 0 } }$ as a function of non-dimensional time $T = U _ { 0 } t _ { l b } / ( 2 r _ { 0 } )$ .


![](images/2440de123d1bc13c7353e453b22ad4bbae8e529e827ce88f8db22fe58ff3e1ed.jpg)


![](images/13ef26785374fe6afa3832f9d2edc473a86d5e2712d3cedf020d7e9a060b2774.jpg)


![](images/a787d59fcdfc7b1162cbf42afd37570bc2cf2fb9e3ce8b7359c3679e5aa58680.jpg)



Fig. 22: Temporal evolution of the 3D droplet splashing over a thin film for $R e = 1 0 0 .$ .


![](images/f8207a5b21a7e79aacb0ffb28711c97b5c1418d58a735d5798ca1a6d25a52943.jpg)


![](images/b6b9ca7a6e2772375ef37850db5cae9ffa5134678684dd83d3d2c2a686fe50c7.jpg)


![](images/53129b128615566dce89f4030ea1ea011350c40be3a3b22441181723fa2dabad.jpg)



Fig. 23: Temporal evolution of the 3D droplet splashing over a thin film for $R e = 5 0 0$ .


## 5.2 Binary collision of droplets

The intriguing problem of binary droplet collision has been extensively studied in the experimental works of Ashghriz [5] and Qian and Law [33] as well as in the numerical works of Schelkle [37], Inamuro et al. [19], Moqaddam et al. [27] and Wang et al. [45]. In particular, Ashgriz and Poo [4] characterized the collisions using the Weber W e and Reynolds Re numbers as in the case of droplet splashing, along with two additional parameters to specify the relative size and the impact angle of of the two droplets 

$$
\mathcal {D} = \frac {d _ {1}}{d _ {2}} \qquad , \qquad B = \frac {2 X}{d _ {1} + d _ {2}}
$$

where $d _ { 1 }$ and $d _ { 2 }$ are the diameters of the two droplets and X is the distance between the center of one droplet to the relative velocity vector placed on the center of the other droplet [19]. Here, it is assumed that the droplets are of the same size, i. e., $\mathcal { D } = 1$ and $d _ { 1 } = d _ { 2 } = D$ . The domain size is $3 D \times 3 D \times 6 D$ and the droplets are initially placed 2D apart. The density and viscosity ratios are 1000 and 100, respectively, and the droplets collide with identical, but opposite velocities of $U _ { 0 } / 2$ . Periodic boundary conditions are considered on all sides. For the case of a head-on collision, i. e., $B = 0$ , the symmetry allows to solve for only $1 / 4$ of the problem. As outlined in [4], for $5 0 0 < R e < 4 0 0 0$ the Reynolds number has no significant role in the dynamics of the impact and it is thus fixed at $R e = 5 0 0$ such that the velocity magnitudes will remain in the low M a regime on moderately fine lattices. The progress in time is measured in the units of $U _ { 0 } / D ;$ with D in lattice units 

![](images/a624706a8e7eb63899f6ac84dd17317e8c7824a0e49ed6f5e167cf708ab98ef0.jpg)


![](images/c9532e11e491b9ca4481ddcd20b3e3bba9bff3a16719d9e0ab41a9571701e96e.jpg)


![](images/09c2da12e5b011de9d8e27a65f979c4a55396cb167f8877796fdcbe51f4b7fb9.jpg)


![](images/144eee9419902dec20389996d7a8e45930db4eae9fcb0c9efd48a7b937b9c0f8.jpg)



(a)


![](images/2f6f68c3c2f56af35c26d0a2e0857a5ede50da77027c3ec41b3994f30def0ed1.jpg)


![](images/317eecff98b75f852b525740fa1e21461e50d21f6e67cbf1ae1a6b17cefe8e7d.jpg)


![](images/30d392ce69d91c404356f33f03e4f7d1c53f9b405c3293b33cc89ff5d9dd8555.jpg)


![](images/2cdc2b60b13e926c50e2b3fc18a7a8a2e1baba4afb125debc4338f509d591e2c.jpg)



(b)


![](images/5dc3d4ce7c84ec1638b35de0d6a8a2fa82bf27704c4a1039af51ceae0dca0880.jpg)


![](images/4ae57eb51729b86992c611a45e171463b4f8b08dfdc92ecbadeadeeafb8ab052.jpg)


![](images/f7dcefdb370f4b7c88dd7d4095bc38cbb6bbbc12c7e9abae4a8e4eb9cca97cb1.jpg)


![](images/f0565f2b9539b80dd54178061cbdf4980d28dd16bd792d39da38090041c0ca65.jpg)



(c)



Fig. 24: Head-on collision of micro-droplets with $B = 0 , R e = 5 0 0$ and (a) $W e = 1 5 ,$ , (b) $W e = 2 5 ,$ , (c) $W e = 4 0 .$ .


![](images/1347c8b530505af58094f44cf9771720e575b1325659c9f43d47eddcfa117618.jpg)



(a)


![](images/94b659de9e38598399eb7e4c86f53537d145e7d41a0509e1957433f0f6a3a726.jpg)



(b)



Fig. 25: Experimental snapshots of head-on collision of micro-droplets with (a) $W e = 2 3 ,$ (c) $W e = 4 0$ , adopted from [4]


![](images/972d4d9fecf52c7e7f2cd114a3f7f0e53f4fa04827085049c1b3fb34a0edf708.jpg)


![](images/057c11cebb10a77547e984e8ff051791902a96ed9fcb0f5f5eb60574027508f7.jpg)


![](images/5519e567f25e8b33124fe52ed815cdb1b904bdeb3c5f3966a0542defb848ce0e.jpg)


![](images/df29bcd2b30a749eacbdf33535810eab7c37e8fb0748bdba7ce47299c182151b.jpg)



Fig. 26: Oblique collision of micro-droplets at Re = 500, W e = 80 and $B = 0 . 8$ .


Figure 24 shows the evolution of the colliding droplets for $W e = 1 5 , 2 5 , 4 0$ on a domain size of $6 4 \times 6 4 \times 2 5 6$ . In particular, the evolution of the droplets for the cases with higher Weber numbers of $W e = 2 5 , 4 0$ exhibits a very close resemblance to the experimental data for $W e = 2 3$ , 40 adopted from [4] in figure 25. For the high surface tension case of $W e = 1 5$ , the collision pattern falls into the so-called coalescence regime. By increasing the Weber number to $W e \ : = \ : 2 5$ the weaker surface tension allows for an eventual spin-off of the merging droplets, after which, the droplets travel away from each other. This regime is therefore called reflective separation. $\mathrm { B y }$ further increasing the Weber number to $W e = 4 0$ the enhanced stretching causes the two round ends to detach from the middle section and break-up, forming a third satellite bubble. which remains in the middle and oscillates until becoming perfectly circular. 

To demonstrate the stability of the scheme at an even more non-regular arrangement, the simulation result for the case of an oblique collision with $W e = 8 0$ and $B = 0 . 8$ are illustrated in figure 26. The droplets partially collide, but keep traveling close to their initial path. The oblique collision, however, causes a twist and rotary motion which creates a long connecting filament. The filament then detaches from the main droplets and further breaks up into two satellite droplets. The obtained results closely repeat the numerical ones in the work of Inamuro et. al [19] for the exact same case. 

## 6 Computational performance

following the general optimization rules in [14, 21, 29, 42] for single phase flow solvers, analysis of instruction and memory workload of the GPGPU parallel code implies a number of low-effort optimization measures as follows: 

 Generate and use the directional forces on-the-fly and avoid saving and writing to and from the global memory. 

 Limit the calculations for directional force as well as for $\mathbf { u } \cdot \nabla \psi$ to central type $\mathrm { i f } \ | \nabla \psi | < 1 0 ^ { - 9 }$ for double precision (DP) and $| \nabla \psi | < 1 0 ^ { - 6 }$ for single precision (SP) and hence avoid the large spatial support of upwind and WENO schemes far from the interface. 

 In order to ease the instruction bottleneck caused by the the MRT collision, The MRT computatoins can be limited to $\psi > ( 1 - 1 0 ^ { - 9 } )$ and $\psi > ( 1 - 1 0 ^ { - 6 } )$ regions in case of DP and $\mathrm { S P }$ , respectively, thus leaving the rest of the domain to rely upon the low-cost SRT model. 

Similar to the so-called adaptive mesh refinement methods [10, 12], the last two points comprise an adaptive computation refinement approach which comes as an effective tool in the present case to ease the computations. Parallel computations are carried out on two high-end GPU-enabled compute machines. The first one is equipped with a Kepler K20Xm GPGPU, having a peak DP performances of 1.32 TFLOPS. For pure SP codes, a second machine equipped with GeForce GTX980 Ti gaming GPU with 5.6 TFLOPS in $\mathrm { S P }$ is employed. A comprehensive study on the errors associated with $\mathrm { S P }$ computations is done in [34], where the resulting maximum $L _ { 2 }$ error in rising bubble quantities is reported to be limited to $\approx 1 \%$ for a lattice as fine as $1 / h = 6 4 0$ . 

A summary of the simulation statistics using different stencils is presented in table 2 for the rising bubble TC1 until $T = 3$ , where the total memory size occupied on GPGPU in GB is denoted by MGPU, number of time steps by NTS and the total simulation time in seconds by TGPU. Note that the simulation times also include the intermediate post-processing periods. It is also noteworthy that the the present GPGPU implantations rely upon two sets of distribution functions for the flow LBE. Although this choice is not optimal with regards to memory consumption, to our knowledge, it guarantees a more comprehensive implementation outline and guarantees the maximum memory access bandwidth in the GPGPU code. Yet, efficient GPGPU implementations using only one set of distribution functions are pursued in [38] for simple single phase flow solutions. 


Table 2: Simulation statistics the 3D rising bubble using different lattice stencils


<table><tr><td>1/h</td><td>MGPU</td><td>NTS</td><td>TGPU</td><td>TGPU/NTS</td></tr><tr><td colspan="5">D3Q15</td></tr><tr><td>96</td><td>0.17</td><td>26508</td><td>170</td><td>0.006</td></tr><tr><td>192</td><td>1.43</td><td>108300</td><td>4290</td><td>0.039</td></tr><tr><td>256</td><td>3.40</td><td>193548</td><td>17797</td><td>0.091</td></tr><tr><td colspan="5">D3Q19</td></tr><tr><td>96</td><td>0.20</td><td>26508</td><td>207</td><td>0.007</td></tr><tr><td>192</td><td>1.65</td><td>108300</td><td>5329</td><td>0.049</td></tr><tr><td>256</td><td>3.93</td><td>193548</td><td>20774</td><td>0.107</td></tr><tr><td colspan="5">D3Q27</td></tr><tr><td>96</td><td>0.26</td><td>26508</td><td>292</td><td>0.011</td></tr><tr><td>192</td><td>2.1</td><td>108300</td><td>7018</td><td>0.064</td></tr><tr><td>256</td><td>5.0</td><td>193548</td><td>29550</td><td>0.152</td></tr></table>

The performance of the 3D simulations in Million lattice update per second (MLUPS) for the D3Q19 stencil is given in figure 27 for full and adaptive computations. We again refer to [34] for elaborate discussions on the eventual computational cost of individual GPU kernels of the coupled solver. In general, one should note that incorporating two-phase forcing terms along with interface capturing can greatly degrade the performance as compared to single phase LB implementation. As such, the present performances are still considered to be satisfactory on a single GPU, although one may achieve even higher performances through more in-depth performance optimization studies for two-phase flows. Moreover, apart from the overall fine MLUPS values, the performance gain caused by the adaptive computations could be as high as 100% for the SP code which shows the effective impact of the employed techniques. 

![](images/cdc433536075baa9e3432f6118bc940be78356206fe4199e2287679c8496efd6.jpg)



Fig. 27: Performance of the 3D coupled scheme using the D3Q19 velocity stencil.


## 7 Summary and conclusions

The coupled LB-level set two-phase method was utilized for the 3D benchmarking problems. Stability was improved through using MRT collision scheme and the LBE is coupled with a mass conserving interface capturing scheme. Compared to conventional two-LBE schemes the current method asks for considerably less amount of memory due while the computational performance is still satisfactory considering the complex force term calculations. From the numerical point of view, the method appeared quite successful in both qualitative and quantitative aspects compared to the reference solutions. In particular, simulations for the high Eo number TC2 emphasize the fact that LB-based two phase flow solvers could reach the accuracy of sharp interface Navier-Stokes-based solution if the interface region is properly refined. Since such coupled schemes do not need any iterative solver for the system of equations and could be efficiently parallelized, employing higher resolutions than those common in NSE solvers would certainly pay off in terms of memory and simulation time. Application of different discrete velocity models for the LBE revealed that the D3Q19 stencil provides the required isotropy at a reasonable extra cost compared to D3Q15, while the D3Q27 stencil asks for more than 50% longer simulation times and no noticeable accuracy improvement upon the D3Q19 stencil, at least for the present range of Re and Eo numbers. Finally, considering the general unique feature of LB-based solutions, e. g. suitability for irregular and micro-scale geometries and its nice parallel scalability, the new coupled LB-LS scheme further proves the method’s feasibility to tackle non-conventional problems in two-phase flow systems, e. g. liquid movement in underground reservoirs and in micro porous layers of fuel cells. In particular, simulations for the latter interesting case of water saturation and phase change in fuel cells are in order and would be the subject of our upcoming publications. 

This work was supported by the German Graduate School of Energy Efficient Production and Logistics of the state of Nordrhein-Westfalen at TU Dortmund and by the DFG through the grant TU102/53 (SPP 1740). 

## References



29. C. Obrecht, F. Kuznik, B. Tourancheau, and J.-J. Roux. A new approach to the lattice boltzmann method for graphics processing units. Computers and Mathematics with Applications, 61(12):3628–3638, 2011. 





30. E. Olsson and G. Kreiss. A conservative level set method for two phase flow. Journal of Computational Physics, 210(1):225– 246, 2005. 





1. Featflow software, 2016. http://www.featflow.de/en. 





31. Fedkiw R. Osher, S. Level Set Methods and Dynamic Implicit Surfaces. Springer, New York, First edition, 2003. 





2. J Adelsberger, P. Esser, M. Griebel, S. Groß, M Klitz, and A. R¨uttgers. 3d incompressible two-phase flow benchmark computations for rising droplets. pages 5274–5285, 2014. 





32. S. Osher and J.A. Sethian. Fronts propagating with curvature-dependent speed: Algorithms based on hamilton-jacobi formulations. Journal of Computational Physics, 79(1):12–49, 1988. 





3. L. Amaya-Bower and T. Lee. Single bubble rising dynamics for moderate reynolds number using lattice boltzmann method. Computers and Fluids, 39(7):1191–1207, 2010. 





33. J. Qian and C.K. Law. Regimes of coalescence and separation in droplet collision. Journal of Fluid Mechanics, 331:59–80, 1997. 





4. N. Ashgriz and J. Y. Poo. Coalescence and separation in binary collisions of liquid drops. Journal of Fluid Mechanics, 221:183–204, 12 1990. 





34. M. A. Safi. Efficient computations for multiphase flow problems using coupled lattice Boltzmann-level set methods. PhD thesis, Fakult¨at f¨ur Mathematik, Technischen Universit¨at Dortmund, 2016. https://eldorado.tu-dortmund.de/handle/ 2003/35115. 





5. Nasser Ashgriz. Handbook of Atomization and Sprays. Springer, 2011. 





35. M.A. Safi and S. Turek. Efficient computations for high density ratio rising bubble flows using a diffused interface, coupled lattice boltzmann-level set scheme. Computers and Mathematics with Applications, 70(6):1290–1305, 2015. 





6. A. Banari, C. Janßen, S.T. Grilli, and M. Krafczyk. Efficient gpgpu implementation of a lattice boltzmann model for multiphase flows with high density ratios. Computers and Fluids, 93:1–17, 2014. 





36. M.A. Safi and S. Turek. Gpgpu-based rising bubble simulations using a mrt lattice boltzmann method coupled with level set interface capturing. Computers and Fluids, 124:170–184, 2016. 





7. A. Banari, C.F. Janßen, and S.T. Grilli. An efficient lattice boltzmann multiphase model for 3d flows with large density ratios at high reynolds numbers. Computers and Mathematics with Applications, 68(12):1819–1843, 2014. 





37. M. Schelkle, M. Rieber, and A. Frohn. Comparison of lattice boltzmanim and navier-stokes simulations of three-dimensional free surface flows. volume 236, pages 207–212, 1996. 





8. M. Cheng, J. Hua, and J. Lou. Simulation of bubble-bubble interaction using a lattice boltzmann method. Computers and Fluids, 39(2):260–270, 2010. 





38. M. Schnherr, K. Kucher, M. Geier, M. Stiebler, S. Freudiger, and M. Krafczyk. Multi-thread implementations of the lattice boltzmann method on non-uniform grids for cpus and gpus. Computers and Mathematics with Applications, 61(12):3730– 3743, 2011. 





9. Grace J. R. Clift, R. and M. E. Weber. Bubbles, Drops and Particles. Academic Press, New York, First edition, 1978. 





39. K. Suga, Y. Kuwata, K. Takashima, and R. Chikasue. A d3q27 multiple-relaxation-time lattice boltzmann method for turbulent flows. Computers and Mathematics with Applications, 69(6):518–529, 2015. 





10. B. Crouse, E. Rank, M. Krafczyk, and J. T¨olke. A lb-based approach for adaptive flow simulations. International Journal of Modern Physics B, 17(1-2):109–112, 2003. 





40. M. Sussman and E. Fatemi. Efficient, interface-preserving level set redistancing algorithm and its application to interfacial incompressible fluid flow. SIAM Journal on Scientific Computing, 20(4):1165–1191, 1999. 





11. D. D’Humi`eres, I. Ginzburg, M. Krafczyk, P. Lallemand, and L.-S. Luo. Multiple-relaxation-time lattice boltzmann models in three dimensions. Philosophical Transactions of the Royal Society A: Mathematical, Physical and Engineering Sciences, 360(1792):437–451, 2002. 





41. N. Takada, M. Misawa, A. Tomiyama, and S. Fujiwara. Numerical simulation of two- and three-dimensional two-phase fluid motion by lattice boltzmann method. Computer Physics Communications, 129(1):233–246, 2000. 





12. A. Fakhari and T. Lee. Finite-difference lattice boltzmann method with a block-structured adaptive-mesh-refinement technique. Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 89(3), 2014. 





42. J. T¨olke and M. Krafczyk. Teraflop computing on a desktop pc with gpus for 3d cfd. Int. J. Comput. Fluid Dyn., 22(7):443–456, 2008. 





13. S. Geller, S. Uphoff, and M. Krafczyk. Turbulent jet computations based on mrt and cascaded lattice boltzmann models. Computers and Mathematics with Applications, 65(12):1956–1966, 2013. 





43. S. Turek, O. Mierka, S. Hysing, and D. Kuzmin. Numerical study of a high order 3D FEM–level set approach for immiscible flow simulation. In S. Repin, T. Tiihonen, and T. Tuovinen, editors, Numerical methods for differential equations, optimization, and technological problems, Computational Methods in Applied Sciences, Vol. 27, pages 65–70. Springer, 2012. 





14. M. Geveler, D. Ribbrock, D. G¨oddeke, and S. Turek. Lattice-boltzmann simulation of the shallow-water equations with fluid-structure interaction on multi- and manycore processors. Lecture Notes in Computer Science (including subseries Lecture Notes in Artificial Intelligence and Lecture Notes in Bioinformatics), 6310 LNCS:92–104, 2010. 





44. Hakon Wadell. Volume, shape, and roundness of quartz particles. The Journal of Geology, 43(3):pp. 250–280, 1935. 





15. X. He, S. Chen, and G.D. Doolen. A novel thermal model for the lattice boltzmann method in incompressible limit. Journal of Computational Physics, 146(1):282–300, 1998. 





45. Y. Wang, C. Shu, and L.M. Yang. An improved multiphase lattice boltzmann flux solver for three-dimensional flows with large density ratio and high reynolds number. Journal of Computational Physics, 302:41–58, 2015. 





16. X. He and L.-S. Luo. Lattice boltzmann model for the incompressible navier-stokes equation. Journal of Statistical Physics, 88(3-4):927–944, 1997. 





46. D. Zhang, K. Papadikis, and S. Gu. Three-dimensional multi-relaxation time lattice-boltzmann model for the drop impact on a dry surface at large density ratio. International Journal of Multiphase Flow, 64:11–18, 2014. 





17. X. He, X. Shan, and G.D. Doolen. Discrete boltzmann equation model for nonideal gases. Physical Review E - Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary Topics, 57(1):R13–R16, 1998. 





47. H.W. Zheng, C. Shu, and Y.T. Chew. A lattice boltzmann model for multiphase flows with large density ratio. Journal of Computational Physics, 218(1):353–371, 2006. 





18. S. Hysing, S. Turek, D. Kuzmin, N. Parolini, E. Burman, S. Ganesan, and L. Tobiska. Quantitative benchmark computations of two-dimensional bubble dynamics. International Journal for Numerical Methods in Fluids, 60(11):1259–1288, 2009. 





19. T. Inamuro, T. Ogata, S. Tajima, and N. Konishi. A lattice boltzmann method for incompressible two-phase flows with large density differences. Journal of Computational Physics, 198(2):628–644, 2004. 





20. C. Josserand and S. Zaleski. Droplet splashing on a thin liquid film. Physics of Fluids, 15(6):1650–1657, 2003. 





21. F. Kuznik, C. Obrecht, G. Rusaouen, and J.-J. Roux. Lbm based flow simulation using gpu computing processor. Computers and Mathematics with Applications, 59(7):2380–2392, 2010. 





22. P. Lallemand and L.-S. Luo. Theory of the lattice boltzmann method: Dispersion, dissipation, isotropy, galilean invariance, and stability. Physical Review E - Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary Topics, 61(6 B):6546– 6562, 2000. 





23. T. Lee and C.-L. Lin. A stable discretization of the lattice boltzmann equation for simulation of incompressible two-phase flows at high density ratio. Journal of Computational Physics, 206(1):16–47, 2005. 





24. T. Lee and L. Liu. Lattice boltzmann simulations of micron-scale drop impact on dry surfaces. Journal of Computational Physics, 229(20):8045–8063, 2010. 





25. M. Mehravaran and S. Kazemzadeh Hannani. Simulation of buoyant bubble motion in viscous flows employing lattice boltzmann and level set methods. Scientia Iranica, 18(2):231 – 240, 2011. 





26. M. Mehravaran and S.K. Hannani. Simulation of incompressible two-phase flows with large density differences employing lattice boltzmann and level set methods. Computer Methods in Applied Mechanics and Engineering, 198(2):223–233, 2008. 





27. A.M. Moqaddam, S.S. Chikatamarla, and I.V. Karlin. Simulation of binary droplet collisions with the entropic lattice boltzmann method. Physics of Fluids, 28(2), 2016. 





28. S. Mukherjee and J. Abraham. A pressure-evolution-based multi-relaxation-time high-density-ratio two-phase latticeboltzmann model. Computers and Fluids, 36(6):1149–1158, 2007. 

