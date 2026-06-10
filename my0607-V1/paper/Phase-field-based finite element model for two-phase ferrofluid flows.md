RESEARCH ARTICLE | FEBRUARY 21 2024 

# Phase-field-based finite element model for two-phase ferrofluid flows 

Pengfei Yuan (袁鹏飞) ; Qianxi Cheng (程仟禧) ; Yang Hu (胡洋) ; Qiang He (何强)  ; Weifeng Huang (黄伟峰)  ; Decai Li (李德才) 

![](images/921fb916f36fb323cb2d109b9fded8049ddf0784bd4eadd152d4dbfb92c8f0ce.jpg)


Check for updates 

Physics of Fluids 36, 022016 (2024) 

https://doi.org/10.1063/5.0185949 

![](images/8ca087aa269c096be60ba8000b43373d1f523098443f58afc5e9492ed99bb224.jpg)



View Online


![](images/666b21503f08075c1e536c3585c3714a74fbb2754e498e1d453a6a35260bb43e.jpg)



Export Citation


# Articles You May Be Interested In

On the Rosensweig instability of ferrofluid-infused surfaces under a uniform magnetic field 

Physics of Fluids (November 2023) 

Rosensweig instability in ferrofluids 

Low Temp. Phys. (October 2011) 

An experimental study on Rosensweig instability of a ferrofluid droplet 

Physics of Fluids (May 2008) 

![](images/255b158547839f510ac302239a08d6eaccb48903e0d3fab7d5757b2f3a5854b1.jpg)


Physics of Fluids 

Special Topics Open for Submissions 

Learn More 

![](images/dd86966763a0c0d5005107d9132855cacebc38699455304884ef0dcf0657525e.jpg)


AIP Publishing 

# Phase-field-based finite element model for two-phase ferrofluid flows

Cite as: Phys. Fluids 36, 022016 (2024); doi: 10.1063/5.0185949 Submitted: 3 November 2023 . Accepted: 25 January 2024 . Published Online: 21 February 2024 

![](images/d01911d0c58ee5dd88e0beabd63ad48ef7ac249ff62edbd2bcb1912a408f0589.jpg)


![](images/9424820a1aa869e62c11cb486a654c72dcd9ebd8f263ce85e3630ace38bc7ff0.jpg)


![](images/34a4c9d34867ddda2df6b534c960c34ea1eaa5dabc1a0aaf54125d5cf699d9e3.jpg)


Pengfei Yuan (袁鹏飞),1,2 Qianxi Cheng (程仟禧),1,2 Yang Hu (胡洋),2 Qiang He (何强),1,a) Weifeng Huang (黄伟峰),1,a) and Decai Li (李德才) 1 

# AFFILIATIONS

1 State Key Laboratory of Tribology, Tsinghua University, Beijing 100084, China 

2 School of Mechanical, Electronic and Control Engineering, Beijing Jiaotong University, Beijing 100044, China 

a)Authors to whom correspondence should be addressed: heqiang@tsinghua.edu.cn and huangwf@mail.tsinghua.edu.cn 

# ABSTRACT

In this study, we propose a phase-field-based finite element model to simulate two-phase ferrofluid flows in two and three dimensions. The proposed model combines the Cahn–Hilliard equation to handle the phase field, the Poisson equation to account for magnetics, and the Navier–Stokes equation to characterize fluid flow. To efficiently handle this coupling, we present a linear, totally decoupled numerical scheme, which involves solving four separate equations independently, namely, a linear elliptic system for the phase function, a Poisson equation for the magnetic potential, a linear elliptic equation for the velocity, and a Poisson equation for the pressure. To assess the accuracy, applicability, and numerical stability of the model, we conduct simulations for several typical problems. These include investigating the deformation of a ferrofluid droplet under a two-dimensional uniform magnetic field model, the bubble coalescence in ferrofluids under a three-dimensional uniform magnetic field model, the collision of two ferrofluid droplets under two-dimensional shear flow, and the twodimensional interfacial instability of a ferrofluid. The numerical results confirm the model’s capability to robustly simulate multiphase flow problems involving high-density and high-viscosity ratios, both in two- and three-dimensional problems. Moreover, the model effectively captures fundamental phenomenological features of two-phase ferrofluid flows under large topological changes such as the Rosensweig instability. 

Published under an exclusive license by AIP Publishing. https://doi.org/10.1063/5.0185949 

# I. INTRODUCTION

Ferrofluid is a colloidal suspension composed of a layer of magnetic nanoscale particles coated with surfactants, dispersed in a carrier fluid such as water, oil, or biocompatible fluids.1 The stable suspension of these nanoparticles in a ferromagnetic fluid is achieved through the combined effects of the Brownian motion and surfactants. Ferrofluids can exhibit strong magnetization when subjected to an applied magnetic field. Due to their stability and controllability, ferromagnetic fluids find extensive applications in various engineering fields, including seals, lubricants, electronic cooling, vibration damping, adaptive optics, and microfluidic pumps.2 Hence, acquiring a comprehensive understanding of the behavior displayed by ferrofluid flow under the influence of external magnetic fields is crucial for advancing current applications and driving the development of novel technologies. 

Multiphase ferrofluid flows are complex multiphysical problems involving interfacial interactions between fluids and ferrofluids. Due to the opaqueness of ferrofluid, direct observation of its flow behavior is highly challenging.3 However, with the rapid advancement of computer hardware technology and the availability of efficient numerical methods, numerical simulations have emerged as a crucial alternative tool for gaining a deeper understanding of the multiphase ferrofluid flow phenomena. Currently, numerical models for simulating interfacial ferrofluid flows can be categorized into two main groups: time-independent models and time-dependent models.4 The time-independent models are developed based on the coupled system of the Maxwell equations, the Navier–Stokes equations, and the Young–Laplace equation. The force balance at the unknown free interface is described by the Young–Laplace equation, rendering it suitable solely for stationary free surface problems. Lavrova et al.5 proposed a decoupling time-independent algorithm for solving the Maxwell equation, the Young–Laplace equation, and the Navier–Stokes equation, giving applications of the model to numerical simulations of dissipative systems, rotating shaft seals, equilibrium shapes of ferrofluid droplets, and mode formation in the mode of instability of the ferrofluid layer normal field. Gollwitzer et al.6 used a similar approach to study the surface morphology of ferrofluids. However, their method can only be used for fixed free surface problems. In the context of timedependent models, the primary obstacle lies in accurately tracking the moving interface. Subsequently, numerous time-independent models were developed. The volume-of-fluid method (VOF),7 the level-set method (LSM),8 the finite element method (FEM),9,10 and the lattice Boltzmann method (LBM)4,11–13 have been employed to describe the evolution of the ferrofluid interface. Korlie et al.7 developed a VOFtype model to simulate bubbles of a non-magnetic fluid rising in a ferrofluid and droplets of a ferrofluid falling through a non-magnetic fluid. Then, they further extended the model with a multicolor function scheme to suppress numerical bubble merging and numerically investigated bubble aggregation in ferrofluids.14 Afkhami et al.15 used a VOF algorithm based on segmented linear interface reconstruction to simulate a ferrofluid droplet suspended in a viscous medium under a uniform magnetic field. Zhu et al.8 developed a finite volume model to study the equilibrium shape of ferrofluid droplets and used a levelset approach to track the interface between the fluids. Shi et al.16 proposed a two-dimensional model to simulate the dynamics of ferrofluid droplets. They developed a new interface-tracking method, called VOSET, to capture the evolution of the interface, which is a coupled VOF and level-set method. Habera and Hron17 resolved the challenge of accounting for interfacial magnetic forces in ferrofluids under an external magnetic field by using the characteristic level-set method. This allowed them to include the effects of additional magnetic terms in the Cauchy’s stress tensor within the incompressible Navier–Stokes equations. The researchers also ensured volume conservation by reinitializing the method. Ni et al.18 present a versatile numerical approach to simulating various magnetic phenomena using a level-set method. A novel two-way coupling mechanism between a magnetic field and a magnetizable mechanical system is the heart of their method. Maroofiazar et al.19 used the level-set method to build a model that can accurately determine the interface between two phases at any given moment. Nochetto et al.20 developed a phase-field model describing the flow behavior of a two-phase ferrofluid based on the Cahn– Hilliard equation and gave an energy-stable numerical format for this model. The model can capture the fundamental phenomenological features of ferrofluids, such as the Rosensweig instability. However, their model can only handle two-phase flows of matched density (or almost matched density). Hu et al.11 extended the lattice Boltzmann method (LBM) based on phase field to simulate the flow of multiphase ferrofluids. Subsequently, they proposed multiple conservative phasefield lattice Boltzmann models, which are applicable to immiscible and incompressible flows with arbitrary number of components21 and thermal capillary flows with large density ratios and thermal physical parameters.22 Ghaderi et al.23 used a two-dimensional hybrid approach combined of lattice Boltzmann and finite volume method to study the falling and coalescence of a pair ferrofluid droplets in the uniform magnetic field. Li et al.24 developed a fractional step-based multiphase lattice Boltzmann method (LBM) combined with magnetic field evolution solving to numerically study various surface deformations of ferrofluid droplets immersed in organic oil. Niu et al.25 proposed a method combining simplified lattice Boltzmann method (SLBM) and self-correcting techniques to investigate the motion, deformation, and clustering of ferrofluid droplets suspended in non-magnetic fluids under different magnetic field intensities. Majidi et al.26 numerically studied the dynamic deformation and rupture of composite ferrofluid droplets under shear flow and a uniform magnetic field using a hybrid model of the lattice Boltzmann method and the finite difference method. In general, the VOF method has good mass conservation properties, but it is not very accurate in curvature calculation and difficult to reconstruct the interface. In addition, the primitive level-set method lacks mass conservation property.16 Unlike the VOF method and the level-set method, the finite element method and the lattice Boltzmann method based on the phase-field method have received much attention in recent years for its physical basis and simplified computational procedure. Therefore, the phase-field method is used in this study. 

In simulating ferrofluid dynamics, the finite volume method and the lattice Boltzmann method are commonly employed numerical simulation techniques, while the finite element method is seldom utilized. Nonetheless, in the domains of computer science and mathematics, the finite element method is widely recognized for its exceptional accuracy and efficiency in numerical computations. The finite element method employs higher-order interpolation functions to approximate continuous solution surfaces, enabling the computation of partial differential equations. Due to its capability to handle complex shapes and geometries, the finite element method has demonstrated remarkable accuracy in numerous practical applications and is particularly effective in addressing intricate nonlinear problems.27 However, it exhibits high algorithmic complexity and requires a relatively steep learning curve. To address these challenges, the FEniCS software platform was launched in 2003. FEniCS is an open-source solver and software library for the finite element method, equipped with a comprehensive and consistent set of mathematical modeling and automated analysis tools. This platform simplifies the process of implementing and analyzing finite element methods, enhancing the overall efficiency and convenience of the numerical simulation workflow.28 In recent years, FEniCS has achieved great success in high-performance simulations of turbulent flows,29 as well as single-phase steady-state electrohydrodynamic simulations in nanopores30 and model cracks.31 

In this paper, a finite element model based on the phase-field method is constructed to solve the coupled equations of the two-phase ferrofluid flow in the framework of FEniCS. The accuracy, applicability, and numerical stability of the model are verified through simulations of typical problems such as the deformation of a ferrofluid droplet under a uniform magnetic field, the bubble coalescence in ferrofluids under a three-dimensional uniform magnetic field model, suspension of a ferrofluid droplet under shear flow, and interfacial instability of a ferrofluid. Numerical examples show that the model has good robustness for problems with high-density and high-viscosity ratios as well as large topological variations. 

The paper is structured as follows: Sec. II gives the control equations and numerical methods for ferrofluid systems. Section III gives numerical tests of the methods and simulations for several classical multiphase flow problems with ferromagnetic fluids. Some concluding remarks and a short discussion are provided in Sec. IV. 

# II. MATHEMATICAL MODELS

# A. Governing equations

The interface-tracking equation in this study is built upon the Cahn–Hilliard equation. The phase field / is used to denote the volume ratio of the ferrofluid, assuming two distinct values, namely, /  1:0 and /  1:0, to characterize the two different fluid components in the bulk. The Cahn–Hillard equation, which governs the evolution of the interface between the two fluids, can be expressed as follows:3 2,33 

$$
\partial_ {t} \phi + \mathbf {u} \cdot \nabla \phi = \nabla \cdot (M (\phi) \nabla G), \tag {1}
$$

where M / is the mobility coefficient and G is the chemical potential, which can be written as 

$$
G = \beta F (\phi) - \kappa \nabla^ {2} \phi , \tag {2}
$$

where $F ( \phi ) = 4 \phi ( \phi ^ { 2 } - 1 )$ , b, and j are the parameters associated with the surface tension and the interfacial thickness 

$$
\beta = \frac {3 \sigma}{4 \varepsilon}, \quad \kappa = \frac {3}{8} \varepsilon \sigma . \tag {3}
$$

According to Magaletti et al.,34 the mobility coefficient M / should be properly designed in order to reach the so-called “sharp-interface limit.” In this study, two different phase-field mobility expressions are employed35 

$$
M (\phi) = \varepsilon M _ {0}, \tag {4}
$$

$$
M (\phi) = M _ {0} \cdot \max \big [ (1 - \phi^ {2}), 0 \big ], \tag {5}
$$

where $M _ { 0 }$ represents the initial phase-field mobility. According to Magaletti et al.,34 the dimensional form of ${ \bf \dot { M } } _ { 0 }$ is given by 

$$
M _ {0} = c _ {0} \frac {\varepsilon^ {2} | U |}{\sigma}, \tag {6}
$$

where U is the characteristic velocity or velocity scale of the physical system, which can be determined by calculating the maximum velocity observed across the entire domain. Additionally, the safety factor c0 is introduced empirically to prevent the overestimation of mobility. 

The density q, the dynamic viscosity g, and the magnetic permeability l for the multiphase system can be calculated by a simple linear interpolation as 

$$
\rho (\phi) = \frac {\rho_ {1} + \rho_ {2}}{2} + \frac {\rho_ {1} - \rho_ {2}}{2} \phi , \tag {7}
$$

$$
\eta (\phi) = \frac {\eta_ {1} + \eta_ {2}}{2} + \frac {\eta_ {1} - \eta_ {2}}{2} \phi , \tag {8}
$$

$$
\mu (\phi) = \frac {\mu_ {1} + \mu_ {2}}{2} + \frac {\mu_ {1} - \mu_ {2}}{2} \phi , \tag {9}
$$

where the subscripts 1 and 2 represent the physical parameters of the two fluids, respectively. 

The magnetic governing equations of ferrofluid flows are Maxwell’s equations, which can be written as 

$$
\nabla \times \mathbf {H} = 0, \tag {10}
$$

$$
\nabla \cdot \mathbf {B} = 0, \tag {11}
$$

where is the magnetic field, is the magnetic induction, and, in H Baddition, the magnetic induction is related to the magnetic field as follows: 

$$
\mathbf {B} = \mu_ {0} (\mathbf {H} + \mathbf {M}) = \mu_ {0} (1 + \chi) \mathbf {H} = \mu \mathbf {H}, \tag {12}
$$

where l and l are the vacuum permeability and permeability, respectively, is the magnetization intensity, and v is the magnetic suscepti-Mbility. Based on the Langevin law, v can be expressed as a function of the magnetic field 1 

$$
\chi = \frac {| \mathbf {M} |}{| \mathbf {H} |} = \frac {M _ {s}}{| \mathbf {H} |} L \left(\frac {3 \chi_ {0} | \mathbf {H} |}{M _ {s}}\right), \tag {13}
$$

where Ms is the saturation magnetization intensity, $L ( x ) = \coth ( x ) - x ^ { - 1 }$ is the famous Langevin equation, and $\chi _ { 0 }$ is the initial magnetic susceptibility. According to Ref. 36, when $| \mathbf { H } | \ll M _ { s } ,$ it can be assumed that v  v0. 

As indicated by Eq. (10), to satisfy the conservative (irrotational) field condition of the magnetic field , a scalar potential w is introduced 

$$
\mathbf {H} = - \nabla \psi . \tag {14}
$$

By substituting Eqs. (12) and (14) into Eq. (11), we can derive the magnetic potential equation 

$$
\nabla \cdot (- \mu \nabla \psi) = 0, \tag {15}
$$

where the permeability l is related to the magnetic susceptibility v as $\mu = \mu _ { 0 } ( 1 + \chi ) .$ . 

The two-phase ferrofluid flow can be described by the incompressible Navier–Stokes equation, which incorporates supplementary force terms in addition to the body force resulting from the magnetic field and surface tension 

$$
\rho \left[ \frac {\partial \mathbf {u}}{\partial t} + (\mathbf {u} \cdot \nabla) \mathbf {u} \right] = - \nabla p + \eta \nabla \cdot [ \mathcal {D} \mathbf {u} ] + \mathbf {f} _ {s} + \mathbf {f} _ {b} + \mathbf {f} _ {m}, \tag {16}
$$

$$
\nabla \cdot \mathbf {u} = 0, \tag {17}
$$

where q and g denote the density and dynamic viscosity, respectively; denotes the velocity field; p presents the pressure; D $\begin{array} { r } { \mathbf { \Psi } = \frac { 1 } { 2 } ( \mathbf { V } \mathbf { u } + \mathbf { V } \mathbf { u } ^ { T } ) } \end{array}$ presents the rate of deformation, and s and $\mathbf { f } _ { b }$ uare u u f fthe surface tension and body force, respectively. A potential form of surface tension $\mathbf { f } _ { s } = - \phi \mathbf { V } G$ is employed. From the magnetic stress tensor $\tau _ { m }$ fproposed by Cowley and Rosensweig,37 the magnetic force $\mathbf { f } _ { m }$ is calculated by Hu et al., 

$$
\mathbf {f} _ {m} = \nabla \cdot \tau_ {m} = \frac {\mu - \mu_ {0}}{2} \nabla | \mathbf {H} | ^ {2}. \tag {18}
$$

According to the word,38 the diffusive flow of the Cahn–Hilliard model, i.e., M / $G, is related to the local compositions of the two components; then, the continuity equation for the binary system can be expressed as 

$$
\frac {\partial \rho}{\partial t} + \boldsymbol {\nabla} \cdot (\rho \mathbf {u}) - \frac {(\rho_ {1} - \rho_ {2})}{2} \boldsymbol {\nabla} \cdot (M (\phi) \boldsymbol {\nabla} G) = 0. \tag {19}
$$

Using Eqs. (19) and (17), the momentum equation (16) for multiphase flow can be rewritten as 

$$
\frac {\partial (\rho \mathbf {u})}{\partial t} + \boldsymbol {\nabla} \cdot \left(\rho \mathbf {u u} - \frac {\left(\rho_ {1} - \rho_ {2}\right)}{2} M (\phi) \boldsymbol {\nabla} G \mathbf {u}\right)
$$

$$
= - \boldsymbol {\nabla} p + \eta \boldsymbol {\nabla} \cdot [ \mathcal {D} \mathbf {u} ] + \mathbf {f} _ {s} + \mathbf {f} _ {b} + \mathbf {f} _ {m}. \tag {20}
$$

For the phase field, in the presence of fluid interaction with solid walls, the following boundary condition is adopted to enforce a predetermined contact angle h at a solid boundary:39 

$$
\hat {\mathbf {n}} \cdot \nabla \phi | _ {\Gamma} = - \sqrt {\frac {2 \beta}{\kappa}} \cos \theta_ {c} (1 - \phi^ {2}), \tag {21}
$$

where C is the boundary, and ^ is the unit vector representing the nornmal direction pointing outward from the solid boundary wall. 

For the magnetic potential, we employ the Neumann boundary condition, where the normal derivative of the magnetic field density $H _ { n , \Gamma }$ at the boundary is expressed as follows: 

$$
\hat {\mathbf {n}} \cdot \nabla \psi | _ {\Gamma} = H _ {n, \Gamma}, \tag {22}
$$

the phase-field chemical potential is subject to the no-flux boundary condition, which can be expressed as 

$$
\hat {\mathbf {n}} \cdot \nabla G | _ {\Gamma} = 0. \tag {23}
$$

For the velocity field, the no-slip condition is applied 

$$
\left. \mathbf {u} (t) \right| _ {\Gamma} = 0. \tag {24}
$$

# B. Discretization scheme

Based on the above content, the two-phase ferrofluid flow system is composed of several coupled equations, which we present collectively as follows: 

$$
\partial_ {t} \phi + \mathbf {u} \cdot \nabla \phi = \nabla \cdot (M (\phi) \nabla G),
$$

$$
G = \beta F (\phi) - \kappa \nabla^ {2} \phi ,
$$

$$
\nabla (- \mu \nabla \psi) = 0,
$$

$$
\frac {\partial (\rho \mathbf {u})}{\partial t} + \nabla \cdot \left(\rho \mathbf {u u} - \frac {(\rho_ {1} - \rho_ {2})}{2} M (\phi) \nabla G \mathbf {u}\right) \tag {25}
$$

$$
= - \boldsymbol {\nabla} p + \eta \boldsymbol {\nabla} \cdot [ \mathcal {D} \mathbf {u} ] + \mathbf {f} _ {s} + \mathbf {f} _ {b} + \mathbf {f} _ {m},
$$

$$
\nabla \cdot \mathbf {u} = 0.
$$

It is feasible to solve nonlinear systems using a monolithic nonlinear scheme, and this approach leads to the assembly and iterative solution of a large system matrix, resulting in intensive computations. In this study, a decoupled linear scheme is developed, which divides the processes of phase-field transport, magnetic field calculation, and hydrodynamic flow. This allows the governing equations to be solved sequentially and separately in a more computationally efficient manner. 

Let $\Omega \in R ^ { d }$ be the computation domain and C be the boundary, where d is the dimension. We can define the finite element subspaces: $\mathbf { U } _ { h } = \left\{ { \mathbf { u } } \in H ^ { 1 } ( \Omega ) \right\} ^ { d }$ is the finite element subspace for velocity, $\mathrm { P } _ { h } =$ $\{ p \in \tilde { L _ { 0 } } ^ { 2 } ( \Omega ) \}$ for pressure, $\Phi _ { h } = \{ \boldsymbol { \phi } \in H ^ { 1 } \mathrm { \bar { ( } } \Omega \mathrm { ) } \}$ for phase field, $\tilde { \mathrm { G } } _ { h } = \left\{ G \in \mathring { H } ^ { 1 } ( \Omega ) \right\}$ for chemical potential, and $\dot { \Psi _ { h } } = \left\{ \bar { \psi } \in H ^ { 1 } ( \Omega ) \right\}$ for magnetic potential. In this study, a first-order time scheme is employed in conjunction with the linear splitting scheme. During the computation, after obtaining the values of $( \mathbf { u } ^ { n } , \bar { p ^ { n } } , \phi ^ { n } , G ^ { n } , \psi ^ { n } )$ from uthe previous computational step, the Galerkin method is employed to solve $( { \mathbf { u } } ^ { n + 1 } , p ^ { n + 1 } , \phi ^ { n + 1 } , G ^ { n + 1 } , \phi ^ { n + 1 } )$ /nþ1 ; belonging to the space $\mathbf { U } _ { h } \times \mathbf { P } _ { h }$ $\times \Phi _ { h } \times  { \mathrm { G } } _ { h } \times \Psi _ { h }$ , where n is the time level. 

# 1. Solving the Cahn–Hilliard equation

The computational procedure begins by solving the Cahn– Hilliard equation. We can computer $( \phi ^ { n + 1 } , G ^ { n + 1 } )$ for $n \geq 0$ by 

$$
\frac {1}{\delta t} \left(\phi^ {n + 1} - \phi^ {n}\right) + \boldsymbol {\nabla} \cdot \left(\mathbf {u} ^ {n} \phi^ {n + 1}\right) - M (\phi) ^ {n + 1} \Delta G ^ {n + 1} = 0,
$$

$$
G ^ {n + 1} - \beta F (\phi^ {n}) - \beta F ^ {\prime} (\phi^ {n}) \left(\phi^ {n + 1} - \phi^ {n}\right) - \kappa \nabla^ {2} \phi^ {n + 1} = 0, \tag {26}
$$

$$
\partial_ {n} \phi^ {n + 1} | _ {\partial \Omega} = 0, \partial_ {n} G ^ {n + 1} | _ {\partial \Omega} = 0.
$$

Due to the nonlinear nature of the term $F \bigl ( \phi ^ { n + 1 } \bigr )$ , it is linearized around $F ( \phi ^ { n } )$ using the expression $F ( \phi ^ { n } ) + F ^ { \prime } ( \phi ^ { n } ) ( \phi ^ { n + 1 } - \phi ^ { n } )$ . Therefore, the additional term ${ \bf { \bar { \theta } } } - \beta F ^ { \prime } ( \phi ^ { n } ) ( \phi ^ { n + 1 } - \phi ^ { \dot { n } } )$ is introduced in Eq. (26). After acquiring the value of the phase field, denoted as $\phi ^ { n }$ , the physical parameters corresponding to the phase field at the n-th time step $( \rho ^ { n } , \eta ^ { n } ;$ , and $\mu ^ { n } )$ can be updated utilizing the following equations: 

$$
\rho^ {n} = \frac {\rho_ {1} + \rho_ {2}}{2} + \frac {\rho_ {1} - \rho_ {2}}{2} \hat {\phi} ^ {n},
$$

$$
\eta^ {n} = \frac {\eta_ {1} + \eta_ {2}}{2} + \frac {\eta_ {1} - \eta_ {2}}{2} \hat {\phi} ^ {n}, \tag {27}
$$

$$
\mu^ {n} = \frac {\mu_ {1} + \mu_ {2}}{2} + \frac {\mu_ {1} - \mu_ {2}}{2} \widehat {\phi} ^ {n},
$$

where a cutoff function is adopted to avoid unphysical parameters 

$$
\hat {\phi} = \left\{ \begin{array}{l l} \phi , & | \phi | \leq 1, \\ \operatorname{sign} (\phi), & | \phi | > 1. \end{array} \right. \tag {28}
$$

# 2. Solving the magnetic potential equation

Next, the magnetic potential equation is solved by 

$$
\nabla \left(- \mu^ {n + 1} \nabla \psi^ {n + 1}\right) = 0, \tag {29}
$$

$$
\partial_ {n} \psi^ {n + 1} | _ {\partial \Omega} = H _ {n, \Gamma}.
$$

# 3. Solving the Navier–Stokes equations

A significant challenge in numerically simulating incompressible flows arises from the coupling between velocity and pressure due to the incompressibility constraint. In this study, an operator splitting method is adopted and the fluid flow step into the following three substeps.4 0 

We start with a forward step to solve the predicted velocity $\tilde { \mathbf { u } } ^ { n + 1 }$ as 

$$
\rho^ {n} \frac {\tilde {\mathbf {u}} ^ {n + 1} - \mathbf {u} ^ {n}}{\delta t} - \eta^ {n + 1} \boldsymbol {\nabla} \cdot [ \mathcal {D} \tilde {\mathbf {u}} ^ {n + 1} ] + \boldsymbol {\nabla} (p ^ {n})
$$

$$
+ \boldsymbol {\nabla} \cdot \left(\rho^ {n + 1} \mathbf {u} ^ {n} \tilde {\mathbf {u}} ^ {n + 1} - \frac {\left(\rho_ {1} - \rho_ {2}\right)}{2} M (\phi) ^ {n + 1} \boldsymbol {\nabla} G ^ {n + 1} \tilde {\mathbf {u}} ^ {n + 1}\right)
$$

$$
- \mathbf {f} _ {S} ^ {n + 1} - \mathbf {f} _ {m} ^ {n + 1} - \mathbf {f} _ {b} ^ {n + 1} + \frac {1}{2} \tilde {\mathbf {u}} ^ {n + 1} \frac {\rho^ {n + 1} - \rho^ {n}}{\delta t}
$$

$$
+ \frac {1}{2} \tilde {\mathbf {u}} ^ {n + 1} \boldsymbol {\nabla} \cdot \left(\rho^ {n + 1} \mathbf {u} ^ {n}\right)
$$

$$
+ \frac {1}{2} \tilde {\mathbf {u}} ^ {n + 1} \frac {\left(\rho_ {1} - \rho_ {2}\right)}{2} \boldsymbol {\nabla} \cdot \left(M (\phi) ^ {n + 1} \boldsymbol {\nabla} G ^ {n + 1}\right) = 0,
$$

$$
\tilde {\mathbf {u}} ^ {n + 1} | _ {\partial \Omega} = 0. \tag {30}
$$

It should be noted that the last three terms in Eq. (30) represent a firstorder approximation of the term 

$$
\frac {1}{2} \left[ \frac {\partial \rho}{\partial t} + \boldsymbol {\nabla} \cdot (\rho u) - \frac {(\rho_ {1} - \rho_ {2})}{2} \boldsymbol {\nabla} \cdot (M (\phi) \boldsymbol {\nabla} G) \right] u. \tag {31}
$$

This term can vanish due to Eq. (19), and it is worth mentioning that adding this term, as suggested in Ref. 36, can lead to an improvement in the computational stability. 

We must require $\boldsymbol { \nabla } \cdot \mathbf { u } ^ { n + 1 } = 0 ,$ and this leads to a Poisson equaution for the pressure difference 

$$
\Delta \left(p ^ {n + 1} - p ^ {n}\right) = \frac {\rho_ {0}}{\delta t} \boldsymbol {\nabla} \cdot \tilde {\mathbf {u}} ^ {n + 1}, \tag {32}
$$

$$
\partial_ {n} p ^ {n + 1} | _ {\partial \Omega} = 0,
$$

with $\begin{array} { r } { \rho _ { 0 } = \frac { 1 } { 2 } \operatorname* { m i n } ( \rho _ { 1 } , \rho _ { 2 } ) } \end{array}$ 

After having computed $p ^ { n + 1 }$ from this equation, we can update the velocity 

$$
\rho^ {n + 1} \left(\mathbf {u} ^ {n + 1} - \tilde {\mathbf {u}} ^ {n + 1}\right) = \delta t \boldsymbol {\nabla} \left(p ^ {n + 1} - p ^ {n}\right), \tag {33}
$$

$$
\left. \mathbf {u} ^ {n + 1} \right| _ {\partial \Omega} = 0.
$$

These schemes result in a decoupled Cahn–Hilliard Navier–Stokes system for two-phase ferrofluid flow at each time step, comprising an elliptic system for the phase function, a linear elliptic equation for the velocity, and a Poisson equation for the pressure. 

The mobility coefficient $M ( \phi ) ^ { n + 1 }$ in Eq. (30) should be carefully selected to ensure both the accuracy and stability of the current model. The initial phase-field mobility, denoted as $M _ { 0 } ^ { n + 1 }$ in Eq. (5), is determined by utilizing the maximum velocity of the previous time step $\lvert u _ { \mathrm { m a x } } ^ { n } \rvert$ , 

$$
M _ {0} ^ {n + 1} = c _ {0} \frac {\varepsilon^ {2} | u _ {\max} ^ {n} |}{\sigma}. \tag {34}
$$

In this study, the safety factor is set to 1.0. Therefore, the phase-field mobility in the $( n + 1 )$ -th time step can be expressed as 

$$
M (\phi) ^ {n + 1} = c _ {0} \frac {\varepsilon^ {2} | u _ {\max} ^ {n} |}{\sigma} \cdot \max \big [ (1 - \phi^ {n} \phi^ {n}), 0 \big ]. \tag {35}
$$

It is important to highlight that the set of equations consisting of Eqs. (36) and (37) constitutes a nonlinear system of equations owing to the cubic polynomial characteristics of the function $\bar { \boldsymbol { F } } ( \phi )$ . Consequently, the determination of the phase field and the chemical potential necessitates a coupled solution. Although several linearization methods have been proposed, in this study, we employ an iterative algorithm to simultaneously solve the two nonlinear equations, ensuring enhanced accuracy in computation. The weak form of the equations is presented. Assuming $( \varphi , g _ { \varphi } ) \in \Phi _ { h } \times G _ { h }$ as the test functions, the weak form for the phase-field equation can be formulated as follows: 

$$
\left(\frac {\phi^ {n + 1} - \phi^ {n}}{d t}, \varphi\right) - \left(\mathbf {u} ^ {n} \phi^ {n + 1}, \nabla \varphi\right) + \left(M (\phi) ^ {n + 1} \nabla G ^ {n + 1}, \nabla \varphi\right) = 0, \tag {36}
$$

$$
\begin{array}{l} \left(G ^ {n + 1}, g _ {\varphi}\right) = \beta (F (\phi^ {n}) - F ^ {\prime} (\phi^ {n}) \big (\phi^ {n + 1} - \phi^ {n} \big), g _ {\varphi}) \\ + \kappa (\nabla \phi^ {n + 1}, \nabla g _ {\varphi}) - \frac {3}{4} \sigma \cos \theta_ {c} \int_ {\Gamma} (1 - \phi^ {2}) g _ {\varphi} d \Gamma . \tag {37} \\ \end{array}
$$

The term ${ \begin{array} { l } { { \frac { 3 } { 4 } } \sigma } \end{array} }$ is derived from the expression $\kappa \sqrt { \frac { 2 \beta } { k } }$ based on Eq. (3). 

For the Laplace equation, the discretized weak formulation states that we find $\psi ^ { n + 1 } \in \Psi _ { h } ^ { * }$ such that 

$$
\left(\mu^ {n + 1} \nabla \psi^ {n + 1}, \nabla V\right) = \int_ {\Gamma} \mu^ {n + 1} H _ {n + 1, \Gamma} V d \Gamma . \tag {38}
$$

Here, the $V \in \Psi _ { h }$ is the test function. 

We find $\tilde { u } ^ { n + 1 } \in U _ { h }$ such that for all $\nu \in U _ { h }$ 

$$
\left(\rho^ {n} \frac {\tilde {\mathbf {u}} ^ {n + 1} - \mathbf {u} ^ {n}}{\delta t}, \mathbf {v}\right) + \left(\left(\left(\rho^ {n + 1} \mathbf {u} ^ {n} - \frac {\left(\rho_ {1} - \rho_ {2}\right)}{2} M (\phi) ^ {n + 1} \boldsymbol {\nabla} G ^ {n + 1}\right) \cdot \boldsymbol {\nabla}\right) \tilde {\mathbf {u}} ^ {n + 1}, \mathbf {v}\right) + (2 \eta^ {n + 1} \mathcal {D} \tilde {\mathbf {u}} ^ {n + 1}, \mathcal {D} \mathbf {v}) - (p ^ {n}, \boldsymbol {\nabla} \cdot \mathbf {v})
$$

$$
- \left(\mathbf {f} _ {S} ^ {n + 1} + \mathbf {f} _ {m} ^ {n + 1} + \mathbf {f} _ {b}, \mathbf {v}\right) + \frac {1}{2} \left(\tilde {\mathbf {u}} ^ {n + 1} \frac {\rho^ {n + 1} - \rho^ {n}}{\delta t}, \mathbf {v}\right) - \frac {1}{2} \left(\rho^ {n + 1} \mathbf {u} ^ {n} - \frac {\left(\rho_ {1} - \rho_ {2}\right)}{2} M (\phi) ^ {n + 1} \nabla G ^ {n + 1}, \nabla \left(\mathbf {u} ^ {n + 1} \cdot \mathbf {v}\right)\right) = 0,
$$

with the Dirichlet boundary condition $\tilde { \mathbf { u } } ^ { n + 1 } | _ { \partial \Omega } = 0$ . In so-called pressure correction step, we find $p ^ { n + 1 } \in  { \mathrm { P } } _ { h }$ u such that for all $q \in { \mathrm { P } } _ { h }$ , 

$$
\big (\nabla \big (p ^ {n + 1} - p ^ {n} \big), \nabla q \big) \Delta = \frac {\rho_ {0}}{\delta t} \big (\nabla \cdot \tilde {\mathbf {u}} ^ {n + 1}, q \big). \tag {40}
$$

Finally, we correct the velocity to have zero divergence by finding $ { \mathbf { u } } ^ { n + 1 } \in \mathbf { U } _ { h }$ that for all $\nu \in U _ { h }$ , 

$$
\frac {1}{\delta t} \left(\rho^ {n + 1} (\mathbf {u} ^ {n + 1} - \tilde {\mathbf {u}} ^ {n + 1}), \mathbf {v}\right) = (p ^ {n + 1} - p ^ {n}, \nabla \cdot \mathbf {v}). \tag {41}
$$

The scheme described above involves solving three decoupled subproblems sequentially. These subproblems are all linear and can be efficiently solved using specialized linear solvers. The sparse linear system is tackled using the generalized minimum residual (GMRES) method, which is a preconditioned Krylov method. 

# III. RESULTS AND DISCUSSION

In this section, we discuss several standard numerical tests to validate the accuracy and capabilities of the proposed model for simulating ferrofluid multiphase flow. 

# A. Ferrofluid droplet deformation

Magnetron microfluidics that utilizes micrometer-sized ferrofluid droplets has found a variety of applications in fields such as cancer cell screening, condition diagnosis,44–46 and chemical engineering.47 

This technology’s small size, contactless operation, and programmable functionality make it highly desirable, and research on the deformation of ferrofluid droplets is thus of great academic importance and value. Earlier, Flament et al.3 successfully measured the surface tension of ferrofluid droplets by deforming them using a uniform magnetic field. Recent advancements in numerical computation techniques have prompted researchers to develop various computational models to study the magnetic response properties of ferrofluid droplets. These models include the LBM phase-field model and the mixed immersion interface and phase-field lattice Boltzmann model, which were presented in our previous study,4,11 as well as the level-set finite element model developed by Cunha et al.48 To verify the accuracy of these models, Flament’s experimental results were used. 

A water-based ferrofluid droplet, which is immiscible with the alcohol filling a cavity, is placed at the center of the cavity, as illustrated in Fig. 1, and a uniform magnetic field in the y direction is applied. The combined effect of surface tension and magnetic field force causes the deformation of the ferrofluid droplet. To validate the model, the same parameters are used as in Flament’s experiment.3 The waterbased ferrofluid and alcohol have densities of $\rho _ { 1 } = 1 . 5 8 \times 1 0 ^ { 3 }$ and $\rho _ { 2 } = 0 . 8 \times 1 0 ^ { 3 } \mathrm { k g } / \mathrm { m } ^ { 3 }$ , respectively, while their dynamic viscosities are $\eta _ { 1 } = 1 6 \times 1 0 ^ { - 3 }$ and $\eta _ { 2 } = 0 . 8 \times 1 0 ^ { - 3 }$ Pas, respectively. The surface tension coefficient of the ferrofluid is $3 . 0 7 \mathrm { m N } / \mathrm { m } ^ { 2 }$ , and the initial magnetic susceptibility value is $\chi _ { 0 } = 2 . 2 .$ The saturation magnetization $M _ { s }$ is $4 0 \mathrm { k A } / \mathrm { m }$ . Two dimensionless parameters, the Ohnesorge number and the magnetic Bond number, are used to characterize this multiphase problem. These are defined as $\begin{array} { r } { \mathrm { O _ { h } } = \frac { \eta } { \sqrt { \rho R \sigma } } , \mathrm { B o _ { M } } = \frac { \mu _ { 0 } H _ { 0 } ^ { 2 } R } { 2 \sigma } . } \end{array}$ 

![](images/21d7780b8a1b984243ebfa62b07b28c7acfa3ce6996db84ab2106c7d484eab08.jpg)



FIG. 1. Schematic diagram of ferrofluid droplet deformation.


As the saturation magnetization is much greater than the external magnetic field in the simulation, the magnetic susceptibility is fixed at $\chi _ { 0 } = 2 . 2 $ The computational domain is a 10 mm square area, and the ferrofluid droplet has a radius R of 1 mm. The crossed rectangle grid (a grid composed of intersecting lines and rectangles, where the lines cross over each other and divide the rectangles into smaller subregions) size is set to $2 4 0 \times 2 4 0$ . The interfacial thickness e spans across two grid cells. The top and bottom boundaries of the computational domain have a no-slip solid wall boundary condition, and the uniform magnetic field condition $( \partial \psi / \partial y = H _ { 0 } )$ is applied to both boundaries. The periodic boundary condition is used for the x direction. 

Figure 2 displays the shapes of ferrofluid droplets at equilibrium under six different magnetic field conditions, with $\mathrm { O _ { h } }$ at 0.016, corresponding to a magnetic Bond number of 0, 0.296, 1.182, 1.726, 2.81, and 6.208 for ${ \mathrm { B o } } _ { \mathrm { M } } ,$ , respectively. The upper row of the figure presents the experimental results used for comparison,3 the middle row exhibits the simulation results of our previous study,11 and the bottom row depicts the results of the simulation calculations presented in this paper. Applying a uniform magnetic field in the vertical direction results in the elongation deformation of the ferromagnetic droplet in the y direction due to the magnetic field force, while the deformation of the droplet is prevented by surface tension. Consequently, the shape assumed by the droplet at equilibrium relies on the ratio between the two factors.11 Therefore, as depicted in Fig. 2, a weak applied magnetic field $( H { = } 1 . 2 \mathrm { k A / m } )$ results in an approximately circular ferrofluid droplet. However, as the applied magnetic field strength increases, the droplet shape changes from circular to elliptical. $\mathrm { A t } \ H { = } 5 . 5 \mathrm { k A / m } ,$ , it can be observed that the half-length axis b of the ferrofluid droplet is significantly larger than the half-short axis a. Second, the comparison of droplet shape results at equilibrium reveals good agreement between the simulation calculations presented in this paper and the experimental results. Furthermore, owing to the adaptive phase-field mobility and the inherent high accuracy of the finite element method, in comparison with our prior study,11 the droplet shape computed in this paper is sharper at the ends when subjected to a larger uniform magnetic field, which better matches the experimental droplet shape results. To quantitatively validate the accuracy of the model calculation proposed in this paper, Fig. 3 presents a comparison of the numerical calculation results of the aspect ratio b/a at equilibrium with the experimental results and our prior study11 calculation results. It is observed that the calculation results comply well with the experimental results, thus confirming the accuracy of the proposed model calculations. 

# B. Three-dimensional modeling of bubble coalescence in a ferrofluid under a uniform magnetic field

The formation of gas compartments due to improper injection of the ferrofluid, leading to the mixing of bubbles during the application of a magnetic field, can hinder direct contact between the ferrofluid and the sealed medium, adversely affecting the sealing effect. Therefore, investigating the coalescence deformation of bubbles in ferrofluids holds significant academic merit. Zheng et al.49 utilized numerical simulation methods to investigate the coalescence of two bubbles in conventional fluids. Their calculations indicated that coalescence was caused by surface tension only when the gap between two bubbles was less than twice the thickness of the interface. In contrast, bubble coalescence phenomena in ferrofluids are significantly more likely to occur and are independent of the distance between the bubbles.14 In addition, the fusion of bubbles results in significant interfacial topological changes, especially when dealing with two fluids of significantly different densities. As a result, the computational model faces a greater challenge. To address this, this subsection aims to conduct a numerical analysis of fusion between two bubbles submerged in a ferrofluid under a uniform magnetic field, using the model developed in this study. By doing so, we can further examine the computational stability of our proposed three-dimensional model with regard to problems related to large deformation and density ratios. Two dimensionless parameters, the Ohnesorge number and the magnetic Bond number, are used to characterize this multiphase problem. These are defined as Oh ¼ g p ; $\begin{array} { r } { \mathrm { O _ { h } } = \frac { \eta } { \sqrt { \rho R \sigma } } , \mathrm { B o _ { M } } = \frac { \mu _ { 0 } H _ { 0 } ^ { 2 } R } { 2 \sigma } , } \end{array}$ 

The ferrofluid used in this study is water-based and contains a 2% volume of 10 nm Fe O nanoparticles. The densities of the ferrofluid and gas are $\rho _ { 1 } = 1 1 0 0$ and $\bar { \rho _ { 2 } } = 1 . 2 9 3 \mathrm { k g } / \mathrm { m } ^ { 3 }$ , respectively, while the ferrofluid and gas have dynamic viscosities of $\eta _ { 1 } = 5 . 0 \times 1 0 ^ { - 3 }$ and $\eta _ { 2 } = 1 . 7 9 \times 1 0 ^ { - 5 }$ Pas, respectively. The surface tension coefficient is 5.0 mN=m2, and the magnetic susceptibility is constant at the initial value $\chi _ { 0 } = 1 . 5 1$ . The saturation magnetization $M _ { s }$ is 14 kA/m. Simple calculations give $\rho _ { 1 } / \rho _ { 2 } = 8 5 0 . 7$ and $\eta _ { 1 } / \eta _ { 2 } = 2 7 9 . 3$ . 

![](images/5af61967df7573e5db91ba6c63f71788e82553baae10950dbda359ef346bdcf7.jpg)


![](images/0983085ca5e612f62542f9daf20eadf075ad1b33786eb7b77879f76a007ce13b.jpg)


![](images/fd1be6beb27ce0fcd527266c2e35d05f072eeca482b7a7be1d265c6ea0698f20.jpg)


![](images/397a6d9f5bb945f8eeb8e7848956c2d7d1015bdee77e1ab014c26c295ad551a1.jpg)


![](images/d5b3297320e346c8fafd88803e8c565b37ef18e74e3e417f544ba0467461a55a.jpg)


![](images/790e9fd43d498ca8a39102096d34af9caf0c13decf5d8d51ba1c08625dcc64f5.jpg)


![](images/ee94c3252f2ca6a265efcb2aa46c8648d4503062324faed569792263b7daf41d.jpg)


![](images/6ec99cec368ee3bcdd0c5f68719de3724d09bd480341e29efff6bb5668c012d8.jpg)


![](images/5052f31da8ac32730f83900d6c8f6cf1219aa96cdb4d1f797febb6bef965a8a8.jpg)


![](images/1f9a00be4d462a1702314358ff272381480d7d1dacbf416dbb4d2d95d06ff7d3.jpg)


![](images/cc5d0bbb94ce72441b26b97289a89e5fca765f560808b537785bddb208f774a3.jpg)


![](images/73734f63531e22646bc70c7c3e938773e88ccc9a498333b57ab2d9cf3196a459.jpg)


![](images/f7ac2669dd605d0cdb571c68171e725cc635a77b547804ec6bb608cffb244ebf.jpg)


![](images/5b9d136eebff0c4b8703cfc4ab0275e02b7fc7bdd1a7c895fa9da36165483108.jpg)


![](images/2f42b7a12c17e3b973053f6e4ec9a4ed6bd2d4934317eea87776310a7710ccbf.jpg)


![](images/0ae5d3d79f5ed1f0330dacaaa4d4ee1955edd4d72c1b8400d2f059ea1a457f8a.jpg)


![](images/8a7c820d111aebbddfa9de7e465b358c9482ccfee8e9cbcca519119297958823.jpg)


![](images/48c2a34494c5b0782ba6ce9c2fafeac25701bf1faf4f51fae320dd7a00641d54.jpg)


![](images/c583d322c0253dea04960f778d193d0dd4187f95469aed1132ba4b7fa8de870a.jpg)



FIG. 2. Comparative analysis between experimental findings and numerical simulations depicting the morphology of ferrofluid droplets subjected to varying magnetic fields at Oh ¼ 0.016, BoM ¼ 0, 0.296, 1.182, 1.726, 2.81, and 6.208. The top row represents experimental results of the literature,3 the middle row presents the phase-field simulation results from our previous research,11 and the bottom row displays the phase-field simulation results from this study.


Figure 4 shows the schematic diagram of the problem with a computational domain of 6 6 6 mm3. There are two bubbles of radius R 0.5 mm placed at (3,3,2) and (3,3,4), respectively. The interface thickness is set to $\varepsilon = 6 \times 1 0 ^ { - 4 }$ mm. The gap d  1 mm between the two bubbles is significantly greater than twice the interface thickness. The ortho-hexahedral grid size is set to 40 40 40. The noslip boundary condition is implemented on the top and bottom boundaries, while periodic boundary conditions are used on the lateral boundaries. Magnetic field conditions of $\partial \psi / \partial y = H _ { 0 }$ are imposed at the top and bottom of the magnetic field domain. 

Figure 5 clearly illustrates the dynamic evolution process of the merging of two bubbles within a ferrofluid under the influence of a uniform magnetic field H 12 kA/m, with $\mathrm { O _ { h } }$ at 0.095 and $\mathtt { B o } _ { \mathrm { M } }$ at 9.072. During the initial stage, the superposition of magnetic field variations caused by each individual bubble leads to the overall nonuniformity of the magnetic field, resulting in the rapid elongation of the bubbles toward one side. The combined effect of interfacial magnetic forces and surface tension creates a region of low pressure between the bubbles, promoting their mutual approach until the tip of one bubble comes into contact with the other, forming a connecting 

bridge and rapidly merging into a larger bubble that eventually adopts a stable shape. In contrast to the work of Chen et $a l . ^ { 5 0 }$ on nonmagnetic liquid bubble fusion, in the present simulation, once the two smaller bubbles merge into a larger one, high magnetic field regions form at the top and bottom ends of the bubble, effectively suppressing the inertial effects caused by surface tension during the bubble fusion process, thereby eliminating noticeable oscillations. To provide a quantitative representation of this process, Fig. 6 showcases the timeevolution of the system’s total kinetic energy, which can be calculated by $\begin{array} { r } { E = \frac 1 2 \sum _ { i , j } \rho _ { i , j } \dot { u } _ { i , j } ^ { 2 } \Delta x ^ { 2 } } \end{array}$ over the entire computational domain. As depicted, the kinetic energy exhibits a clear peak, reaching its maximum value at approximately 0.044 s. Furthermore, in comparison with previous studies on two-dimensional droplet merging,4 the merged droplet in this case is significantly influenced by surface tension, resulting in a final shape that is more closely spherical in nature. 

In addition, it should be noted that the merging process leads to a significant topological change. The conservation of mass within the model is confirmed by evaluating the total mass of the system. The total mass of the system, $M _ { t }$ is calculated using the following equation: 

![](images/560c0d03b0691f91687e15a6d1cd8f3f1c4296950ea255e97836d583ec71ed5f.jpg)



FIG. 3. Comparison between the present numerical calculations, the previous numerical calculations reported by Hu et $a / . ^ { \dag \dag }$ and experimental results3 for the aspect ratio (b/a) of ferrofluid droplets exposed to magnetic fields of $H { = } 1 . 2 , 2 . 4 ,$ , 2.9, 3.7, and 5.5 kA/m.


![](images/7ec09caebdc6cfd3809af0065709c8403c05726327d0dcc48e17776a651978af.jpg)



FIG. 4. Schematic diagram of bubble fusion.


$$
M _ {t} = \sum_ {i, j} \phi_ {i, j} \Delta x ^ {2}. \tag {42}
$$

Figure 6 displays the evolution of the dimensionless total mass $\mathbf { M } / \mathbf { M } _ { 0 }$ at a magnetic field strength of $H { = } 1 2 \mathrm { k A / m }$ , where $\mathrm { { M } } _ { 0 }$ is the initial total mass. As can be seen, the total mass loss of the present model is within 0.004%. In cases where bubbles merge and cause a significant topological change, the mass curve does not exhibit a noticeable abrupt change. Instead, the observed mass variations can be attributed to minor decreases resulting from unavoidable computational errors. Thus, this model has a considerable advantage in dealing with problems that involve high viscosity, high mass ratios, and substantial topological changes. 

# C. Rheological properties of ferrofluid droplets suspended in shear flow

Recent research indicates that ferrofluid droplets can be employed as magnetically controllable materials to modulate suspension rheology, showing high potential for applications in magnetothermal therapy, convective heat transfer fluxes, and enhanced oil recovery.51 This section investigates the rheological properties of a dilute ferrofluid droplet solution under shear flow using the current model. Figure 7 illustrates a computational domain consisting of a viscous ferrofluid droplet located at the center, which is subject to a uniform magnetic field. The droplet is surrounded by another viscous fluid medium, where the upper and lower two planar walls move in opposite directions at the same velocity. A shear flow with a constant shear rate, $\gamma ,$ is applied in this system. The flat-plate moving velocity is obtained as $\begin{array} { r } { u _ { t o p } = \frac { 1 } { 2 } \gamma \mathrm { L } , 0 \big ) } \end{array}$ . The left and right boundaries of the computational domain adopt the periodic boundary condition. The orientation angle, $\theta ,$ of the droplet is defined as the positive angle between the x axis and the droplet’s main axis. In this shear flow, the suspension’s relative viscosity is calculated as 

$$
\eta_ {r} = \frac {\eta_ {s}}{\eta} = \frac {S / \gamma}{\eta}, \tag {43}
$$

where $\eta _ { s }$ is the effective viscosity of the suspension, g is the viscosity of the continuous phase, and S is the average shear stress. According to the previous studies,4,52,53 the streamlines of the Couette flow are in x direction and the velocity gradient is in y direction, they will vary depending on the viscosity ratios, and the average shear stress can be obtained by averaging the shear stress acting on the moving flat wall 

$$
S = \eta \big (\partial_ {y} u + \partial_ {x} v \big), \tag {44}
$$

where u and v are the components of the velocity vector in the x and y directions. Three dimensionless parameters, Reynolds number, magnetic Bond number, and capillary number, are used to characterize this multiphase problem. These are defined as $\begin{array} { r } { = \frac { \mu _ { 0 } H _ { 0 } ^ { 2 } R } { 2 \sigma } } \end{array}$ l0 H20 R , and Ca  gcR . 。 $\mathrm { C a } = { \frac { \eta \gamma R } { \pi } }$ r $\begin{array} { r } { \mathrm { R } _ { \mathrm { e } } = \frac { \rho \gamma R ^ { 2 } } { 4 \eta } } \end{array}$ 4g ; BoM 

To validate the proposed model for calculating the relative viscosity, the data of droplets subjected to shear flow in the absence of a magnetic field were compared with the previous work of Ghigliotti et $a l . , ^ { { \breve { 5 } } 4 }$ Cunha et $a l . , ^ { 4 8 }$ and He et $a l . ^ { 4 }$ The variation of relative viscosity with viscosity ratio a for $\mathrm { R e } = 0 . 0 1$ and $\mathrm { C a } = 0 . 3$ is plotted in Fig. 8. The present results are in good agreement with them, and the reason for the slight difference is due to the difference between the phase-field method and the level-set method, and the difference in the choice of the thickness of the phase-field interface. 

In this study, the viscosities of the ferrofluid droplet and the surrounding fluid are $\eta _ { 1 } = 1 \times 1 0 ^ { - 3 }$ and $\eta _ { 2 } = 8 \times 1 0 ^ { - 3 } \mathrm { P a } { \mathrm { s } } ,$ the densities are $\rho _ { 1 } = 1 5 8 0$ and $\rho _ { 2 } = 8 0 0 \mathrm { k g } / \mathrm { m } ^ { 3 }$ , respectively, and the magnetic susceptibility v is 2.2 and 0, respectively. The surface tension coefficient between the two phases is $\sigma = \mathrm { \bar { 3 } . 0 7 } \mathrm { m N } / \mathrm { m } ^ { 2 }$ . The size of the computational domain is set to $\mathrm { L } = 1 0 \mathrm { m m } .$ , and the initial radius of the ferrofluid droplet is R 1.0 mm. The computational domain is discretized by $2 4 0 \times 2 4 0$ crossed rectangle grids. The interfacial thickness e spans across two grid cells. To study the effect of magnetic field on the rheological properties of ferrofluid droplet suspension, several cases with different magnetic field strengths were studied and simulated for 1.0, 1.5, 2.0, 2.5, and 3.0 kA/m, with Re and Ca fixed at 0.375 and 0.39, respectively, corresponding to a Bond number of 0.205, 0.462, 0.82, 1.28, and 1.845 for ${ \mathrm { B o } } _ { \mathrm { M } } ,$ respectively. In this simulation, we introduce the concept of effective flow resistance height h, which is the 

![](images/80728f0b61b77aab62ed1d1ac296db57e8d2352da21034be3510bc21e7c45235.jpg)



FIG. 5. Two bubbles merging within a ferrofluid under a uniform magnetic field $( H = 1 2 \mathsf { k A } / \mathsf { m } )$ with $0 _ { h }$ at 0.095 and $B o _ { M }$ at 9.072. The left halves of the figure represent the velocity vectors, and the right halves display the distributions of the magnetic field intensity H.


![](images/33cf0ab2e87ae8bd5fd2594fcd0e40e3d8fceed65164fab90fafe772019cf9e2.jpg)



FIG. 6. Investigation of the time-dependent characteristics of total mass and kinetic energy during bubble coalescence in a ferrofluid with $0 _ { \mathfrak { h } } = 0 . 0 9 5$ and ${ \sf B o } _ { \sf M } = 9 . 0 7 2 .$ . The inserted contour plot represents the phase boundary at $\phi { = } 0 ,$ , effectively highlighting the bubble state at a specific time point.


![](images/bcbc0fabd0766939b4960587af6c8e852881a0c363690125b3a8c58214c2b290.jpg)



FIG. 7. Schematic diagram of droplet under shear flow.


elongation of a magnetic droplet in the direction normal to the translational velocity of the plate. 

Figure 9 depicts the impact of magnetic fields on droplet rheology. The graph reveals that an increase in magnetic Bond number leads to more droplet stretching, which was also consistent with previous research.48 More droplet stretching means that a larger effective flow resistance height is produced, which leads to an increase in the mean shear stress and thus an increase in relative viscosity. As shown in the chart, raising the magnetic field concentration causes an increase in both larger effective flow resistance height h and relative viscosity $\eta _ { r }$ . For the magnetic field strength surpassing $1 . 5 \mathrm { k A } / \mathrm { m } ,$ droplet effective flow resistance height and relative viscosity curves increase even more. However, the droplet orientation angle initiates an opposite shift. These observations suggest that at sufficient magnetic field strength, the orientation of the droplet changes slightly with magnetic strength, yet significant effective flow resistance height increases arise. Longer droplets produce significantly increased impedance to the shear flow, thus causing a greater effective viscosity. Conclusively, these findings indicate that external magnetic forces causing topological changes in ferromagnetic fluid droplets can efficiently control the rheological properties of emulsions. 

![](images/d20d0ec5a5a35914ecea4894cd98e0d342be8aa635694af27e87029559a9933f.jpg)



FIG. 8. Variation of relative viscosity gr with viscosity ratio a for $\mathsf { R e } = 0 . 0 1$ and Ca  0.3 without magnetic field.


Another simulation is conducted under a steady magnetic field of 2 kA/m, where shear rates c of 5, 10, 15, 20, 25, and 30 were chosen to study their influence on droplet’s rheological properties. The results, as indicated in the Fig. 10, revealed a decrease in the droplet orientation angle and an increase in both effective flow resistance height and relative viscosity as the shear rate escalated. Contradictory to the impact of magnetic field strength on droplets’ rheology, the progression of the relative viscosity curve tends to plateau with increasing shear rate. Figures 9 and 10 illustrate a strong congruity between the trends of relative viscosity and effective flow resistance height curves. On this basis, we infer that the rheological characteristics of an emulsion can be efficiently regulated by adjusting the effective flow resistance height of its ferrofluid droplets. 

To further examine the rheology of multi-droplet collision under magnetic field, the viscosity of ferrofluid droplet and surrounding fluid is set to $\eta _ { 1 } = 5 \times 1 0 ^ { - 2 }$ and $\eta _ { 2 } = 8 \times 1 0$ -2Pa s, respectively; their density is $\rho _ { 1 } = 1 5 8 0$ and $\rho _ { 2 } = 8 0 0 \mathrm { k g / m ^ { 3 } } ;$ ; and the magnetic susceptibility $\chi$ is 1.51 and 0, respectively. The surface tension coefficient between the two phases is $\sigma = 5 . 0 \mathrm { m N } / \mathrm { m } ^ { 2 }$ . The schematic diagram is shown in Fig. 11, the size of the computational domain is set to 12  6 mm, and the computational domain is discretized by $6 0 0 \times 3 0 0$ crossed rectangle cells. The interfacial thickness e spans across two grid cells. At the initial moment, two bubbles of radius $\mathrm { R } = 0 . 5$ mm are placed at (4, 3.69) and (8, 2.31), respectively. 

In this simulation, the farthest distance between two droplets is defined as the effective flow resistance height h. The inclination angle h represents the average of the inclination angles of the two droplets when they are not fused, and it also represents the inclination angle of the larger droplet after fusion. Re, Ca, and $\mathtt { B o } _ { \mathrm { M } }$ are chosen to be fixed as 0.003, 0.04, and 3.08, respectively. Figure 12 shows the change in relative viscosity throughout the process, which can be divided into three stages. During the first stage, which occurs between 0.1 and $0 . 6 s ,$ the droplets move smoothly in opposite directions. At this point, the effective flow resistance height h, relative viscosity $\eta _ { r } ,$ , and angle h of the two droplets remain almost unchanged. The second stage occurs 

![](images/977146f5d6d72b20c9ee6b27061e0ba3ecea5f2ec1372034d7cd633baa7fede4.jpg)


![](images/17f4a3f52335c1f876c1d82f098e632fb93ccec9a87e20d5805e322f00e574b8.jpg)



FIG. 9. (a) Relative viscosity $\eta _ { r } ,$ droplet orientation angle h, and effective flow resistance height h vs the magnetic field H(H for 1.0, 1.5, 2.0, 2.5, and 3.0 kA/m) at Re 0.375 and Ca  0.39. (b) Flow field schematic at $H = 2 . 0$ kA/m and $\gamma = 1 5 .$ . The left side shows velocity field and its streamlines, while the right side displays the magnetic field and its magnetic field lines. The white line in the middle represents the contour of phase field $\phi = 0 .$


![](images/3315b631b5beeb7c19b7b350c9affd7a2bb9b84a1d01622536081a1b3b8884bc.jpg)


![](images/959f8ffeb6bc0938782759d6d37ecdf62949b4a3d3b7ce6b63985c4005834785.jpg)



FIG. 10. (a) Relative viscosity ${ \ : \eta _ { r } , }$ , droplet orientation angle $\theta ,$ and effective flow resistance height h vs the shear rate $\gamma \left( \gamma \right)$ for $5 ,$ 10, 15, 20, 25, and 30). (b) Flow field schematic $\mathsf { a t } H = 2 . 0$ kA/m and $\gamma = 2 0 . \qquad $ The left side shows velocity field and its streamlines, while the right side displays the magnetic field and its magnetic field lines. The white line in the middle represents the contour of phase field /¼0.


![](images/4e9dcb6ec08a5345e95d7d62250ca956222f41b3fa342d6f29bf77af1925213b.jpg)



FIG. 11. Schematic diagram of droplet shear collision.


between 0.6 and $0 . 8 s ,$ during which the droplets move closer together and their long axes elongate away from the y-axis midline (x, 3) due to the field inhomogeneity of a single droplet affecting the other droplet. The effective flow resistance height of the two droplets increases, the angle h becomes larger, and the relative viscosity $\eta _ { r }$ increases. This phenomenon is consistent with the previously established conclusion that longer droplets tend to produce greater resistance to shear flow, and thus, greater effective viscosity as magnetic field strength increases. The third stage, occurring between 0.8 and 1.5 s, under the combined effect of surface tension and magnetic force, the droplets fuse and return to the equilibrium state, where the shape can be maintained, and the effective flow resistance height h, relative viscosity $\eta _ { r } ,$ , and angle h gradually return to values similar to those observed during the first stage. 

# D. Rosensweig’s instability of ferrofluid

This section utilizes a phase-field model to simulate the Rosensweig instability of a ferrofluid under the influence of a uniform magnetic field. Rosensweig’s instability in ferrofluids refers to small random deformations of the fluid surface due to environmental vibrations, thermal motions, and other perturbations, which can be amplified when a strong magnetic field is applied to the paramagnetic fluid along the normal direction. The magnetic flux is stronger at the top and gradually decreases toward the bottom, leading to continued deformation growth until a feedback loop is generated that modifies the total magnetic field. These results are in the emergence of regular peak and valley structures on the surface of the ferromagnetic fluid. The growth of these patterns will stop when magnetic force is in equilibrium with gravity and surface tension. 

![](images/47357aec8063c01dace58cbd47528ec195cd6384b8981b2c5884c0f91cf9c6a6.jpg)



FIG. 12. Time evolution of the relative viscosity $\eta _ { r }$ and droplet contact angle h at $\mathsf { R e } = 0 . 0 0 3$ and $\mathsf { C a } = 0 . 0 4$ . The insets in the figure illustrate the droplets states and effective flow resistance height h at corresponding times.


![](images/019f7b6108751f8cbdc2e9aa26ee271d74c0d97531c3913b5a3c0fee67ccb0c2.jpg)



FIG. 13. Schematic diagram of Rosensweig’s instability of magnetic liquids.


Figure 13 illustrates the setup for this study. A water-based ferrofluid droplet immiscible with alcohol is placed in the center of a cavity filled with alcohol. The densities of the ferrofluid and alcohol are $\rho _ { 1 } = 1 . 5 8 \times 1 0 ^ { 3 }$ and $\rho _ { 2 } = 0 . 8 \times 1 0 ^ { 3 } \mathrm { k g } / \mathrm { m } ^ { 3 }$ , respectively, while their viscosities are $\eta _ { 1 } = 3 2 \times 1 0 ^ { - 3 }$ and $\eta _ { 2 } = 0 . 8 \times 1 0 ^ { - 3 }$ Pas. The surface tension coefficient between the two phases is $\sigma = 3 . 8 9 \mathrm { m N } / \mathrm { m } ^ { 2 }$ . The initial value of the magnetic susceptibility is $\chi _ { 0 } = 2 . 2$ , and saturation magnetization $M _ { s }$ is 40 kA/m. As the saturation magnetization is much greater than the external magnetic field considered in this simulation, the magnetic susceptibility $\chi = \chi _ { 0 }$ . The computational domain is a rectangular container with length 20.5 mm and width 7 mm, where one-third of the container is filled with ferrofluid and the rest is filled with alcohol. The crossed rectangle grid size in the vessel is $6 0 0 \times 2 0 0 .$ . The interfacial thickness e spans across two grid cells. No-slip solid wall boundary conditions are imposed on the top and bottom boundaries of the computational domain, while the uniform magnetic field condition $\partial \psi / \partial y = H _ { 0 }$ is utilized for both the top and bottom boundaries. Periodic boundary conditions are used in the x direction. Wetted boundaries are also imposed on the left and right of the phase-field domain to provide initial conditions for the simulation. 

Figure 14 displays the shapes of the magnetic liquids in equilibrium at three different magnetic field strengths, with the reference experimental results3 shown on the left and the simulated results on the right. Under the influence of a uniform magnetic field $H _ { 0 }$ applied in the vertical direction, the ferrofluid’s shape reaches a delicate equilibrium over time, influenced by Maxwell stress, surface tension, and gravity. $\mathrm { A } s H _ { 0 }$ increases, the interfacial perturbation grows from a fine wavy shape to a comb-like structure with a single line peak. No hysteresis was observed in the simulations, and the magnetic susceptibility $\chi _ { 0 } < 2 . 5 4$ indicates that the hysteresis instability is a second-order jump.55 The critical values of magnetic field and wavelength are $H _ { c } = 4 . 7$ kA/m and $\lambda _ { c } = 2 . 9 2$ mm, respectively. The demagnetization factor (C) is negligible (in equation $\bar { H ^ { } } = H _ { 0 } \dot { - } \Gamma M _ { 0 } )$ since the magnetic field is only applied along the vertical axis. Therefore, the critical magnetization intensity Mc 10.3 kA/m. 

Krakov et $a l . ^ { 5 6 }$ suggested that the grid, magnetic field, and diffusion peak thickness have significant impacts on instability. In this study, an analytic expression for linear stability analysis is utilized to gain more insight into the dependence of Rosensweig’s instability on each parameter. This expression, provided by Cowley and Rosensweig pioneers,37 is valid only for small asymptotic values of the magnetic susceptibility $\chi _ { 0 } .$ By using this expression to construct physical parameters, we can study the effects of the grid, magnetic field, and diffusion peak on the instability. The critical magnetic field $H _ { c }$ and critical wavelength $\lambda _ { c }$ are determined by this expression 

$$
H _ {c} = \left(\frac {2}{\mu_ {0}} \frac {\mu_ {0} / \mu + 1}{(\mu_ {0} / \mu - 1) ^ {2}}\right) ^ {1 / 2} (\sigma g \Delta \rho) ^ {1 / 4}, \quad \lambda_ {c} = 2 \pi \left(\frac {\sigma}{g \Delta \rho}\right) ^ {1 / 2}. \tag {45}
$$

Here, $\Delta \rho = \rho _ { 1 } - \rho _ { 2 } ,$ and Equation (45) was used as an educated guess to predict the number of peaks in the vessel. We take a small magnetic susceptibility $\chi _ { 0 } = 0 . 9 , \ : \rho _ { 1 } = 1 . 5 8 \times 1 0 ^ { 3 } \mathrm { k g / m ^ { 3 } }$ and surface tension $\sigma = 6 . 0 \mathrm { m N } / \mathrm { m } ^ { 2 }$ for the ferrofluid. For the nonmagnetic fluid, $\chi _ { 0 } = 0 , \ : \rho _ { 2 } = 0 . 8 \times 1 0 ^ { 3 } \mathrm { k g / m ^ { 3 } }$ . Artificially assuming $\lambda _ { c } = 3 . 5$ mm, this gives a gravitational acceleration $\mathrm { g } = 2 \dot { 5 } \mathrm { m } / s ^ { 2 }$ and $H _ { c } { = } 1 0 . 8 \mathrm { k A / m }$ based on Eq. (45). In this simulation, the computational domain is chosen as a rectangular container with a length of 28 mm and a width of 7 mm, where one-third of the container is filled with ferrofluid and the rest with organic liquid. No-slip solid wall boundary conditions are imposed on the top and bottom boundaries of the computational domain, while the uniform magnetic field condition $( \partial \psi / \partial y = H _ { 0 } )$ is utilized for both the top and bottom boundaries. Periodic boundary conditions are used in the x direction. Wetted boundaries are also imposed on the left and right of the phase-field domain to provide initial conditions for the simulation. 

![](images/85ae7b3228a664808d9ec1402cbfdda0397db286ba92052af34050cc9678a66c.jpg)



FIG. 14. Comparison between experimental results and numerical simulations of the shape of ferromagnetic fluid liquid under the action of different magnetic fields, with the experimental results of the literature15 in the left column and the numerical simulation results in the right column.


![](images/063188409161fbb5eec5c2a55f8c8643c94f2017276f310fd7193fd70911187d.jpg)



FIG. 15. Characterization of ferrofluid morphology under steady-state conditions utilizing varied grid sizes.


![](images/665bcb4aaa247ea1d6cb4374f9d57f2ad4c3302b9d352de34a92fbe770b55fc2.jpg)



FIG. 17. Trend graph of peak heights and valley depths under different magnetic field intensities. The inserted illustration represents contour lines of the phase field at $\phi = 0 ,$ , depicting the profiles of peaks and valleys under different magnetic field conditions.


Rosensweig’s instability is a well-known phenomenon in ferrofluid, but its physical mechanism is complex as it involves a delicate balance between magnetic field force, gravity, and surface tension. Therefore, selecting an appropriate grid is crucial for the accuracy of numerical simulations. Four different crossed rectangle grid sizes $( \mathrm { M _ { 1 } } = 2 2 4 \times 5 6 , \mathrm { M _ { 2 } } = 3 3 6 \times 8 4 , \mathrm { M _ { 3 } } = 4 8 4 \times 1 1 2 .$ ; and $\dot { \mathbf { M } } _ { 4 } = 6 0 0$ 140) were simulated in this study. Under each mesh, the interfacial thickness e spans across two grid cells. Figure 15 presents a comparison of the peak shape at steady state for different grid sizes, under a magnetic field $H _ { 0 } = 1 5 \mathrm { k A / m }$ . The results indicate that grid size has a significant impact on the uniformity and wavelength of the peak. As the grid is refined, the wavelength $\lambda _ { c }$ decreases from 4 mm to a predicted value of 3.5 mm, which is in good agreement with the simulation results of the M3 grid (484 112) presented in (c) and (d). Therefore, to optimize simulation time and reduce cost, further simulations were conducted using the M3 grid. 

![](images/a495030accc9f357d42de301fbc373df77e89e30d42d356a8f9ff52370148dc1.jpg)



FIG. 16. Evolution of the peaks at different time steps. (Velocity field and its flow lines on the left, magnetic field and its magnetic field lines on the right, the white wave peaks represent the evolution of the phase field.).



TABLE I. Variation of wavelength with diffusion layer thickness.


<table><tr><td>Different diffusion layer thicknesses ε (mm)</td><td>Wavelength λc (mm)</td></tr><tr><td>1</td><td>5.6</td></tr><tr><td>0.75</td><td>4.9</td></tr><tr><td>0.5</td><td>4.3</td></tr><tr><td>0.25</td><td>3.8</td></tr><tr><td>0.125</td><td>3.5</td></tr></table>


TABLE II. CPU runtime (in seconds) for calculating the Rosensweig’s instability problem (with H0 15 kA/m) over 3000 time steps.


<table><tr><td>CPUcore count</td><td>Crossedrectangle grid</td><td>Totalcomputing time</td><td>Accelerationefficiency</td></tr><tr><td>1</td><td><eq>484 \times 112</eq></td><td>109 250.4(s)</td><td>1</td></tr><tr><td>2</td><td><eq>484 \times 112</eq></td><td>16 544.7(s)</td><td>3.3</td></tr><tr><td>4</td><td><eq>484 \times 112</eq></td><td>12 198.3(s)</td><td>2.23</td></tr><tr><td>8</td><td><eq>484 \times 112</eq></td><td>9829.2(s)</td><td>1.39</td></tr></table>

The property of ferrofluid being more easily magnetized than the surrounding fluid concentrates the magnetic field at the tips of the wave structure spikes, causing a reduction in magnetic energy and resulting in the movement of the ferrofluid spikes along the magnetic force lines.57 Figure 16 illustrates the spiking process under a magnetic field strength of $H _ { 0 } = 1 5 \mathrm { k A / m }$ . The left column shows the velocity field and its flow lines, while the right column shows the magnetic field and its magnetic field lines. The white wave crests represent the evolution of the phase field. The entire instability development process can be divided into three stages. The first stage is the starting spurring stage at t 0.02 s. Fine protrusions appear on both sides of the wetting boundary after applying the magnetic field, resulting in the initial spurring. The magnetic field concentrates with the plane fine protrusions and the magnetic energy at the peak increases, driving the local kinetic energy to increase and forming a vortex. According to Krakov et al.,56 the size of the vortex determines the wavelength of the wave. The second stage is the traction stage shown at $\mathrm { t } = 0 . 0 4 \ : s .$ As the first peak grows, it forms a barrier that blocks the flow of the internal liquid and allows a new vortex to form next to it. This forms a new fine protrusion in the plane, providing the conditions for the concentration of the magnetic field. The third stage is the stabilization stage. After the first and second stages, each peak’s kinetic energy gradually transfers to the next peak and stabilizes to reach the equilibrium state. 

According to the experiment of Krakov et al.,56 the wavelength k increases with an increase in the thickness of the diffusion front. In this section, the wavelengths are simulated for different diffusion layer thicknesses $\varepsilon _ { 1 } = 1 , \varepsilon _ { 2 } = 0 . 7 5 , \varepsilon _ { 3 } = 0 . 5 , \varepsilon _ { 4 } = 0 . 2 5 ,$ and $\varepsilon _ { 5 } = 0 . 1 2 5$ mm (Table I). 

To investigate the effect of magnetic field on instability, simulations were carried out for $H _ { 0 } = 1 \dot { 1 } _ { : }$ , 12, 14, 15, 16, and 18 kA/m. According to Fig. 17, it was found that magnetic field strength has almost no effect on the wavelength, but significantly influences the peak height and trough depth. Specifically, a larger magnetic field corresponds to a higher peak height and a deeper trough depth. However, as the magnetic field strength continues to increase, the peak height and trough depth approach each other. By fitting a LogNormal function to the peak height and trough depth data at different magnetic field strengths, it can be observed that the growth of these features stops beyond a certain threshold magnetic field strength. 

Finally, let us provide a brief overview of the simulation performance. Our simulations were conducted on a custom desktop computer equipped with an AMD Epyc 7532 32-core processor. Additionally, we utilized the built-in MPI module in FEniCS for parallel computation. In order to assess the parallel performance, we compared the time required for computing the Rosensweig instability problem over 3000 time steps (with $H _ { 0 } ^ { - } = 1 5 \mathrm { k A / m } )$ across different configurations, including single-core, dual-core, quad-core, and octacore setups. Table II presents the analysis of computation time results and acceleration ratios. 

# IV. CONCLUSION

In this work, we propose a phase-field-based finite element model for two-phase ferrofluid flows. This model combines the Cahn– Hilliard equation for the phase field, the Poisson equation for the magnetics, and the Navier–Stokes equations for fluid flow, which are strongly coupled. In this study, we present an efficient linear, decoupled numerical scheme. This led to, at each time step, a linear elliptic system for the phase function, a Poisson equation for the magnetic potential, a linear elliptic equation for the velocity, and a Poisson equation for the pressure. These schemes are efficient and easy to implement. This model is implemented via the FEniCS framework. 

Through an analysis of the deformation of ferrofluid droplets under various uniform magnetic fields, the accuracy and numerical stability of the proposed model have been verified. The equilibrium state of the droplets was then compared with the published experimental results, revealing a strong agreement between the two in terms of aspect ratio and curvature shape. 

Furthermore, to verify the ability of this method to simulate three-dimensional multiphase flow problems in ferrofluids with highdensity ratios, the merging of two bubbles in a ferrofluid with a density ratio of 850 and a viscosity ratio of 280 was simulated. A detailed analysis of the different phenomena and reasons for each stage of the merging process is presented. The results also show that the smaller the magnetic Bond number, the less the droplet oscillation. 

We examine the rheological properties of single and multiple dilute ferrofluid droplet suspensions under shear flow conditions. The notion of the “effective flow resistance height” is introduced, proposing mechanisms by which effective modifications to the flow resistance height of ferrofluid droplets can efficiently regulate the rheological behavior of the emulsion. 

Finally, this paper investigates the interfacial instability of ferrofluids under the action of a uniform magnetic field. The wavelength decreases with grid refinement, increases with increasing interface thickness, and exhibits insensitivity to the magnetic field response. The peak height and valley depth are significantly affected by the magnetic field, but once saturation is reached, the peak no longer undergoes further growth. Furthermore, we conduct a detailed analysis of the formation and development of peaks within the system. Our findings reveal that this phenomenon is a coupled process involving both magnetic energy concentration and the generation and transfer of magnetic flow vortices. 

# ACKNOWLEDGMENTS

This work was supported by the National Natural Science Foundation of China (Grant Nos. 12102228, 52275201, and 12172039). 

# AUTHOR DECLARATIONS

# Conflict of Interest

The authors have no conflicts to disclose. 

# Author Contributions

Conceptualization (equal); Data curation (equal); Pengfei Yuan:Formal analysis (equal); Investigation (equal); Methodology (equal); Project administration (equal); Resources (equal); Software (equal); Validation (equal); Visualization (equal); Writing – original draft (equal); Writing – review & editing (equal).   Data cura-Qianxi Cheng:tion (equal); Formal analysis (equal); Visualization (equal); Writing – original draft (equal); Writing – review & editing (equal). Yang Hu:Data curation (equal); Formal analysis (equal); Investigation (equal); Methodology (equal); Software (equal); Supervision (equal); Visualization (equal); Writing – original draft (equal). Qiang He:Conceptualization (equal); Formal analysis (equal); Funding acquisition (equal); Investigation (equal); Methodology (equal); Software (equal); Supervision (equal); Visualization (equal); Writing – original draft (equal); Writing – review & editing (equal). Weifeng Huang:Formal analysis (equal); Methodology (equal); Project administration (equal); Supervision (equal).  Formal analysis (equal); Decai Li:Methodology (equal); Supervision (equal); Visualization (equal). 

# DATA AVAILABILITY

The data that support the findings of this study are available from the corresponding author upon reasonable request. 

# NOMENCLATURE

# SYMBOL

a Half-short axis of the ferrofluid droplet (m) 

Magnetic induction (T) 

$\mathtt { B o } _ { \mathrm { M } }$ Magnetic Bond number 

b Half-length axis of the ferrofluid droplet (m) 

$\mathrm { C _ { a } }$ Capillary numberc0Safety factor 

E System’s total kinetic energy (J) 

$\mathbf { f } _ { b }$ Body force (N) 

$\mathbf { f } _ { m }$ Magnetic force (N) 

$\mathbf { f } _ { s }$ Surface tension (N) 

$G$ Chemical potential (J) 

g Gravitational acceleration (m/s2) 

$\mathbf { H }$ Magnetic field (A/m) 

Hc Critical values of magnetic field (A/m) 

h Effective flow resistance height (m) 

Magnetization intensity (A/m) 

$M _ { s }$ Saturation magnetization intensity (A/m) 

$M ( \phi )$ Mobility coefficient 

$M _ { 0 }$ Initial phase-field mobility 

$\mathrm { L }$ Computational domain length (m) 

$M _ { t }$ Total mass of the system 

$\mathrm { M _ { c } }$ Critical magnetization intensity (A/m) 

$\mathrm { O _ { h } }$ Ohnesorge number 

$\boldsymbol { p }$ Pressure (Pa) 

$\mathrm { R }$ Radius of ferrofluid droplet (m) 

$\mathrm { R _ { e } }$ Reynolds number 

$S$ Average shear stress (Pa) 

t Time (s) 

U Characteristic velocity (m/s) 

Velocity field (m/s) 

uutop $u _ { \mathrm { t o p } }$ Flat plate moving velocity (m/s) 

a Viscosity ratio 

C Demagnetization factor (kJ/m3) 

c Shear rate (s-1 ) 

gr Suspension’s relative viscosity 

$\eta _ { s }$ Effective viscosity of the suspension (Pa s) 

g1 Viscosity of ferrofluid (Pa s) 

$\eta _ { 2 }$ Viscosity of nonmagnetic fluid (Pa s) 

e Interfacial thickness (m) 

h Orientation angle (
) 

hc Predetermined contact angle (
) 

$\lambda _ { c }$ Critical values of wavelength (m) 

k Wavelength (m) 

l Permeability (N/A2) 

l0 Vacuum permeability (N/A2) 

q1 Density of ferrofluid (kg=m3) 

q2 Density of nonmagnetic fluid (kg=m3) 

r Surface tension (N/m2) 

v Magnetic susceptibility 

v0 Initial magnetic susceptibility 

/ Phase field 

w Scalar potential (A) 

# REFERENCES



1 M. D. Cowley, Magnetic Fluids: Engineering Applications, edited by B. M. Berkovsky, V. F. Medvedev and M. S. Krakov (Oxford University Press, 1993), p. 243; V. Bashtovoy and M. Berkovsky, Magnetic Fluids and Applications Handbook, edited by B. M. Berkovsky (Begell House Inc., 1996). 





2 I. Torres-Díaz and C. Rinaldi, “Recent progress in ferrofluids research: Novel applications of magnetically controllable and tunable fluids,” Soft Matter (43), 8584–8602 (2014). 





103 C. Flament, S. Lacis, J. C. Bacri, A. Cebers, S. Neveu, and R. Perzynski, “Measurements of ferrofluid surface tension in confined geometry,” Phys. Rev. E , 4801–4806 (1996). 





534 Q. He, W. Huang, J. Xu, Y. Hu, and D. Li, “A hybrid immersed interface and phase-field-based lattice Boltzmann method for multiphase ferrofluid flow,” Comput. Fluids , 105821 (2023). 





2555 O. Lavrova, G. Matthies, T. Mitkova, V. Polevikov, and L. Tobiska, “Numerical treatment of free surface problems in ferrohydrodynamics,” J. Phys. , S2657 (2006). 





6 C. Gollwitzer, G. Matthies, R. Richter, I. Rehberg, and L. Tobiska, “The surface topography of a magnetic fluid: a quantitative comparison between experiment and numerical simulation,” J. Fluid Mech. , 455–474 (2007). 





5717 M. S. Korlie, A. Mukherjee, B. G. Nita, J. G. Stevens, A. D. Trubatch, and P. Yecko, “Modeling bubbles and droplets in magnetic fluids,” J. Phys. , 204143 (2008). 





8 G.-P. Zhu, N.-T. Nguyen, R. V. Ramanujan, and X.-Y. Huang, “Nonlinear deformation of a ferrofluid droplet in a uniform magnetic field,” Langmuir , 14834–14841 (2011). 





9 Y. Cao and Z. J. Ding, “Formation of hexagonal pattern of ferrofluid in magnetic field,” J. Magn. Magn. Mater. , 93–99 (2014). 





35510M. Hamid, M. Usman, Z. H. Khan, R. Ahmad, and W. Wang, “Dual solutions and stability analysis of flow and heat transfer of Casson fluid over a stretching sheet,” Phys. Lett. A , 2400–2408 (2019). 





38311Y. Hu, D. Li, and X. Niu, “Phase-field-based lattice Boltzmann model for multiphase ferrofluid flows,” Phys. Rev. E , 033301 (2018). 





9812Y. Li, X.-D. Niu, A. Khan, D.-C. Li, and H. Yamaguchi, “A numerical investigation of dynamics of bubbly flow in a ferrofluid by a self-correcting procedurebased lattice Boltzmann flux solver,” Phys. Fluids , 082107 (2019). 





3113S.-T. Zhang, X.-D. Niu, Q.-P. Li, A. Khan, Y. Hu, and D.-C. Li, “A numerical investigation on the deformation of ferrofluid droplets,” Phys. Fluids , 012102 (2023). 





14W. K. Lee, R. Scardovelli, A. D. Trubatch, and P. Yecko, “Numerical, experimental, and theoretical investigation of bubble aggregation and deformation in magnetic fluids,” Phys. Rev. E , 016302 (2010). 





8215S. Afkhami, A. J. Tyler, Y. Renardy, M. Renardy, T. G. S. Pierre, R. C. Woodward, and J. S. Riffle, “Deformation of a hydrophobic ferrofluid droplet suspended in a viscous medium under uniform magnetic fields,” J. Fluid Mech. , 358–384 (2010). 





66316D. Shi, Q. Bi, and R. Zhou, “Numerical simulation of a falling ferrofluid droplet in a uniform magnetic field by the VOSET method,” Numer. Heat Transfer, Part A , 144–164 (2014). 





6617M. Habera and J. Hron, “Modelling of a free-surface ferrofluid flow,” J. Magn. Magn. Mater. , 157–160 (2017). 





43118X. Ni, B. Zhu, B. Wang, and B. Chen, “A level-set method for magnetic substance simulation,” ACM Trans. Graphics , 29:1–29:13 (2020). 





39 19R. Maroofiazar, M. Daryani, and A. R. Vakhshouri, “Numerical investigation of ferrofluid sloshing by applying MHD magnetic field: Using level set method,” Eur. J. Comput. Mech. , 351–372 (2019). 





28 20R. H. Nochetto, A. J. Salgado, and I. Tomas, “A diffuse interface model for twophase ferrofluid flows,” Comput. Methods Appl. Mech. Eng. , 497–531 (2016). 





21Y. Hu, D. Li, and Q. He, “Generalized conservative phase field model and its lattice Boltzmann scheme for multicomponent multiphase flows,” Int. J. Multiphase Flow , 103432 (2020). 





132 22Y. Hu, D. Li, X. Niu, and S. Shu, “A diffuse interface lattice Boltzmann model for thermocapillary flows with large density ratio and thermophysical parameters contrasts,” Int. J. Heat Mass Transfer , 809–824 (2019). 





13823A. Ghaderi, M. H. Kayhani, and M. Nazari, “Simulation of the co-axial ferrofluid droplets interaction under uniform magnetic field,” Therm. Sci. , 1027 (2019). 





24X. Li, Z.-Q. Dong, P. Yu, X.-D. Niu, L.-P. Wang, D.-C. Li, and H. Yamaguchi, “Numerical investigation of magnetic multiphase flows by the fractional-stepbased multiphase lattice Boltzmann method,” Phys. Fluids , 083309 (2020). 





3225X.-D. Niu, A. Khan, Y. Ouyang, M.-F. Chen, D.-C. Li, and H. Yamaguchi, “A simplified phase-field lattice Boltzmann method with a self-corrected magnetic field for the evolution of spike structures in ferrofluids,” Appl. Math. Comput. , 127503 (2023). 





43626M. Majidi, M. A. Bijarchi, A. G. Arani, M. H. Rahimian, and M. B. Shafii, “Magnetic field-induced control of a compound ferrofluid droplet deformation and breakup in shear flow using a hybrid lattice Boltzmann-finite difference method,” Int. J. Multiphase Flow , 103846 (2022). 





14627F. Bassi and S. Rebay, “A high-order accurate discontinuous finite element method for the numerical solution of the compressible Navier–Stokes equations,” J. Comput. Phys. , 267–279 (1997). 





13128A. Logg, K.-A. Mardal, and G. N. Wells, Automated Solution of Differential Equations by the Finite Element Method: The FEniCS Book (Springer, 2012). 





29M. Mortensen and K. Valen-Sendstad, “Oasis: A high-level/high-performance open source Navier–Stokes solver,” Comput. Phys. Commun. , 177–188 (2015). 





30G. Mitscha-Baude, A. Buttinger-Kreuzhuber, G. Tulzer, and C. Heitzinger, “Adaptive and iterative methods for simulations of nanopores with the PNP– Stokes equations,” J. Comput. Phys. , 452–476 (2017). 





33831A. Bolet, G. Linga, and J. Mathiesen, “Electrohydrodynamic channeling effects in narrow fractures and pores,” Phys. Rev. E , 043114 (2018). 





9732Y. Hu, Q. He, D. Li, Y. Li, and X. Niu, “On the total mass conservation and the volume preservation in the diffuse interface method,” Comput. Fluids , 104291 (2019). 





33G.-D. Zhang, X. He, and X. Yang, “Reformulated weak formulation and efficient fully discrete finite element method for a two-phase ferrohydrodynamics Shliomis model,” SIAM J. Sci. Comput. , B253–B282 (2023). 





45   34F. Magaletti, F. Picano, M. Chinappi, L. Marino, and C. M. Casciola, “The sharp-interface limit of the Cahn–Hilliard/Navier–Stokes model for binary fluids,” J. Fluid Mech. , 95–126 (2013). 





71435E. Campillo-Funollet, G. Grun, and F. Klingbeil, € “On modeling and simulation of electrokinetic phenomena in two-phase flow with general mass densities,” SIAM J. Appl. Math. , 1899–1925 (2012). 





7236A. O. Ivanov and O. B. Kuznetsova, “Nonmonotonic field-dependent magnetic permeability of a paramagnetic ferrofluid emulsion,” Phys. Rev. E , 041405 (2012). 





37M. D. Cowley and R. E. Rosensweig, “The interfacial stability of a ferromagnetic fluid,” J. Fluid Mech. , 671–688 (1967). 





3038H. Ding, P. D. M. Spelt, and C. Shu, “Diffuse interface model for incompressible two-phase flows with large density ratios,” J. Comput. Phys. , 2078– 2095 (2007). 





39D. Jacqmin, “Contact-line dynamics of a diffuse fluid interface,” J. Fluid Mech. , 57–88 (2000). 





40240H. P. Langtangen, K.-A. Mardal, and R. Winther, “Numerical methods for incompressible viscous flow,” Adv. Water Resour. , 1125–1146 (2002). 





2541S. Nagrath, L. V. Sequist, S. Maheswaran, D. W. Bell, D. Irimia, L. Ulkus, M. R. Smith, E. L. Kwak, S. Digumarthy, A. Muzikansky, P. Ryan, U. J. Balis, R. G. Tompkins, D. A. Haber, and M. Toner, “Isolation of rare circulating tumour cells in cancer patients by microchip technology,” Nature , 1235–1239 (2007). 





42J. Pipper, Y. Zhang, P. Neuzil, and T.-M. Hsieh, “Clockwork PCR including sample preparation,” Angew. Chem., Int. Ed. , 3900–3904 (2008). 





4743N. Pamme and C. Wilhelm, “Continuous sorting of magnetic cells via on-chip free-flow magnetophoresis,” Lab Chip (8), 974–980 (2006). 





644L. Hajba and A. Guttman, “Circulating tumor-cell detection and capture using microfluidic devices,” TrAC Trends Anal. Chem. , 9–16 (2014). 





5945M. A. M. Gijs, “Magnetic bead handling on-chip: New opportunities for analytical applications,” Microfluid. Nanofluid. , 22–40 (2004). 





1 46H.-H. Chen and D. Gao, “Particle enrichment employing grooved microfluidic channels,” Appl. Phys. Lett. , 173502 (2008). 





92   47P. K. Yuen, L. J. Kricka, P. M. Fortina, N. J. Panaro, T. Sakazume, and P. Wilding, “Microchip module for blood sample preparation and nucleic acid amplification reactions,” Genome Res. (3), 405–412 (2001). 





1148L. H. P. Cunha, I. R. Siqueira, T. F. Oliveira, and H. D. Ceniceros, “Fieldinduced control of ferrofluid emulsion rheology and droplet break-up in shear flows,” Phys. Fluids , 122110 (2018). 





30   49H. W. Zheng, C. Shu, and Y. T. Chew, “Lattice Boltzmann interface capturing method for incompressible flows,” Phys. Rev. E , 056705 (2005). 





72 50P. Chen, M. P. Dudukovi-c, and J. Sanyal, “Three-dimensional simulation of bubble column flows with bubble coalescence and breakup,” AIChE J. , 696– 712 (2005). 





51L. H. P. Cunha, I. R. de Siqueira, F. R. Cunha, and T. F. Oliveira, “Effects of external magnetic fields on the rheology and magnetization of dilute emulsions of ferrofluid droplets in shear flows,” Phys. Fluids , 073306 (2020). 





32 52W. Peng, Y. Hu, D. Li, and Q. He, “Full-scale simulation of the fluid-particle interaction under magnetic field based on IIM–IBM–LBM coupling method,” Front. Mater. , 932854 (2022). 





9 53A. Xu, L. Shi, and T. S. Zhao, “Lattice Boltzmann simulation of shear viscosity of suspensions containing porous particles,” Int. J. Heat Mass Transfer , 969–976 (2018). 





54G. Ghigliotti, T. Biben, and C. Misbah, “Rheology of a dilute two-dimensional suspension of vesicles,” J. Fluid Mech. , 489–518 (2010). 





65355V. M. Zaitsev and M. I. Shliomis, “Nature of the instability of the interface between two liquids in a constant field,” Soviet Phys. Dokl. , 1001 (1970). 





56M. S. Krakov, A. R. Zakinyan, and A. A. Zakinyan, “Instability of the miscible magnetic/non-magnetic fluid interface,” J. Fluid Mech. , A30 (2021). 





91357D. Andelman and R. E. Rosensweig, “The phenomenology of modulated phases: From magnetic solids and fluids to organic films and polymers,” in Polymers, Liquids and Colloids in Electric Fields (World Scientific, 2009), pp. 1–56. 

