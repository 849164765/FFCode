# Phase-field model for the Rayleigh–Taylor instability of immiscible fluids

By A N T O N I O C E L A N I1 , A N D R E A M A Z Z I N O2 , P A O L O M U R A T O R E - G I N A N N E S C H I3 and L A R A V O Z E L L A2,3 



1Institut Pasteur, CNRS, URA 2171, 25 Rue du docteur Roux, 75015 Paris, France 





2Department of Physics - University of Genova, and CNISM & INFN - Genova Section, via Dodecaneso 33, 16146 Genova, Italy 





3Department of Mathematics and Statistics - University of Helsinki, P.O. Box 4, 00014 Helsinki, Finland 



(Received 31 October 2018) 

The Rayleigh–Taylor instability of two immiscible fluids in the limit of small Atwood numbers is studied by means of a phase-field description. In this method the sharp fluid interface is replaced by a thin, yet finite, transition layer where the interfacial forces vary smoothly. This is achieved by introducing an order parameter (the phase field) whose variation is continuous across the interfacial layers and is uniform in the bulk region. The phase field model obeys a Cahn–Hilliard equation and is two-way coupled to the standard Navier–Stokes equations. Starting from this system of equations we have first performed a linear analysis from which we have analytically rederived the known gravity-capillary dispersion relation in the limit of vanishing mixing energy density and capillary width. We have performed numerical simulations and identified a region of parameters in which the known properties of the linear phase (both stable and unstable) are reproduced in a very accurate way. This has been done both in the case of negligible viscosity and in the case of nonzero viscosity. In the latter situation only upper and lower bounds for the perturbation growth-rate are known. Finally, we have also investigated the weaklynonlinear stage of the perturbation evolution and identified a regime characterized by a constant terminal velocity of bubbles/spikes. The measured value of the terminal velocity is in perfect agreement with available theoretical prediction. The phase-field approach thus appears to be a valuable tecnhique for the dynamical description of the stages where hydrodynamic turbulence and wave-turbulence enter into play. 

# 1. Introduction

The Rayleigh-Taylor (RT) instability is a fluid-mixing mechanism occurring when a heavy, denser, fluid is pushed into a lighter one. For a fluid in a gravitational field, such a mechanism was first discovered by Lord Rayleigh in the 1880s (Rayleigh 1883) and later applied to all accelerated fluids by Sir Geoffrey Taylor in 1950 (Taylor 1950). The relevance of this mixing mechanism embraces many different phenomena occurring in completely different contexts. We just mention, among the many, astrophysical supernova explosions and geophysical formations like salt domes and volcanic islands (Di Prima & Swinney 1981; Dimonte & Schneider 2000), continental magmatism caused by lithospheric gravitational instability (Lee, Rudnick & Brimhall Jr. 2001; Ducea & Saleeby 

1998), inertial confinement fusion (Cook & Zhou 2002) and cloud formation in atmospheric sciences (Schultz et al. 2006). 

Back to classical fluids applications, RT instability is the first step eventually leading to a fully developed turbulent regime. A deeper understanding of the mechanism of flows driven by RT instability thus would shed light on the many processes that underpin fully developed turbulence. 

The difficulty inherent in sustaining an unstable density stratification has challenged experimentalists for over half a century. Several innovative approaches have been recently developed (see e.g., Ramaprabhu & Andrews 2004). 

With the advent of supercomputers, high-resolution numerical simulations of RT at high Reynolds numbers have become a reality. However, simulations using many different benchmark codes and experiments disagree already on apparently innocent observables like, for instance, the value of the growth constant, α, associated to the spread of the turbulent mixing zone (see, e.g., Di Prima & Swinney 1981). The differences can be as high as 100%. 

Despite the long history of RT turbulence, a consistent phenomenological theory has been presented only very recently by Chertkov (2003) for the miscible case. The theoretical predictions by Chertkov have been verified by Celani, Mazzino & Vozella (2006) exploiting numerical simulations in two spatial dimensions. For the three-dimensional miscible case we refer, e.g, to Young et al. (2001). 

In many of the aforementioned situations where the RT instability has an important role, the two fluids are immiscible owing to a non negligible surface tension. At level of linear analysis the role played by a non zero surface tension was addressed by Chandrasekhar (1961). The successive dynamics falling in a turbulent regime has been recently analyzed by Chertkov, Kolokolov & Lebedev (2005). Using a phenomenological approach, the authors suggest the existence of a Kolmogorov cascade between the integral scale and a time-dependent scale related to the typical drop size. Below the latter scale, associated to an emulsion-like region, a wave energy cascade takes in. This is mediated by weakly interacting capillary waves propagating on top of the drop surface. Eventually, the energy is dissipated by viscous forces. 

RT instability and RT turbulence of immiscible fluids thus appear richer than the corresponding miscible situations. The existence of two different cascades poses a serious challenge to numerical investigations of the immiscible RT problems. The emulsion-like phase indeed occurs at very small scales and the energy transfer takes place on the interfaces. These are geometrical objects close to singularities and thus difficult to describe appropriately in a numerical scheme. Accuracy and efficiency are thus fundamental requirements to reproduce the correct statistical features characterizing immiscible RT turbulence. 

Our aim here is to perform a first step along this direction by focusing on direct numerical simulations of immiscible RT instability. The numerical strategy we exploit here is known as phase-field model (Bray 2002; Cahn & Hilliard 1958; Badalassi, Ceniceros & Banerjee 2003; Ding, Spelt & Shu 2007). The main idea of the method is to treat the interface between two immiscible fluids as a thin mixing layer across which physical properties vary steeply but continuously. The evolution of the mixing layer is ruled by an order parameter (the phase field) that obeys a Cahn-Hilliard equation (Cahn & Hilliard 1958). The method permits to avoid a direct tracking of the interface and easily produces the correct interfacial tension from the mixing layer free energy. 

We present here an accurate numerical study that validates the phase–field approach by testing known results of immiscible RT instability both at level of linear and weakly nonlinear analysis. From our results, it turns out that this strategy is a valuable option for a quantitative treatment of the turbulent regime characterized by the interplay between hydrodynamic and interface degrees of freedom. 

![](images/caef344f616b6e2fa3ae4a117fd7cd842b9e0de0c40ea7de2cdf7403865e4860.jpg)



Figure 1. Fluids configuration corresponding to a heavier fluid of density $\rho _ { 2 }$ placed above a lighter one of density $\rho _ { 1 } < \rho _ { 2 }$ .


The paper is organized as follows. In Sec. 2 we introduce the Rayleigh–Taylor problem and discuss the related phase-field approach. A detailed analysis of the energy balance between purely hydrodynamic degrees of freedom and interface degrees of freedom is presented. Finally, the dispersion relation for gravity-capillary waves is obtained by analytical calculations starting from the phase-field equation coupled to the Navier–Stokes equations. 

In Sec. 3 the results from the direct numerical simulations are presented and compared with known results for the linear analysis. We focus both on the case of zero viscosity and on that of negligible viscosity. Both stable and unstable configurations are considered. Finally, the weakly nonlinear regime is considered and the resulting terminal velocity of bubbles/spikes compared with existing theoretical predictions. 

Sec. 4 is devoted to some conclusions and perspectives. 

# 2. System configuration and phase-field model

Our system consists of two immiscible, incompressible fluids (labeled by 1 and 2) having different densities, $\rho _ { 1 }$ and $\rho _ { 2 } > \rho _ { 1 }$ , with the denser fluid placed, e.g., above the less dense one (see Fig. 1). In the absence of gravity, this flow configuration is stable. In presence of the gravitational force, surface tension may be able to keep the system in equilibrium, provided the density contrast is not too large. 

Let us start by describe the equilibrium configuration and then pass to the evolution (RT instability) that occurs when a perturbation is imposed to the interface separating the two fluids. 

# 2.1. Equilibrium state

Let us consider an equilibrium state where fluid 1 is placed below fluid 2 and they are separated by a sharp interface. The fact that the interface is sharp (i.e. a discontinuity in the fluid properties) poses a serious challenge to numerical simulations. Indeed, for sharp interfaces, the evolution equations are obtained by following fluid 1 and 2 separately with the appropriate boundary condition at interface (see, for instance, Smolianski, Haario & Luukka 2005; Sethian 1999). Other approaches follow the interface alone. In this latter case, the movement of the interface is naturally amenable to a Lagrangian description, while the bulk flow is conventionally solved in an Eulerian framework. These approaches employ a mesh that has grid points on the interfaces and deforms according to the flow. A major shortcoming of these approaches is in that they cannot handle properly topological changes such as breakup, coalescence and reconnections (see Yue et al. 2004, and references therein). In this respect, the phase–field method is, by far and large, more effective, at the expense of a larger number of grid points required. 

The idea of the phase–field method is to replace the sharp interface with a diffuse one in such a way that the numerical computation of interface movement and deformation can be carried out on fixed grids (Anderson, McFadden & Wheeler 1998; Jacqmin 1999). More quantitatively, this amounts to assigning to the system a Ginzburg–Landau free energy, F, espressed in term of the order parameter $\phi$ as (Cahn & Hilliard 1958; Bray 2002; Yue et al. 2004): 

$$
\mathcal {F} [ \phi ] = \int_ {\Omega} \frac {\Lambda}{2} | \boldsymbol {\partial} \phi (\boldsymbol {x}) | ^ {2} + \frac {\Lambda}{4 \epsilon^ {2}} (\phi^ {2} - 1) ^ {2} d \boldsymbol {x}, \tag {2.1}
$$

where Ω is the region of space occupied by the system, Λ is a mixing energy density and ǫ is the capillary width, representative of the interface thickness. The order parameter φ is a field which serves to identify fluid 1 and 2. We assume $\phi = 1$ in the region occupied by fluid 1 and $\phi = - 1$ in those where fluid 2 is present. 

The equilibrium state is the minimizer of the free energy . The mechanism which keeps the system in this configuration is the competition between two effects due to the two addends in (2.1). The first term favours a perfect mixing (i.e. $\Lambda | \partial \phi | ^ { 2 } / 2 = 0$ in F , this term being the interface energy contribution) whereas the second one one drives the system towards demixing (the associated term in F, the bulk contribution, has indeed a minimum for $\phi = \pm 1 )$ ). The nontrivial final equilibrium state is just the results of this competition. More quantitatively, the final state is obtained by minimizing the free-energy functional with respect to variations of the function φ, i.e., solving: 

$$
\mu \equiv \delta \mathcal {F} / \delta \phi = 0 \Leftrightarrow - \partial^ {2} \phi + \frac {\phi^ {3} - \phi}{\epsilon^ {2}} = 0, \tag {2.2}
$$

where $\mu$ is the so-called chemical potential (see, for instance, Cahn & Hilliard 1958; Bray 2002; Yue et al. 2004). If one considers an one-dimensional interface, varying along the gravitational direction y, one easily finds the solution of Eq. (2.2) as (Cahn & Hilliard 1958; Bray 2002; Yue et al. 2004): 

$$
\phi (y) = \pm \tanh \left(\frac {y}{\sqrt {2} \epsilon}\right). \tag {2.3}
$$

This solution exists and is stable in all dimensions although the decay rate of perturbations depends upon the dimensionality (Korvola, Kupiainen & Taskinen 2005). From (2.2) one immediately realizes that the sharp-interface limit is obtained for ǫ → 0: in this case tanh $\left( y / ( \sqrt { 2 } \epsilon ) \right) \to \mathrm { s i g n } ( y )$ . Moreover, the surface tension σ is equal to the integral of the free-energy density along the interface (see, for example, Landau & Lifshitz 2000). For a plane interface, this integral yields (Cahn & Hilliard 1958; Bray 2002; Yue et al. 2004): 

$$
\sigma = \frac {2 \sqrt {2}}{3} \frac {\Lambda}{\epsilon}. \tag {2.4}
$$

It is now easy to verify how the sharp interface limit is obtained: it suffices to take the limits Λ and ǫ to zero keeping σ fixed to the value prescribed by surface tension (Liu & Shen 2003). 

# 2.2. Perturbation evolution

Let us now suppose to impose a small perturbation on the (finite thickness) interface separating the two fluids. Such perturbation will displace the phase field from the previous equilibrium configuration, which minimized the free-energy $\mathcal { F }$ , to a new configuration for which in general, $\mu \neq 0$ . The system will react so as to try to reach again an equilibrium configuration. In formulae: 

$$
\frac {\partial \phi}{\partial t} + \pmb {v} \cdot \pmb {\partial} \phi = \gamma \pmb {\partial} ^ {2} \mu = \gamma \Lambda \pmb {\partial} ^ {2} \left[ - \partial^ {2} \phi + \frac {(\phi^ {3} - \phi)}{\epsilon^ {2}} \right] \quad , (2. 5)
$$

γ being the so-called mobility (see, for instance, Bray 2002; Yue et al. 2004). Notice the presence of the Laplacian operator in front of $\mu .$ . Notice that the mass of each fluid is conserved, as imposed by the physics of the problem under consideration. 

The dynamics of the velocity field is governed by the usual Boussinesq Navier-Stokes equations (Kundu & Cohen 2001) plus an additional stress contribution arising at the interface where the effect of surface tension enters into play (Bray 2002; Yue et al. 2004; Berti et al. 2005). The equations of motion are: 

$$
\left(\partial_ {t} v _ {\alpha} + \boldsymbol {v} \cdot \boldsymbol {\partial} v _ {\alpha}\right) = - \frac {\partial_ {\alpha} p}{\rho_ {o}} + \nu \partial^ {2} v _ {\alpha} - \frac {\phi}{\rho_ {o}} \partial_ {\alpha} \frac {\delta \mathcal {F}}{\delta \phi} + \frac {\rho^ {\prime}}{\rho_ {o}} g _ {\alpha} \tag {2.6}
$$

$$
\boldsymbol {\partial} \cdot \boldsymbol {v} = 0. \tag {2.7}
$$

In the first equation $\rho _ { o } = ( \rho _ { 1 } + \rho _ { 2 } ) / 2$ and ν is the kinematic viscosity. The quantity $- \phi \pmb { \partial } ( \delta \mathcal { F } / \delta \phi ) / \rho _ { o }$ is the coupling term that accounts for capillary forces. It is easy to verify that it can be rewritten as $- \Lambda \left( \partial ^ { 2 } \phi \partial \phi \right) / \rho _ { o }$ plus a gradient term which can be absorbed into the pressure term. Finally, $\rho ^ { \prime } g _ { \alpha } / \rho _ { o }$ is the buoyancy contribution, $\rho ^ { \prime }$ being the deviation of the actual density, $\rho ,$ from the mean density $\rho _ { o } \colon$ 

$$
\rho^ {\prime} = \rho - \rho_ {o}.
$$

The buoyancy contribution can be rewritten in terms of $\rho _ { 1 } , \rho _ { 2 }$ and $\phi$ as: 

$$
\begin{array}{l} \frac {\rho^ {\prime}}{\rho_ {o}} g _ {\alpha} = \frac {\rho - \rho_ {o}}{\rho_ {o}} g _ {\alpha} = \\ = \frac {\rho_ {1} \left(\frac {1 + \phi}{2}\right) + \rho_ {2} \left(\frac {1 - \phi}{2}\right) - \rho_ {o}}{\rho_ {o}} g _ {\alpha} \\ = - \mathcal {A} \phi g _ {\alpha} \tag {2.8} \\ \end{array}
$$

where $\mathcal { A } \equiv ( \rho _ { 2 } - \rho _ { 1 } ) / ( \rho _ { 2 } + \rho _ { 1 } )$ is the Atwood number. 

# 2.3. Energetics

Let us define the kinetic energy (per unit volume), $E _ { K }$ , and the potential energy (per unit volume), $E _ { P }$ , for our system ruled by Eqs. (2.5), (2.6) and (2.7). 

By definition of potential energy, we have: 

$$
\begin{array}{l} E _ {P} = \frac {1}{\Omega} \int \int d x d y \rho_ {2} g y \frac {1 - \phi}{2} + \frac {1}{\Omega} \int \int d x d y \rho_ {1} g y \frac {1 + \phi}{2} + E _ {P} ^ {o} = \\ = - \frac {1}{2} \langle y \phi \rangle (\rho_ {2} - \rho_ {1}) g = - \rho_ {o} \mathcal {A} g \langle y \phi \rangle , \tag {2.9} \\ \end{array}
$$

Ω being the total volume occupied by the fluids and brackets, $\langle \cdots \rangle$ , denote spatial averages. In Eq. (2.9) the constant $E _ { P } ^ { o }$ is chosen such to set the potential energy to zero for 

vanishing Atwood number. 

In a similar way, one can define the kinetic energy per unit volume as: 

$$
\begin{array}{l} E _ {K} = \frac {1}{\Omega} \int \int d x   d y   \rho_ {2}   \frac {1 - \phi}{2}   \frac {\boldsymbol {v} ^ {2}}{2} + \frac {1}{\Omega} \int \int d x   d y   \rho_ {1}   \frac {1 + \phi}{2}   \frac {\boldsymbol {v} ^ {2}}{2} = \\ = \rho_ {2} \left\langle \left(\frac {1 - \phi}{2}\right) \frac {\boldsymbol {v} ^ {2}}{2} \right\rangle + \rho_ {1} \left\langle \left(\frac {1 + \phi}{2}\right) \frac {\boldsymbol {v} ^ {2}}{2} \right\rangle = \\ = \rho_ {o} \left\langle \frac {\boldsymbol {v} ^ {2}}{2} \right\rangle - \rho_ {o} \mathcal {A} \left\langle \phi \frac {\boldsymbol {v} ^ {2}}{2} \right\rangle . \tag {2.10} \\ \end{array}
$$

From Eqs. (2.5) and (2.6) we immediately realize that such equations are left invariant under the simultaneous transformation ${ \pmb { g } }  - { \pmb { g } } , \phi  - \phi$ . As a consequence, $\langle \phi v ^ { 2 } / 2 \rangle = 0$ and the resulting kinetic energy simply reads: 

$$
E _ {K} = \rho_ {o} \left\langle \frac {\boldsymbol {v} ^ {2}}{2} \right\rangle . \tag {2.11}
$$

By defining $E _ { \mathcal { F } } \equiv \mathcal { F } / \Omega$ the total energy of the two-fluid system is 

$$
E = E _ {P} + E _ {K} + E _ {\mathcal {F}}.
$$

The equation for $E _ { K }$ is obtained by multiplying Eq. (2.6) by $\rho _ { o } v _ { \alpha }$ and then taking spatial average. We easily get: 

$$
d E _ {K} / d t = \rho_ {o} \partial_ {t} \langle \frac {\boldsymbol {v} ^ {2}}{2} \rangle = - \rho_ {o} \nu \langle (\partial_ {\alpha} \boldsymbol {v}) ^ {2} \rangle + \rho_ {o} \mathcal {A} g \langle v \phi \rangle - \Lambda \langle v _ {\alpha} (\partial_ {\alpha} \phi) (\partial^ {2} \phi) \rangle . \tag {2.12}
$$

Let us now take $\operatorname { E q . }$ . (2.5), multiply it by y, and take the average: 

$$
\partial_ {t} \langle y \phi \rangle + \langle y \partial_ {y} (v \phi) \rangle = \gamma \Lambda \langle y \partial^ {2} \left(- \partial^ {2} \phi + \frac {\phi^ {3} - \phi}{\epsilon^ {2}}\right) \rangle = 0, \tag {2.13}
$$

by translational invariance and Leibniz rule. We thus have: 

$$
d E _ {P} / d t = - \partial_ {t} \left(\rho_ {o} \mathcal {A} g \langle y \phi \rangle\right) = - \rho_ {o} \mathcal {A} g \langle v \phi \rangle , \tag {2.14}
$$

where we have used the fact that $\langle y \partial _ { y } ( v \phi ) \rangle = - \langle ( \partial _ { y } y ) v \phi \rangle = - \langle v \phi \rangle$ . The free–energy variation is 

$$
\begin{array}{l} \partial_ {t} \mathcal {F} = \int \int \frac {\delta \mathcal {F}}{\delta \phi} \frac {\partial \phi}{\partial t} d x d y = \\ = \int \int \frac {\delta \mathcal {F}}{\delta \phi} \left[ - \boldsymbol {v} \cdot \boldsymbol {\partial} \phi + \gamma \partial^ {2} \left(\frac {\delta \mathcal {F}}{\delta \phi}\right) \right] d x d y = \\ = - \gamma \left\langle \left[ \partial \left(\frac {\delta \mathcal {F}}{\delta \phi}\right) \right] ^ {2} \right\rangle \Omega - \int \int \left(\frac {\delta \mathcal {F}}{\delta \phi}\right) v _ {i} \partial_ {i} \phi d x d y = \\ = - \gamma \langle \left[ \partial \left(\frac {\delta \mathcal {F}}{\delta \phi}\right) \right] ^ {2} \rangle \Omega - \Lambda \int \int \left[ (- \partial^ {2} \phi) + \frac {\phi^ {3} - \phi}{\epsilon^ {2}} \right] v _ {i} \partial_ {i} \phi d x d y = \\ = - \gamma \left\langle \left[ \partial \left(\frac {\delta \mathcal {F}}{\delta \phi}\right) \right] ^ {2} \right\rangle \Omega - \Lambda \left\langle e _ {\alpha \beta} \left(\partial_ {\alpha} \phi\right) \left(\partial_ {\beta} \phi\right) \right\rangle \Omega . \tag {2.15} \\ \end{array}
$$

i.e., 

$$
\partial_ {t} E _ {\mathcal {F}} = - \gamma \left\langle \left[ \partial_ {\alpha} \left(\frac {\delta \mathcal {F}}{\delta \phi}\right) \right] ^ {2} \right\rangle - \Lambda \left\langle e _ {\alpha \beta} (\partial_ {\alpha} \phi) (\partial_ {\beta} \phi) \right\rangle . \tag {2.16}
$$

where we have introduced the strain tensor $e _ { \alpha \beta } \equiv \left( \partial _ { \alpha } v _ { \beta } + \partial _ { \beta } v _ { \alpha } \right) / 2$ and assumed boundary conditions suitable to justify integrations by parts. 

The energy balance takes then the form: 

$$
\partial_ {t} (E _ {K} + E _ {P} + E _ {\mathcal {F}}) = - \rho_ {o} \nu \langle (\partial_ {\alpha} \pmb {v}) ^ {2} \rangle - \gamma \langle \left[ \partial_ {\alpha} \left(\frac {\delta \mathcal {F}}{\delta \phi}\right) \right] ^ {2} \rangle . (2. 1 7)
$$

The global system in thus intimately dissipative, even for a vanishing kinetic viscosity. It is worth emphasizing the cancellation of $\Lambda \langle e _ { \alpha \beta } ( \partial _ { \alpha } \phi ) ( \partial _ { \beta } \phi ) \rangle$ by the kinetic and the free–energy contributions, due to exchanges between the velocity field and the interface. 

# 2.4. Dispersion relation for the phase-field model

The aim of this section is to show that the well-known dispersion relation for gravitycapillary waves (Chandrasekhar 1961) can be easily obtained within the phase-field formalism. To do that, let us concentrate our attention on a two-dimensional problem and indicate by y the gravity direction. Moreover, we will assume heavier fluid to be placed below the lighter one, in a way to have a stable situation. For a given perturbation imposed to the interface, the problem is to determine how the perturbation evolves in time. 

Denoting by $h ( x , t )$ a small perturbation imposed to a planar interface, we can rewrite $\phi$ as: 

$$
\phi = f \left(\frac {y - h (x , t)}{\epsilon}\right) \quad , \tag {2.18}
$$

where h can be larger than ǫ, yet it has to be smaller than the scale of variation of h (small amplitudes). 

Locally, the interface is in equilibrium, i.e.: 

$$
f ^ {\prime \prime} = V ^ {\prime} (f), \tag {2.19}
$$

where $V ( \phi ) = ( \phi ^ { 2 } - 1 ) ^ { 2 } / 4 \epsilon ^ { 2 }$ . In this limit we have: 

$$
\mu = - \Lambda \frac {\partial^ {2} f}{\partial x ^ {2}} = \frac {\Lambda}{\epsilon} \left[ f ^ {\prime} \frac {\partial^ {2} h}{\partial x ^ {2}} - \frac {f ^ {\prime \prime}}{\epsilon} \left(\frac {\partial h}{\partial x}\right) ^ {2} \right] . (2. 2 0)
$$

Linearizing Eq. (2.6) for small interface velocity we have, neglecting the viscous term: 

$$
\rho_ {o} \partial_ {t} v = - \partial_ {y} p - \phi \partial_ {y} \mu - \mathcal {A} g \rho_ {o} \phi . \tag {2.21}
$$

The integration in the vertical direction interpreted in the principle value sense 

$$
q _ {y} := \lim _ {L \uparrow \infty} \int_ {- L} ^ {L} v d y, \tag {2.22}
$$

$$
\rho_ {o} \partial_ {t} q _ {y} := \lim _ {L \uparrow \infty} \left\{\frac {\Lambda}{\epsilon} \int_ {- L} ^ {L} \left[ f f ^ {\prime \prime} \frac {\partial^ {2} h}{\partial x ^ {2}} - \frac {1}{\epsilon} f f ^ {\prime \prime \prime} \left(\frac {\partial h}{\partial x}\right) ^ {2} \right] d (y / \epsilon) - \mathcal {A} g \rho_ {o} \int_ {- L} ^ {L} f d y \right\} (2. 2 3)
$$

yields: 

$$
\rho_ {o} \partial_ {t} q _ {y} = \sigma \frac {\partial^ {2} h}{\partial x ^ {2}} - 2 \mathcal {A} g \rho_ {o} h, (2. 2 4)
$$

having used the relations $\textstyle \int ( f ^ { \prime } ) ^ { 2 } d y = 2 { \sqrt { 2 } } / 3 , \int f f ^ { \prime \prime \prime } d y = 0$ and 

$$
\lim _ {L \uparrow \infty} \int_ {- L} ^ {+ L} f d y = + 2 h. \tag {2.25}
$$

The height variation of the interface has to match the vertical fluid velocity, thus giving: 

$$
\partial_ {t} h = v (x, h (x, t), t) \equiv v ^ {(i n t)} (x, t). \tag {2.26}
$$

The last step is to relate the velocity at the interface with the integral $q _ { y } .$ . This is done by restricting to potential flows: 

$$
\boldsymbol {v} = \partial \psi \quad \partial^ {2} \psi = 0. \tag {2.27}
$$

For $y > 0$ , denoting with “ˆ” the Fourier Transform, we have: 

$$
\psi (x, y, t) = \int_ {0} ^ {\infty} e ^ {- k y + i k x} \hat {\psi} (k, t) d k + \text { c.c. } \tag {2.28}
$$

$$
v (x, y, t) = - \int_ {0} ^ {\infty} k e ^ {- k y + i k x} \hat {\psi} (k, t) d k + \text { c.c. } \tag {2.29}
$$

$$
q _ {y} (x, t) = - 2 \int_ {0} ^ {\infty} e ^ {i k x} \hat {\psi} (k, t) d k + \text { c.c. } \tag {2.30}
$$

$$
v ^ {(i n t)} = - \int_ {0} ^ {\infty} k e ^ {i k x} \hat {\psi} (k, t) d k + \text { c.c. } \tag {2.31}
$$

Therefore: 

$$
\hat {v} ^ {(i n t)} = \frac {k \hat {q} _ {y}}{2}, \tag {2.32}
$$

so that in k−space we have: 

$$
\partial_ {t} \hat {h} = \frac {k \hat {q} _ {y}}{2} \rho_ {o} \partial_ {t} \hat {q} _ {y} = (- \sigma k ^ {2} - 2 \mathcal {A} g \rho_ {o}) \hat {h}. (2. 3 3)
$$

From these two equations we immediately get: 

$$
\partial_ {t} ^ {2} \hat {h} + \omega^ {2} \hat {h} = 0, \tag {2.34}
$$

with: 

$$
\omega^ {2} (k) = + \mathcal {A} g k + \frac {\sigma}{2 \rho_ {o}} k ^ {3} \tag {2.35}
$$

that is the expected dispersion relation (Chandrasekhar 1961). For the stable configuration we have, for all values of σ: $\mathcal { A } g k + \sigma / \left( 2 \rho _ { o } \right) k ^ { 3 } > 0$ , i.e. any initially imposed perturbation will not grow indefinitely. 

From Eq. (2.34) and the initial condition: 

$$
\partial_ {t} \hat {h} (k, t) = 0 \quad \text { at } t = 0 \quad , \tag {2.36}
$$

we immediately have: 

$$
\hat {h} (k, t) = \hat {h} (k, 0) \cos (\omega t) \tag {2.37}
$$

and the velocity at the interface reads: 

$$
\hat {v} _ {y} ^ {i n t} (k, t) = - \hat {h} (k, 0) \omega \sin (\omega t). \tag {2.38}
$$

Assuming an initial perturbation of the form $h ( x , 0 ) = h _ { 0 } \cos \left( \bar { k } x \right)$ , from Eqs. (2.31) and (2.38) we obtain: 

$$
\hat {\psi} (\bar {k}, t) = \frac {1}{\bar {k}} \hat {h} (\bar {k}, 0) \omega \sin (\omega t), \tag {2.39}
$$

and the velocity components, for $y > 0$ , read: 

$$
v ^ {\uparrow} (x, y, t) \equiv v (x, y, t) = - \cos (\bar {k} x) e ^ {- \bar {k} y} h _ {o} \omega \sin (\omega t) \tag {2.40}
$$

$$
u ^ {\uparrow} (x, y, t) \equiv u (x, y, t) = \sin (\bar {k} x) e ^ {- k y} h _ {o} \omega \sin (\omega t), \tag {2.41}
$$

where we used the relation $h _ { 0 } = 2 \hat { h } ( k , 0 )$ . 

For $y < 0$ , in a similar way we obtain the velocity field components: 

$$
v ^ {\downarrow} (x, y, t) \equiv v (x, y, t) = - \cos (\bar {k} x) e ^ {+ \bar {k} y} h _ {o} \omega \sin (\omega t) \tag {2.42}
$$

$$
u ^ {\downarrow} (x, y, t) \equiv u (x, y, t) = - \sin {(\bar {k} x)} e ^ {+ \bar {k} y} h _ {o} \omega \sin {(\omega t)} \quad . \tag {2.43}
$$

When in the initial configuration the heavier fluid placed above the lighter one, the dispersion relation (2.35) trasforms in: 

$$
\omega^ {2} (\bar {k}) = - \mathcal {A} g \bar {k} + \frac {\sigma}{2 \rho_ {o}} \bar {k} ^ {3}, \tag {2.44}
$$

which is readly obtained by flipping the sign of $g .$ For $\sigma < \sigma _ { c } \equiv 2 \rho _ { o } / ( A g \bar { k } ^ { 2 } )$ surface tension is not able to contrast gravity-induced vertical motion with the final result that amplitude perturbations grows exponentially: the flow is unstable. More precisely, from relation (2.44) and for $\sigma < \sigma _ { c }$ we have: 

$$
\omega (\bar {k}) = \sqrt {- \mathcal {A} g \bar {k} + \frac {\sigma}{2 \rho_ {o}} \bar {k} ^ {3}} \equiv i \alpha (\bar {k}), \tag {2.45}
$$

and Eqs. (2.40) - (2.43) transform in: 

$$
v ^ {\uparrow} (x, y, t) \equiv v (x, y, t) = \cos (\bar {k} x) e ^ {- \bar {k} y} h _ {0} \alpha \sinh (\alpha t) \tag {2.46}
$$

$$
u ^ {\uparrow} (x, y, t) \equiv u (x, y, t) = - \sin {(\bar {k} x)} e ^ {- \bar {k} y} h _ {0} \alpha \sinh {(\alpha t)} \quad , \tag {2.47}
$$

for $y > 0$ , and: 

$$
v ^ {\downarrow} (x, y, t) \equiv v (x, y, t) = \cos (\bar {k} x) e ^ {+ k y} h _ {0} \alpha \sinh (\alpha t) \tag {2.48}
$$

$$
u ^ {\downarrow} (x, y, t) \equiv u (x, y, t) = \sin (\bar {k} x) e ^ {+ \bar {k} y} h _ {0} \alpha \sinh (\alpha t), \tag {2.49}
$$

for $y < 0 .$ 

# 3. Numerical investigation

In this section we report results we have obtained exploiting direct numerical simulations (DNS) of the phase-field model for the Rayleigh–Taylor problem described in the preceeding sections. Our attention will be focused both on the linear phase of the perturbation evolution and on the weakly nonlinear regime governed by plumes, for A ≪ 1. In the present study we will consider initial perturbations imposed to the interface varying along one of the horizontal directions, say the x-axis, and invariant along the other horizontal direction, say the z-axis. The perturbation is thus intimately two-dimensional a fact that allows us to solve the original Navier–Stokes equations coupled to the phase field in two dimensions. This clearly permits to obtain high accuracy and thus to properly test the phase-field approach against known results for both the linear and the nonlinear evolution stage. 

For a two-dimensional flow it is convenient to introduce the vorticity field ω $[ \omega = ( \partial \times \pmb { v } ) _ { z } ]$ ] and study the equations 

$$
\partial_ {t} \omega + \boldsymbol {v} \cdot \boldsymbol {\partial} \omega = + \nu \partial^ {2} \omega - \frac {\Lambda}{\rho_ {o}} \boldsymbol {\partial} \times (\partial^ {2} \phi \boldsymbol {\partial} \phi) - \mathcal {A} (\boldsymbol {\partial} \phi) \times \boldsymbol {g} \tag {3.1}
$$

$$
\partial_ {t} \phi + \boldsymbol {v} \cdot \boldsymbol {\partial} \phi = \gamma \boldsymbol {\partial} ^ {2} \mu = \gamma \Lambda \boldsymbol {\partial} ^ {2} \left[ - \partial^ {2} \phi + \frac {(\phi^ {3} - \phi)}{\epsilon^ {2}} \right] \quad . \tag {3.2}
$$

In order to efficiently and accurately solve those equations we exploit a pseudospectral method . Accordingly, periodic boundary conditions have to be assumed along the two directions. For the horizontal direction it is a natural choice (see e.g. Cabot & Cook (2006); Liu & Shen (2003)) while along the vertical one this choice deserves some comments. As initial condition we started from the hyperbolic-tangent profile, Eq. (2.3), for φ with the interface placed in the middle of the domain. The fact that we have periodic boundary conditions along y simply means that far from the middle of the domain the hyperbolic-tangent profile has to be distorted in order to satisfy periodic boundary conditions. However, both in the linear and in the weakly nonlinear regimes the amplitude of the interface perturbation is always much smaller than the vertical size of the box, so that the actual choice of boundary conditions at the top and bottom can be safely neglected. 

Such a strategy has been already exploited for the miscible case by Celani, Mazzino & Vozella (2006). 

The box has a horizontal to vertical aspect ratio $L _ { x } / L _ { y } = 1$ for the linear analysis stage and $L _ { x } / L _ { y } = 1 / 2$ for the weakly nonlinear evolution. In the latter case we take a smaller aspect ratio owing to the fact that the perturbation can reach a higher amplitude (with respect to case of the linear analysis). 

In both cases the resolution is 1024 × 1024 collocation points. We need such a high resolution (despite the fact that we focus on a linear and weakly nonlinear study) in order to have a well described interface separating the two phases. In our simulations the mixing width $( \sim 4 \epsilon )$ is 6 mesh points. 

The time evolution is implemented by a standard second-order Runge–Kutta scheme. 

The physically relevant parameters in the present problem are the kinematic viscosity $\nu ,$ the buoyancy intensity $\mathcal { A } g$ and the surface tension σ. Both $\mathcal { A } g$ and ν will be varied in our study, while σ will be kept fixed to a fixed value (see below). The surface tension is related to the ratio $\Lambda / \epsilon$ with ǫ (and thus $\Lambda )$ sufficiently small in order to have a finite value for the surface tension and, at the same time, to reproduce the correct sharp-interface limit. Finally, the parameter γ appearing in the relaxation term in $\operatorname { E q . }$ (3.2) must satisfy the requirement that $\gamma \Lambda$ be small, so as to enforce ‘istantaneous’ local equilibrium between flow and interface. Here we used the value (model units) $\gamma \Lambda = 1 0 ^ { - 8 }$ . 

All simulations presented here start from an initial condition corresponding to an equilibrium configuration: velocity identically zero and hyperbolic tangent profile for the phase field φ, expressed by the relation of the form: tanh $( ( y - h ( x , t = 0 ) ) / c )$ with 

$$
h (x, \mathrm{t} = 0) = h _ {0} \sin (k x).
$$

For a given k we choose the initial amplitude $h _ { 0 }$ in a way that $h _ { 0 } / \lambda ( \mathrm { w h e r e } \lambda \equiv 2 \pi / k )$ 号 is sufficiently small to fall in the linear phase $\left( \mathrm { i . e . ~ } h _ { 0 } / \lambda \ll 1 \right)$ and $h _ { 0 }$ is sufficiently large for the wave disturbance to see an almost infinitesimal mixing width $\left( \mathrm { i . e . ~ } h _ { 0 } / \epsilon \gg 1 \right)$ . Specific numerical values are reported in the next sections. 

# 3.1. Linear instability for negligible viscosity

The aim of this section is to verify the growth-rate (2.45) which holds in the linear phase when the viscosity is negligible. 

In order to do so, we take a small value of $\nu \ ( \nu = 1 0 ^ { - 5 }$ in the model units) and vary k (up to $k _ { c } \equiv ( 2 \mathcal { A } g \rho _ { o } / \sigma ) ^ { 1 / 2 }$ , the critical wave-number separating unstable from stable wave-modes) and $A g$ and take a fixed value of $\sigma .$ . The ratio $h _ { 0 } / \lambda = 0 . 0 6$ while $h _ { 0 } / \epsilon$ ranges from  10 to  40 in the range of k considered. 

The behavior of the square growth-rate $\alpha ^ { 2 }$ is shown in dimensionless form in Fig. 2 as a function of k for three different values of $k _ { c }$ (obtained by varying $A g )$ and in Fig. 3 by varying $A g$ for three different values of $k < k _ { c }$ . In both figures, symbols refer to the numerical results and the dashed line is the theoretical expectation given by (2.45). 

![](images/74dfa235f68faf3f8d0f3db2a74f56195e594c2b8621b732c946a4cc31f9c10c.jpg)



Figure 2. The square growth-rate $\alpha ^ { 2 }$ (see $\operatorname { E q } .$ (2.45)) for three different values of Ag corresponding to three different values of the critical wave number $k _ { c } \equiv ( 2 . A g \rho _ { o } / \sigma ) ^ { 1 / 2 } \colon k _ { c } = 3 . 4$ (solid circle), $k _ { c } = 4 . 7$ (solid triangle) and $k _ { c } = 5 . 7$ (solid rhombus). The dashed line is the linear-theory prediction expressed by the relation (2.45).


![](images/40074b4688e86dd0e1ac433fc317450696f45d707a0d9e7b42fdab87d2045497.jpg)



Figure 3. The square growth-rate $\alpha ^ { 2 }$ for $k = 1$ (solid circles), $k = 2$ (solid triangles) and $k = 3$ (solid rhombus), all smaller than $k _ { c } ,$ for six different values of Ag ranging from 0.11 to 0.61. The dashed line corresponds to the linear-theory prediction.


The numerical data in Figs. 2 and 3 have been obtained via best-fit of $\langle v ^ { 2 } \rangle$ , the spatial average of $v ^ { 2 }$ as a function of time. The latter average is computed over a horizontal strip containing the interface (placed in the middle of the computational domain) and having an extension of $a _ { y }$ above and below the interface. This has been done to avoid spurious contaminations coming from the upper and lower domain regions affected by the boundary conditions. In formulae: 

$$
\begin{array}{l} \langle v ^ {2} \rangle = \frac {1}{2 a _ {y}} \frac {1}{L _ {x}} \int_ {- a _ {y}} ^ {0} d y \int_ {0} ^ {L _ {x}} d x (v ^ {\downarrow}) ^ {2} + \frac {1}{2 a _ {y}} \frac {1}{L _ {x}} \int_ {0} ^ {a _ {y}} d y \int_ {0} ^ {L _ {x}} d x (v ^ {\uparrow}) ^ {2} \\ = \frac {1}{2 a _ {y} k} \left[ - e ^ {- 2 k a _ {y}} + 1 \right] \alpha^ {2} h _ {0} ^ {2} \sinh^ {2} (\alpha t), \tag {3.3} \\ \end{array}
$$

where we used the expression (2.46) and (2.48) for $v ^ { \uparrow }$ and $v ^ { \downarrow }$ , respectively. 

The best fit has been done with α as unique free parameter and its high accuracy can be verified in Fig. 4 where we show the time evolution of $\langle v ^ { 2 } \rangle$ for $k _ { c } = 4 . 7$ (solid triangles in Fig. 2) and for four values of k smaller than $k _ { c } . ~ \mathrm { A t }$ tα > 1.5 nonlinear effects start to enter into play giving rise to corrections to the linear analysis (see Sec. 3.4). Up to that time, linear theory is very accurate as one can also realize by looking at the insets of Fig. 2 where the sinusoidal form of $h ( x , t )$ is reported for $t \alpha = 1 . 5$ . 

![](images/62af291469cbec557c1d8c2e3e8366c5f6f057e90a7aa6ffcd4657a888647740.jpg)



Figure 4. Time behavior of $\langle v ^ { 2 } \rangle$ for $k _ { c } = 4 . 7$ (in Fig. 2 corresponding to the solid triangle) and for four values of $k < k _ { c } . ~ ( \mathrm { a } ) \stackrel { . } { k } = 1$ , (b) k = 2, (c) k = 3 and (d) $k = 4 .$ The numerical results (symbols) are compared with the corresponding best fit expressions (see the text for details). In the insets the interface perturbation, $h ( x , t )$ , is plotted at $t \alpha = 1 . 5$ revealing a very accurate linear analysis prediction.


# 3.2. Linear instability for finite viscosity

The aim of this section is to investigate numerically how the growth-rate, α, is modified by viscosity. As discussed in Appendix A, both an upper and a lower bound for the perturbation growth-rate are known (see Eqs. (A 1) and (A 2)) and we want to assess how the actual growth-rates compare with those. 

For such purpose, we choose a surface tension, σ, and $\mathcal { A } g$ in such a way to obtain instability for few (unstable) wavenumbers. Our choice was $k _ { c } = 5 . 7$ (see Sec. 3.1) thus corresponding to 5 unstable wavenumbers. 

As far as the initial perturbation is concerned, we report here the case corresponding to $k = 1$ . Initial perturbations with a larger wavenumber simply need an initial smaller amplitude (and eventually a larger numerical resolution) in order to satisfy $h _ { 0 } \gg \epsilon$ and $h _ { 0 } \ll \lambda$ . Here, we have $h _ { 0 } / \lambda = 0 . 0 3$ and $h _ { 0 } / \epsilon \sim 2 0$ . Such ratios turned out to be sufficiently ‘asymptotic’ to produce accurate results. The effect of viscosity is studied by considering twelve values of viscosity in the range $1 0 ^ { - 5 } \le \nu \le 5 ~ 1 0 ^ { - 2 }$ (model units). 

The results of our simulations are summarized in Fig. 5 where the behavior of the square perturbation growth-rate, $\alpha _ { \nu } ^ { 2 } .$ , is shown as a function of viscosity. The numerical predictions have been compared with the available theoretical bounds (dashed lines). 

Note that the numerical points are always in between the two bounds and also how the relative differences between the upper bound and the numerical values $\mathrm { a r e } < 1 1 \%$ . This latter fact is compatible, for example, with the results of Menikoff et al. (1977). 

![](images/00061d705552c508aa4b2d49c8a0aac76d5ede5c61b837b4b67f9385455967fb.jpg)



Figure 5. Behavior of the dimensionless perturbation growth-rate, $\alpha _ { \nu } .$ , for $k = 1$ and $\mathcal { A } g$ corresponding to $k _ { c } = 5 . 7$ . Dotted lines correspond to upper and lower bounds for the growth-rate (see Eqs. (A 2) and $\left( \mathrm { A } 1 \right) )$ . The arrow selects a value of the viscosity for which the time evolution of $\langle v ^ { 2 } \rangle$ is reported in the inset. The continuous line is the best fit slope (see text).


The value of the growth-rates have been obtained via best of $\langle v ^ { 2 } \rangle$ (see Eq. (3.3)). Unlike what we did in previous section, here we perform the fit within the exponential region. The reason is that the non-asymptotic form of the perturbation time-evolution is unknown in the present case. 

The fit accuracy can be appreciated in the inset of Fig. 5 where the temporal evolution of the pertubation for $\nu = 0 . 3$ (model units) is shown together with the best fit slope (dashed line) from which $\alpha _ { v }$ is determined. Error bars, estimated by looking at the fit sensitivity by varying the length of the fit interval, are of the order of the symbol sizes. 

# 3.3. Stable configuration: gravity-capillary waves

The performance of the phase-field approach in the unstable regime predicted by linear theory both in the presence and in the absence of viscosity proved to be very good. As discussed in Sec. 2.4, for sufficiently large surface tensions and/or sufficiently small differences between fluids density, a perturbation initially imposed to the fluid interface may maintain its initial amplitude giving rise to the dispersion relation (2.35). The waves resulting from the balance between gravity and surface tension are known as gravitycapillary waves. Our aim here is to verify their dispersion relation. 

To do that, we have fixed the parameters to obtain a critical wavenumber of order one. For $\mathcal { A } g = 0 . 0 0 8$ (model units) and the same σ as in the unstable case, one has $k _ { c } = 0 . 9$ . The first accessible wavenumber is thus stable and should evolve in time according to (2.35). However, the geometrical/computational configuration used in the unstable case did not produce sufficiently accurate results. In particular, using the same domain aspect ratio $L _ { x } / L _ { y } = 1$ and the same ratio between perturbation amplitude and perturbation wave-length we found a dynamics too dissipative with respect to what is expected. In the absence of viscosity, dissipation arises in the phase field formulation due to the sole contribution proportional to γ in Eq. (2.17). The latter parameter has been taken sufficiently small to ensure a negligible effects inside a period of oscillation. The specific value was $\gamma = 6 . 2 5 \times 1 0 ^ { - 5 }$ . To avoid spurious dissipation, as that induced by nonlinear effects, we reduced the amplitude of the initial perturbation with respect to the unstable case. Also, we increased the size of the periodicity box along the gravitational direction in a way to reduce possible spurious contribution arising from the upper/lower part of the computational domain where instabilities, not present in the unstable case, might now develop. The above choice on the amplitude of the inital perturbation implies a consequent reduction of $\epsilon .$ . The following set of parameters have been used: $\epsilon = 0 . 0 0 8$ , $L _ { x } / L _ { y } = 1 / 4$ and a resolution $N x \times N y$ of $2 5 6 \times 4 0 9 6$ . For an initial perturbation on $k = 1$ , its initial amplitude $h _ { 0 }$ has been chosen to have $h _ { 0 } / \lambda = 0 . 0 1 2$ and $h _ { 0 } / \epsilon \sim 1 0$ . The behavior of the maximum, $\eta ( t )$ , of the initial perturbation is shown as a function of time in Fig. 6. The continuous line is relative to a sinusoidal with pulsation ω obtained from (2.35). The agreement between theory and numerics is satisfactory both for the amplitude and for the pulsation. Note the small reduction of $\eta ( t )$ , in one oscillation period: only 1 grid box over 4096. 

![](images/51abfe128ba973e800945915b05a8d4143619af913a154b41171764bc0adea7d.jpg)



Figure 6. Time behavior of the perturbation maximum, $\eta ( t )$ , for $k = 1$ and $h _ { 0 } / \lambda = 0 . 0 1 2$ . The critical wave number is $k _ { c } = 0 . 9$ . Numerical results (symbols) are compared with the prediction from linear theory (see Eq. (2.37)).


# 3.4. Weakly non-linear stage

In this section we investigate the early stages of the nonlinear dynamics. We focus on the rising/falling velocity of plumes in the limit of small Atwood numbers when spikes and bubbles are known to coincide. The theoretical prediction for the terminal velocity is reported in Appendix B. Our aim here is both to verify the existence of a regime characterized by a costant ‘terminal’ velocity and, secondly, to compare the prediction (B 1) for such terminal velocity with our numerical data. 

The physical parameters are chosen to magnify the effect of the surface tension on the terminal velocity. This happens when the wavenumber k of the initial pertubation (still supposed unimodal) is slightly below $k _ { c }$ . Here we choose $\mathcal { A } g$ and σ such that $k _ { c } ~ =$ 4.004 and thus look at the dynamics associated to the wavenumber $k = 4$ . The initial perturbation has an amplitude $h _ { 0 } / \lambda = 0 . 0 6 ;$ the initial dynamics is thus linear. Although we are interested to investigate the case of zero viscosity, in order to prevent numerical instabilities we add a small viscosity $\nu = 2 \times 1 0 ^ { - 5 }$ (model units). In $\mathrm { F i g . 7 }$ the perturbation amplitude is shown as a function of time: symbols correspond to our numerical data and the dashed line is the prediction (B 1). A good agreement is found between numerics and theory in the range $1 . 2 < t U / \lambda < 1 . 8$ . At larger times, neighboring plumes start to interact and the arguments leading to (B 1) do not apply any longer. In Fig. 8 we show some snapshots of the evolution of the two fluids. Figures are equally spaced in time in the interval $1 . 2 < t U / \lambda < 1 . 8$ . Black corresponds to $\phi = - 1 ;$ ; white to $\phi = 1$ . Their shape is similar to that experimentally observed. Note the aforementioned spike/bubble symmetry corresponding to the up-down symmetry of our original evolution equations. by Waddell, Nieserhaus & Jacobs (2001). 

![](images/ced2218ff056415f5b013174c76e62aa03347d0209b5e4ea99155abb757610fa.jpg)



Figure 7. Time evolution of amplitude perturbation η(t). The dots are our numerical results, the dashed line is the prediction by Eq. (B 1).


![](images/a802e84715644113e69ef1f7a8a9133c22fe567d482c9fa8b02a88925857c5c2.jpg)



Figure 8. Two-color snapshots of the phase field. Black (white) corresponds to $\phi = - 1$ $( \phi = 1 )$ . Frames are equally spaced in time in the interval $1 . \dot { 2 } < t \dot { U } / \lambda < \bar { 1 . 8 }$ (see also Fig. 7).


# 4. Conclusions and perspectives

In this paper we showed that the phase–field model provides a valuable numerical instrument for the study of immiscible, convective hydrodynamics. As a testground for this model, we have considered the Rayleigh–Taylor instability. Numerical results compare very well with known analytical results both for the linearly stable and unstable case, and for the weakly nonlinear stages of the latter. 

All these results are very encouraging in view of the next important step that is the the numerical simulation of immiscible RT turbulence. There, the interplay of all the fundamental mechanisms that we have illustrated here (instabilities and wave propagation) is expected to give rise to a small-scale emulsion-like phase dominated by gravity-capillary waves and by a large-scale hydrodynamic range of scales where classical Kolmogorov turbulence should appear. This theoretical suggestion still awaits numerical confirmation, and the phase–field model provides the appropriate method to pursue this goal. 

We acknowledge useful discussions with Hekki Haario. AM and LV have been partially supported by PRIN 2005 project n. 2005027808 and by CINFAI consortium (AM). LV acknowledges support from From Discrete to Continuous models for Multiphase Flows TEKES project n. 40289/05. 

# Appendix A. Bounds for the perturbation growth-rate in the presence of viscosity

The effect of viscosity is to reduce the perturbation growth-rate. However it does not remove the instabilities. Analytically, it is more difficult to consider the effect of viscosity with respect to surface tension (see Eq. (115) at page 443 of Chandrasekhar 1961). Nonetheless, it is possible to determine a lower and an upper bound to the growth-rate $\alpha _ { \nu }$ . These bounds are the solutions to the following equations (Menikoff et al. 1977): 

$$
\alpha_ {\nu} ^ {4} + 2 \nu k ^ {2} \alpha_ {\nu} ^ {3} + (\nu^ {2} k ^ {3} - \frac {\alpha^ {2}}{k}) k \alpha_ {\nu} ^ {2} - (\nu^ {2} k ^ {3} + \frac {\alpha^ {2}}{k}) \nu k ^ {3} \alpha_ {\nu} - (\nu^ {4} k ^ {6} - \frac {\alpha^ {4}}{k ^ {2}}) k ^ {3} = 0 \quad (A 1)
$$

$$
\alpha_ {\nu} ^ {2} + 2 \nu k ^ {2} \alpha_ {\nu} - \alpha^ {2} = 0. (A 2)
$$

where α is the growth-rate in the inviscid case (see Eq. (2.45)). The solution of Eq. (A 2) is: 

$$
\alpha_ {\nu} = - k ^ {2} \nu + \sqrt {k ^ {4} \nu^ {2} + \alpha^ {2}} \tag {A3}
$$

while only a numerical solution is available for $\mathrm { E q . \ ( A 1 ) }$ . 

The goodness of those upper and lower bounds are numerically investigated in Sec. 3.2 by means of the phase-field model. 

# Appendix B. Models for the terminal bubbles/spike velocities in the weekly nonlinear regime

Substantial deviations from the linear theory are observed when the perturbation amplitude reaches a size of the order of $0 . 1 \lambda - 0 . 4 \lambda$ (Sharp 1984). 

In that case the perturbation evolution is nonlinear. Then the disturbance grows nonlinearly and the interface starts to deform. Indeed, at least for finite values of A, the interface can be divided into spikes corresponding to the regions where the heavier fluid penetrates into the lighter one, and bubbles associated to those regions where lighter fluid rises in the heavier one. The roll-up of vortices produces a mushroom-type shape for bubbles and spikes (see, for instance, Waddell, Nieserhaus & Jacobs 2001). When the fluid densities are similar (corresponding to our case $\mathcal { A } \ll 1 )$ spikes and bubbles coincide and approach a constant and equal velocity. In both cases, the exponential growth of the velocity perturbation amplitude characterizing the linear phase of the evolution is replaced by a linear-in-time behavior (Waddell, Nieserhaus & Jacobs 2001). Two models are available to describe this stage: the drag-buoyancy model (Alon et al. 1995) and the “Layzer model” (Layzer 1955; Goncharov 2003; Young & Ham 2006). The former model describes bubble and spike motion by balancing the buoyancy and drag forces and it assumes that this velocities reach a constant values for sufficiently long times. The latter model uses an expansion of the perturbation amplitudes and conservation equations near the tip of bubbles and spikes. This approach has been first applied to the fluidvacuum interface (A = 1) (Layzer 1955) and then extended to arbitrary Atwood numbers (Goncharov 2003) and to include the surface tension contribution (Young & Ham 2006). According to the latter study, in our case (bidimensional flow, immiscible fluids and small Atwood number) one expects that the terminal bubble and spike velocity be equal to (Young & Ham 2006): 

$$
U (t \rightarrow \infty) = \sqrt {\frac {2}{3} \mathcal {A} \frac {g}{k} - \frac {1}{9} \frac {\sigma}{\rho_ {2} + \rho_ {1}} k}. \tag {B1}
$$

This expectation is numerically tested, in Sec. 3.4, by exploiting the phase-field method. 

# REFERENCES



Alon, U., Hecht, J., Ofer, D. & Shvarts, D. 1995 Power laws and similarity of Rayleigh– Taylor and Richtmyer–Meshkov mixing fronts at all density ratio. Phys. Rev. Lett. 74(4), 534–537 





Anderson, D. M., McFadden, G. B. & Wheeler, A. A. 1998 Diffuse-interface methods in fluid mechanics. Annu. Rev. Fluid Mech. 30, 139–165 





Badalassi, V. E., Ceniceros & H. D., Banerjee, S. 2003 Computation of multiphase systems with phase field models. J. Comput. Phys. 190, 371–397 





Berti, S., Boffetta, G., Cencini, M. & Vulpiani, A. 2005 Turbulence and coarsening in active and passive binary mixtures. Phys. Rev. Lett. 95, 224501-1–224501-4 





Bray, A. J. 2002 Theory of phase-ordering kinetics. Advances in Physics 51(2), 481–587 





Cabot, W. H. & Cook, A. W. 2006 Reynolds number effects on Rayleigh–Taylor instability with possible implications for type-Ia supernovae. nature physics 2, 562–568 





Cahn, J. W. & Hilliard, J. E. 1958 Free energy of a non uniform system. I. Interfacial free energy. J. Chem. Phys. 28, 258–267 





Canuto, C., Hussaini, M. Y., Quarteroni, A. & Zang, T. A. 1988 Spectral Methods in Fluid Dynamics Springer Series in Computational Physics. Springer-Verlag 





Celani, A.,Mazzino, A. & Vozella, L. 2006 Rayleigh–Taylor turbulence in two-dimensions. Phys. Rev. Lett. 96, 134504-1–134504-4 





Chandrasekhar, S. 1961 Hydrodynamic and Hydromagnetic Stability. New York: Dover 





Chertkov, M. 2003 Phenomenology of Rayleigh–Taylor turbulence. Phys. Rev. Lett. 91, 115001-1–115001-4 





Chertkov, M., Kolokolov, I. & Lebedev, V. 2005 Effects of surface tension on immiscible Rayleigh–Taylor turbulence. Phys. Rev. E 71, 055301-1–055301-4 





Cook, A. W. & Zhou, Y. 2002 Energy transfer in Rayleigh–Taylor instability. Phys. Rev. E 66, 026312-1–026312-12 





Dimonte, G. & Schneider, M. 2000 Density ratio dependence of Rayleigh–Taylor mixing for sustained and impulsive acceleration histories. Phys. Fluids 12, 304–321 





Ding, H., Spelt, P. D. M. and Shu, Chang 2007 Diffuse interface model for incompressible two-phase flows with large density ratios. J. Comput. Phys. 226, 2078–2095 





Di Prima, R. C. & Swinney, H. L. 1981 Hydrodynamic Instabilities and the Transition to Turbulence. eds. Swinney, H. L. & Gollup, J. P. Springer, Berlin 





Ducea, M. & Saleeby, J. 1998 A case for delamination of the deep batholithic crust beneath the Sierra Nevada, California. Int. Geology Rev. 40, 78–93 





Goncharov, V. N. 2003 Analytical model of nonlinear, single-mode, classical Rayleigh–Taylor instability at arbitrary Atwood numbers. Phys. Rev. Lett. 88(13), 134502-2–134502-4 





Korvola, T., Kupiainen, A. & Taskinen, J. 2005 Anomalous scaling for three-dimensional Cahn-Hilliard fronts. Comm. Pure Appl. Math. 58(8), 1077–1115 





Kull, H. J. 1991 Theory of the Rayleigh–Taylor instability. Phys. Rep. 206(5), 197–325 





Kundu, P. K. & Cohen, I. M. 2001 Fluids Mechanics - Second Edition Academic Press 





Jacqmin, D. 1999 Calculation of two-phase Navier–Stokes flows using Phase-Field modeling. J. Comp. Phys. 155, 96–127 





Layzer, D 1955 On the instability of superposed fluids in a gravitational field. Astrophys. J. 122, 1–12 





Lee , C.-T., Rudnick, R. L. & Brimhall Jr., G. H. 2001 Deep lithospheric dynamics beneath the Sierra Nevada during the Mesozoic and Cenozoic as inferred from xenolith petrology. Geochem. Geophys. Geosys. 2, 2001GC000152 





Liu, C. & Shen, J. 2003 A phase field model for the mixture of two incompressible fluids and its approximation by a Fourier-spectral method. Physica D 179, 211–228 





Menikoff, R., Mjolsness, R. C., Sharp, D. H. & Zemach, C. 1977 Unstable normal mode for Rayleigh–Taylor instability in viscous fluids. Phys. Fluids 20(12), 2000–2004 





L.D. Landau & E.M. Lifshitz 2000 Fluid Mechanics Volume 6 of Course of Theoretical Physics Second Edition, Revised Butterworth Heinemann 





Ramabrabhu, P. & Andrews, M. J. 2004 Experimental investigation of Rayleigh–Taylor mixing at small Atwood numbers. J. Fluid Mech. 502, 233–271 





Lord Rayleigh 1883 Investigation of the caracter of the equilibrium of an incompressible heavy fluid of variable density. Proc. London Math. Soc. 14, 170. 





Schultz, D.M., Kanak, K. M., Straka, J. M., Trapp, R. J., Gordon, B. A., Zrnic, D.´ S., Bryan, G. H., Durant, A. J., Garrett, T. J., Klein, P. K. & Lilly, D. K. 2006 The mysteries of Mammatus clouds: observations and formation mechanisms. J. Atmos. Sci. 10, 2409–2435 





Sethian, A.J. 1999 Level Set methods and Fast Marching Methods: Evolving Interfaces in Computational Geometry, Fluid Mechanics, Computer Vision and Materials Science. Cambridge University Press: Cambridge 





Sharp, D. H. 1984 An overview of Rayleigh–Taylor instability. Physica D 12, 3–18 





Smolianski, A., Haario, H. & Luukka, P. 2005 Vortex shedding behind a rising bubble and two-bubble coalescence: A numerical approach. Appl. Math. Model 29, 615–632 





Taylor, G. I. 1950 The instability of liquid surfaces when accelerated in a direction perpendicular to their planes I. Proc. R. Soc. A 201, 192–197 





Waddell, J. T., Niederhaus, C. E. & Jacobs, J. W. 2001 Experimental study of Rayleigh– Taylor instability: low Atwood number systems with single-mode initial perturbations. Phys. Fluids 13(5), 1263–1273 





Young, Y. N., Tufo, H., Dubey, A. & Rosner, R. 2001 On the miscibile Rayleigh–Taylor instability: two and three dimensions. J. Fluid Mech. 447, 377–408 





Young, Y. N.& Ham, F. E. 2006 Surface tension in incompressible Rayleigh–Taylor mixing flow. J. Turbul. 71(7), 1–23. 





Yue, P., Feng, J.J., Liu, C. & Shen, J. 2004 A diffuse-interface method for simulating twophase flows of complex fluids. J. Fluid Mech. 515, 293–317 

