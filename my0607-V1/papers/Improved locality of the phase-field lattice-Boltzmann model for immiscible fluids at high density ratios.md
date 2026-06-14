# Improved locality of the phase-field lattice-Boltzmann model for immiscible fluids at high density ratios

Abbas Fakhari,1,2 Travis Mitchell,3 Christopher Leonardi,3 and Diogo Bolster1 

1Department of Civil and Environmental Engineering and Earth Sciences, University of Notre Dame, Notre Dame, Indiana 46556, USA 

2Department of Chemical and Biomolecular Engineering, University of Pennsylvania, Pennsylvania 19104, USA 

3School of Mechanical and Mining Engineering, The University of Queensland, St Lucia, QLD 4072, Australia 

(Received 10 April 2017; revised manuscript received 10 August 2017; published 1 November 2017) 

Based on phase-field theory, we introduce a robust lattice-Boltzmann equation for modeling immiscible multiphase flows at large density and viscosity contrasts. Our approach is built by modifying the method proposed by Zu and He [Phys. Rev. E 87, 043301 (2013)] in such a way as to improve efficiency and numerical stability. In particular, we employ a different interface-tracking equation based on the so-called conservative phase-field model, a simplified equilibrium distribution that decouples pressure and velocity calculations, and a local scheme based on the hydrodynamic distribution functions for calculation of the stress tensor. In addition to two distribution functions for interface tracking and recovery of hydrodynamic properties, the only nonlocal variable in the proposed model is the phase field. Moreover, within our framework there is no need to use biased or mixed difference stencils for numerical stability and accuracy at high density ratios. This not only simplifies the implementation and efficiency of the model, but also leads to a model that is better suited to parallel implementation on distributed-memory machines. Several benchmark cases are considered to assess the efficacy of the proposed model, including the layered Poiseuille flow in a rectangular channel, Rayleigh-Taylor instability, and the rise of a Taylor bubble in a duct. The numerical results are in good agreement with available numerical and experimental data. 

DOI: 10.1103/PhysRevE.96.053301 

## I. INTRODUCTION

Numerical modeling of multiphase flows remains a challenging subject in fluid mechanics. Despite significant advances in computational fluid dynamics (CFD), configurations featuring high density ratios and/or high Reynolds numbers remain intractable [1]. In addition, the interfacial region between immiscible fluids is typically of the order of nanometers, which makes it impractical for macroscopic CFD techniques to resolve these regions. 

As an alternative, diffuse-interface modeling represents a compelling approach for the numerical simulation of multiphase flows [2]. In diffuse-interface methods, the sharp interface between different fluids is replaced with a smooth transition region across which fluid properties change continuously, thereby removing abrupt jumps and potential singularities at the interface. Several diffuse-interface models exist [1]. In this study, we will use the so-called conservative phase-field model [3], which is a subclass of diffuse-interface models [4], for interface tracking purposes. 

Given the mesoscopic nature of interfacial flows, the lattice-Boltzmann method (LBM) stands out as a natural candidate, and a well-established tool, with which the governing equations can be solved [5–7]. Historically, there are four major classes of lattice-Boltzmann (LB) models for multiphase flows. These are the chromodynamic or colorgradient model [8], the pseudopotential model [9,10], the free-energy model [11], and the mean-field model [12]. For the most part, the primitive forms of these models suffer from numerical artefacts and other restrictions such as the lack of Galilean invariance, large spurious velocities, and inability to model multiphase flows with large density contrasts [13]. Consequently, these models have been incrementally improved over the past few decades [14–23]. A good review of previous 

LB models and recent advancements in the field can be found in Refs. [21,23]. 

Despite continued progress in the LBM for studying multiphase flows, there remains plenty of scope for further improvement, particularly in situations where the density ratio, viscosity ratio, and/or the Reynolds number is high. In this study, we propose an LB model for direct numerical simulation of multiphase flows at high density ratios. Rather than using the traditional Cahn-Hilliard equation [4], the present model consists of an LB equation (LBE) for interface tracking [24] based on the conservative phase-field equation [3]. We also adapt and build on Zu and He’s [18] LBE for recovering the hydrodynamic properties. Compared with existing LB models based on advanced free-energy [15,17] or phase-field models [16,22], the proposed model is more efficient and more accurate, especially for configurations featuring large density ratios. Our model maintains stability and accuracy at high density ratios without needing to use mixed (combination of central and biased) finite-difference (FD) schemes as is the case in some advanced free-energy models [15]. Aside from the complexity in implementation, using mixed FD schemes is known to potentially compromise mass and momentum conservation [25]. The proposed model is also equipped with a multiple-relaxation time (MRT) collision operator [17,26] to enhance stability when modeling flows with both small viscosities as well as large viscosity contrasts. Moreover, the proposed LBM consists of only one nonlocal variable, i.e., the phase field, for which FDs are required to calculate its derivatives. Limiting the nonlocality of data in the model improves its parallel performance, particularly on GPUs, which in turn makes the model suitable for high-performance computing. 

We examine the accuracy of the proposed LBM by simulating three two-dimensional (2D) benchmark problems. The first is a gravity-driven, two-layer flow in a rectangular channel, which, in the context of color-gradient models, has been used to test corrections and improved accuracy and stability schemes at high density and viscosity ratios [19,27]. The second benchmark case is the well-known Rayleigh-Taylor instability, for which the results are compared with existing numerical simulations. We use the model to simulate fluid properties similar to an air-water system at a relatively high Reynolds number. The third and last benchmark is the buoyancy-driven motion of a planar Taylor bubble in a duct. For this, results of the bubble shape profile are compared to previous numerical findings and the rise velocity is compared with reported results from theoretical and numerical models as well as experimental studies, and good agreement is found. After highlighting the robustness of the model through the benchmark cases above, the computational efficiency is assessed against recent phase-field-based LB models. 

## II. MACROSCOPIC EQUATIONS

## A. Interface tracking equation

The interface-tracking equation in this study is built upon the Allen-Cahn equation [28] as opposed to the commonly used Cahn-Hilliard theory [29]. We use a specific version of the phase-field model [30] that was proposed by Sun and Beckermann [31] and reformulated in conservative form [3] to improve conservation properties. In what follows, we shall refer to this formulation as the conservative phasefield model [3]. In this model, the phase field $\phi$ assumes two extreme values, $\phi _ { \mathrm { L } }$ and $\phi _ { \mathrm { H } } ,$ , in the bulk of the light and heavy fluids, respectively. The phase-field equation governs the evolution of the interface between the two fluids [3], 

$$
\frac {\partial \phi}{\partial t} + \nabla \cdot \phi \boldsymbol {u} = \nabla \cdot M \left[ \left(\nabla \phi - \frac {\nabla \phi}{| \nabla \phi |} \frac {[ 1 - 4 (\phi - \phi_ {0}) ^ {2} ]}{\xi}\right) \right], \tag {1}
$$

where t is the time, u is the macroscopic velocity vector, M is the mobility, ξ is the interfacial thickness, and $\phi _ { 0 } = ( \phi _ { \mathrm { L } } +$ $\phi _ { \mathrm { H } } ) / 2$ indicates the location of the interface. The equilibrium profile of the phase field for an interface located at $\scriptstyle { \boldsymbol { x } } _ { 0 }$ is assumed to vary according to 

$$
\phi (\boldsymbol {x}) = \phi_ {0} \pm \frac {\phi_ {\mathrm{H}} - \phi_ {\mathrm{L}}}{2} \tanh \left(\frac {\boldsymbol {x} - \boldsymbol {x} _ {0}}{\xi / 2}\right), \tag {2}
$$

which is typically used to set the initial condition for the phase field. The sign is chosen such that the minimum value of the phase field is assigned to the light fluid. 

For example, the plus sign is used for initializing an air bubble while the minus sign is used for a liquid drop. 

## B. Navier-Stokes equations

The continuity and momentum equations for incompressible multiphase flows are given by 

$$
\frac {\partial \rho}{\partial t} + \nabla \cdot \rho \boldsymbol {u} = 0, \tag {3a}
$$

$$
\begin{array}{l} \rho \left(\frac {\partial \boldsymbol {u}}{\partial t} + \boldsymbol {u} \cdot \nabla \boldsymbol {u}\right) = - \nabla p + \nabla \cdot (\mu [ \nabla \boldsymbol {u} + (\nabla \boldsymbol {u}) ^ {T} ]) \\ + \boldsymbol {F} _ {\mathrm{s}} + \boldsymbol {F} _ {\mathrm{b}}, \tag {3b} \\ \end{array}
$$

where $\rho$ and $\mu$ are the local fluid density and viscosity, respectively, $p$ is the macroscopic pressure, $F _ { \mathrm { b } }$ is a body force, and $ { \boldsymbol { F } } _ { \mathrm { s } }$ is the surface tension force. In this work, the surface tension force takes the form 

$$
\boldsymbol {F} _ {\mathrm{s}} = \mu_ {\phi} \nabla \phi , \tag {4}
$$

where 

$$
\mu_ {\phi} = 4 \beta (\phi - \phi_ {\mathrm{L}}) (\phi - \phi_ {\mathrm{H}}) (\phi - \phi_ {0}) - \kappa \nabla^ {2} \phi (5)
$$

is the chemical potential for binary fluids. The coefficients $\beta$ and κ are related to the surface tension σ and interface thickness $\xi$ by $\beta = 1 2 \sigma / \xi$ and $\kappa = 3 \sigma \xi / 2$ . 

## III. LATTICE-BOLTZMANN EQUATIONS

## A. LBE for interface tracking

We propose the following LBE for tracking the interface between different fluids [24]: 

$$
\begin{array}{l} h _ {\alpha} (\boldsymbol {x} + \boldsymbol {e} _ {\alpha} \delta t, t + \delta t) = h _ {\alpha} (\boldsymbol {x}, t) \\ - \frac {h _ {\alpha} (\boldsymbol {x} , t) - \bar {h} _ {\alpha} ^ {\mathrm{eq}} (\boldsymbol {x} , t)}{\tau_ {\phi} + 1 / 2} + F _ {\alpha} ^ {\phi} (\boldsymbol {x}, t), \tag {6} \\ \end{array}
$$

in which the forcing term is given by 

$$
F _ {\alpha} ^ {\phi} (\boldsymbol {x}, t) = \delta t \frac {[ 1 - 4 (\phi - \phi_ {0}) ^ {2} ]}{\xi} w _ {\alpha} \boldsymbol {e} _ {\alpha} \cdot \frac {\nabla \phi}{| \nabla \phi |}, \tag {7}
$$

and $h _ { \alpha }$ is the phase-field distribution function, $\tau _ { \phi }$ is the phase-field relaxation time, and $w _ { \alpha }$ and $e _ { \alpha }$ are the weight coefficients and the mesoscopic velocity set, respectively. For the D2Q9 lattice used in this study $w _ { 0 } = 4 / 9 , w _ { 1 - 4 } = 1 / 9$ , w5 $_ { - 8 } = 1 / 3 6 [ 3 2 ]$ , and 

$$
\boldsymbol {e} _ {\alpha} = c \left\{ \begin{array}{l l} (0, 0), & \alpha = 0 \\ (\cos \theta_ {\alpha}, \sin \theta_ {\alpha}), & \theta_ {\alpha} = (\alpha - 1) \pi / 2, \quad \alpha = 1 - 4, \\ (\cos \theta_ {\alpha}, \sin \theta_ {\alpha}) \sqrt {2}, & \theta_ {\alpha} = (2 \alpha - 9) \pi / 4, \quad \alpha = 5 - 8 \end{array} \right. \tag {8}
$$

where $c = \delta x / \delta t$ , and δx and δt are the length scale and time scale of the underlying lattice, respectively. On uniform grids, it is common practice to take $\delta x = \delta t = 1$ lu (lattice units). The equilibrium phase-field distribution function is defined as 

$$
\bar {h} _ {\alpha} ^ {\mathrm{eq}} = h _ {\alpha} ^ {\mathrm{eq}} - \frac {1}{2} F _ {\alpha} ^ {\phi}, \tag {9}
$$

where $h _ { \alpha } ^ { \mathrm { e q } } = \phi \Gamma _ { \alpha }$ and 

$$
\Gamma_ {\alpha} = w _ {\alpha} \left[ 1 + \frac {\boldsymbol {e} _ {\alpha} \cdot \boldsymbol {u}}{c _ {s} ^ {2}} + \frac {(\boldsymbol {e} _ {\alpha} \cdot \boldsymbol {u}) ^ {2}}{2 c _ {s} ^ {4}} - \frac {\boldsymbol {u} \cdot \boldsymbol {u}}{2 c _ {s} ^ {2}} \right] \tag {10}
$$

is the dimensionless distribution function. The speed of sound in the system is defined as $c _ { s } = c / \sqrt { 3 }$ . The mobility M is related to the phase-field relaxation time by 

$$
M = \tau_ {\phi} c _ {s} ^ {2} \delta t. \tag {11}
$$

The phase field is updated by taking the zeroth moment of the phase-field distribution function after the streaming, or propagation, step 

$$
\phi = \sum_ {\alpha} h _ {\alpha}. \tag {12}
$$

Then the density $\rho$ is calculated by linear interpolation 

$$
\rho = \rho_ {\mathrm{L}} + (\phi - \phi_ {\mathrm{L}}) (\rho_ {\mathrm{H}} - \rho_ {\mathrm{L}}), \tag {13}
$$

where ρL and ρH are the densities of the light and heavy fluids, respectively. 

Details of the conservative phase-field model have been previously discussed in the literature [22,24]. However, it is worth noting that the current model for interface tracking is intended for immiscible (multicomponent), incompressible fluids as opposed to various other LB methods that are developed for the study of miscible (single-component) fluids [20]. A comparative study between the Cahn-Hilliard-based and Allen-Cahn-based LB models for the interface-tracking equation was conducted in Ref. [33]; however, no hydrodynamic interactions were considered. Moreover, in contrast to the D2Q9 lattice used in the present study, a less isotropic lattice (D2Q5) was used in Ref. [33]. It has been argued in Ref. [34] and shown in Ref. [35] that this type of lattice structure reduces accuracy in simulations. 

In the current phase-field formulation, we neglect a highorder temporal term for the sake of efficiency and locality of the model. Previous studies have analyzed this term (see Eq. (17) in Ref. [36]), which is related to the temporal derivative of the phase-field flux. Ren et al. [36] compared the complete formulation to the original scheme presented in Ref. [24], neglecting the nonlinear terms in their equilibrium distribution function. Meanwhile, Chai and Zhao [34] argued that using a linear equilibrium distribution function leads to additional numerical diffusion in the recovered advection-diffusion equation. Aside from the fact that having a temporal derivative in the scheme impedes the efficiency and implementation of the algorithm on parallel machines, we did not see any significant improvement in our results after inclusion of the high-order term. 

## B. LBE for hydrodynamics

In this study, we propose some improvements to the velocity-based LB approach proposed by Zu and He [18]. The LBE for hydrodynamics is defined as 

$$
g _ {\alpha} (\boldsymbol {x} + \boldsymbol {e} _ {\alpha} \delta t, t + \delta t) = g _ {\alpha} (\boldsymbol {x}, t) + \Omega_ {\alpha} (\boldsymbol {x}, t) + F _ {\alpha} (\boldsymbol {x}, t), \tag {14}
$$

where the hydrodynamic forcing is 

$$
F _ {\alpha} (\boldsymbol {x}, t) = \delta t w _ {\alpha} \frac {\boldsymbol {e} _ {\alpha} \cdot \boldsymbol {F}}{\rho c _ {s} ^ {2}}, \tag {15}
$$

and $g _ { \alpha }$ is the velocity-based distribution function for incompressible fluids with its modified equilibrium distribution given by 

$$
\bar {g} _ {\alpha} ^ {\mathrm{eq}} = g _ {\alpha} ^ {\mathrm{eq}} - \frac {1}{2} F _ {\alpha}, \tag {16}
$$

where 

$$
g _ {\alpha} ^ {\mathrm{eq}} = p ^ {*} w _ {\alpha} + (\Gamma_ {\alpha} - w _ {\alpha}), \tag {17}
$$

and $p ^ { * } = p / \rho c _ { s } ^ { 2 }$ is the normalized pressure. 

Here, the modified equilibrium distribution function in Eq. (16) is defined by subtracting half of the forcing term (according to trapezoidal rule or Crank-Nicholson discretization) from the regular equilibrium distribution function to simplify the collision step, particularly when the MRT model is used [17]. After substituting Eq. (16) into Eq. (14) and rearranging, we obtain $\begin{array} { r } { F _ { \alpha } ( \pmb { x } , t ) = \delta t ( 1 - \frac { 1 } { 2 \tau } ) w _ { \alpha } \frac { \pmb { e } _ { \alpha } \cdot \mathbf { \hat { F } } } { \rho c _ { s } ^ { 2 } } } \end{array}$ eα · F , which is consistent with Eq. (20) in Guo et al. [25] to leading order in velocity. This is also the same forcing term that was proposed and verified in Ref. [18]. Using this forcing term, Zu and He [18] derived the governing macroscopic equations. We have also examined the higher-order form of the forcing term, but did not observe any noticeable difference in the results. This is likely due to the fact that the external force F in multiphase LB models is small, i.e., $| \boldsymbol { F } | \sim \mathcal { O } ( \mathbf { M } \mathbf { a } ^ { 2 } )$ . 

The collision operator, $\Omega _ { \alpha } ,$ is defined in Eqs. (26) and (27), and the forcing term is [18] 

$$
\boldsymbol {F} = \boldsymbol {F} _ {\mathrm{s}} + \boldsymbol {F} _ {\mathrm{b}} + \boldsymbol {F} _ {\mathrm{p}} + \boldsymbol {F} _ {\mu}, \tag {18}
$$

where $\boldsymbol { F } _ { \mathrm { p } }$ and $\boldsymbol { F } _ { \mu }$ are two additional terms in the velocitybased formulation [18]. The pressure force can be written as 

$$
\boldsymbol {F} _ {\mathrm{p}} = - p ^ {*} c _ {s} ^ {2} \nabla \rho , \tag {19}
$$

and the viscous force is [see Eq. (31) for implementation] 

$$
\boldsymbol {F} _ {\mu} = \nu [ \nabla \boldsymbol {u} + (\nabla \boldsymbol {u}) ^ {T} ] \cdot \nabla \rho , \tag {20}
$$

where ν is the kinematic viscosity, which is related to the hydrodynamic relaxation time τ by 

$$
\nu = \tau c _ {s} ^ {2} \delta t. \tag {21}
$$

Given the link between the relaxation time and fluid properties, there are many ways to calculate the relaxation time from the phase field. First we discuss two of the more popular approaches, and then we propose a technique which will be shown to be more consistent and more accurate (see Sec. IV A). One approach is to use a harmonic interpolation, which favors lower values, to calculate the relaxation time [15] 

$$
\frac {1}{\tau} = \frac {1}{\tau_ {\mathrm{L}}} + (\phi - \phi_ {\mathrm{L}}) \left(\frac {1}{\tau_ {\mathrm{H}}} - \frac {1}{\tau_ {\mathrm{L}}}\right), \tag {22}
$$

where $\tau _ { \mathrm { L } }$ and τH are the relaxation rates for the light and heavy fluids, respectively. Another common approach is to use a linear interpolation, which typically favors larger values, 

$$
\tau = \tau_ {\mathrm{L}} + (\phi - \phi_ {\mathrm{L}}) (\tau_ {\mathrm{H}} - \tau_ {\mathrm{L}}). \tag {23}
$$

This, from Eq. (21), is equivalent to calculating the kinematic viscosity of the fluid using a linear interpolation. Alternatively, here we propose that the dynamic viscosity is first updated using a linear interpolation such that 

$$
\mu = \mu_ {\mathrm{L}} + (\phi - \phi_ {\mathrm{L}}) (\mu_ {\mathrm{H}} - \mu_ {\mathrm{L}}), \tag {24}
$$

where $\mu _ { \mathrm { L } }$ and $\mu _ { \mathrm { H } }$ are the viscosities of the light phase and heavy phase, respectively. After calculating the viscosity of the fluid, we can simply compute the relaxation time via 

$$
\tau = \frac {\mu}{\rho c _ {s} ^ {2}}. \tag {25}
$$

As will be shown in Sec. IV A, Eq. (25) leads to the most accurate results in LB simulations. 

The simplest form commonly used for the collision operator is the single-relaxation-time (SRT) or Bhatnagar-Gross-Krook (BGK) model, 

$$
\Omega_ {\alpha} ^ {\mathrm{BGK}} = - \frac {g _ {\alpha} - \bar {g} _ {\alpha} ^ {\mathrm{eq}}}{\tau + 1 / 2}. \tag {26}
$$

Another popular choice is the more sophisticated multirelaxation-time (MRT) model [26], which has been shown to be more accurate and more stable than the BGK model [17]: 

$$
\Omega_ {\alpha} ^ {\mathrm{MRT}} = - \mathbf {M} ^ {- 1} \hat {\mathbf {S}} \mathbf {M} \big (g _ {\alpha} - \bar {g} _ {\alpha} ^ {\mathrm{eq}} \big), \tag {27}
$$

where M is an orthogonal matrix for transforming the distribution functions from physical space into moment space [26], and Sˆ is a diagonal relaxation matrix, which may take the following form [17]: 

$$
\hat {\mathbf {S}} = \operatorname{diag} (1, 1, 1, 1, 1, 1, 1, s _ {\nu}, s _ {\nu}), \tag {28}
$$

where 

$$
s _ {\nu} = \frac {1}{\tau + 1 / 2}. \tag {29}
$$

One of the benefits of the LBM is that the deviatoric stress tensor can be locally obtained in terms of the hydrodynamic distribution function. For the BGK model, the viscous force in the i direction $( F _ { \mu , i } , i \in x , y )$ , can be obtained from 

$$
F _ {\mu , i} ^ {\mathrm{BGK}} = - \frac {\nu}{(\tau + 1 / 2) c _ {s} ^ {2} \delta t} \left[ \sum_ {\alpha} e _ {\alpha i} e _ {\alpha j} \left(g _ {\alpha} - g _ {\alpha} ^ {\mathrm{eq}}\right) \right] \frac {\partial \rho}{\partial x _ {j}}, \tag {30}
$$

while for the MRT model 

$$
\begin{array}{l} F _ {\mu , i} ^ {\mathrm{MRT}} = - \frac {\nu}{c _ {s} ^ {2} \delta t} \left[ \sum_ {\beta} e _ {\beta i} e _ {\beta j} \right. \\ \left. \times \sum_ {\alpha} \left(\mathbf {M} ^ {- 1} \hat {\mathbf {S}} \mathbf {M}\right) _ {\beta \alpha} \left(g _ {\alpha} - g _ {\alpha} ^ {\mathrm{eq}}\right) \right] \frac {\partial \rho}{\partial x _ {j}}. \tag {31} \\ \end{array}
$$

It is worth highlighting the main differences between the present model and the one put forth by Zu and He [18]. Aside from a major difference in the interface tracking LBEs, in that they use a Cahn-Hilliard type model while we use a conservative phase-field model, there are subtle, but important, differences in the hydrodynamic LBEs. The most notable is that Zu and He [18] used FDs to calculate the forcing term in Eq. (20). This adds the velocity vector to the list of nonlocal variables (i.e., distribution functions and phase field), which can impede optimal parallel computation. Another difference is that our equilibrium distribution function in Eq. (16) is not only modified to make the collision step simpler, but is also calculated all at once as a vector. This is in contrast to the way the equilibrium distribution function was written in Ref. [18], which distinguishes the particle distribution function at rest $( \alpha = 0 )$ with other directions. The third difference is that the velocity and pressure are coupled in Ref. [18] and, consequently, an iterative, predictor-corrector scheme is required to update the hydrodynamic pressure and velocity. In our model, after solving the LBE (14) using a routine collision-streaming sequence, the hydrodynamic properties are updated independently according to 

$$
p ^ {*} = \sum_ {\alpha} g _ {\alpha}, \tag {32a}
$$

$$
\boldsymbol {u} = \sum_ {\alpha} g _ {\alpha} \boldsymbol {e} _ {\alpha} + \frac {\boldsymbol {F}}{2 \rho} \delta t. \tag {32b}
$$

Note that the velocity is updated after the pressure and, as such, there is no need for the predictor-corrector scheme. Additionally, the gradient of density in Eqs. (19) and (30) can be replaced with the gradient of the phase field using Eq. (13), 

$$
\nabla \rho = (\rho_ {\mathrm{H}} - \rho_ {\mathrm{L}}) \nabla \phi , \tag {33}
$$

hence making φ the only nonlocal macroscopic variable in our multiphase LB model. This is beneficial for parallel computations on distributed memory machines. Instead of treating the $\boldsymbol { e } _ { \alpha } \cdot \nabla \phi$ terms as directional derivatives along the lattice links, as was done in Ref. [15], we compute the derivatives of the phase field in Eqs. (4), (5), (7), and (33) using second-order, isotropic centered differences [37–39], and then execute the dot product. Specifically, the gradient of the phase field in Eqs. (4), (7), and (33) is calculated by 

$$
\nabla \phi = \frac {c}{c _ {s} ^ {2} \delta x} \sum_ {\alpha} \boldsymbol {e} _ {\alpha} w _ {\alpha} \phi (\boldsymbol {x} + \boldsymbol {e} _ {\alpha} \delta t, t), \tag {34}
$$

and its Laplacian in Eq. (5) is calculated by 

$$
\nabla^ {2} \phi = \frac {2 c ^ {2}}{c _ {s} ^ {2} (\delta x) ^ {2}} \sum_ {\alpha} w _ {\alpha} [ \phi (\boldsymbol {x} + \boldsymbol {e} _ {\alpha} \delta t, t) - \phi (\boldsymbol {x}, t) ]. \tag {35}
$$

This enhances the computational efficiency of the proposed model when compared with the model proposed by Lee and Liu [15], who utilized central differences in the calculation of the forcing term for the equilibrium distribution function while employing mixed (central and biased) differences in the collision step. Using a combination of central and biased differences has been shown to compromise conservation of mass and momentum [25]. As to the collision model, we use $\Omega _ { \alpha } = \Omega _ { \alpha } ^ { \mathrm { M R T } }$ to obtain stable results at high Reynolds numbers. 

## IV. RESULTS

Before starting rigorous test studies in the following sections, it is worthwhile to discuss the choice of $\phi _ { \mathrm { L } }$ and φH in LB models. Based on our experience, it is better to use $\phi _ { \mathrm { L } } = - 0 . 5$ and $\phi _ { \mathrm { H } } = 0 . 5$ when the density ratio is 1 (or when $\rho _ { \mathrm { L } } \simeq \rho _ { \mathrm { H } } )$ . This leads to perfect symmetry in the results. For example, in Sec. 5.3.2 of Ref. [22], a slight difference in the peak-to-peak values of the lift coefficient was observed for two antisymmetric setups (see Table 2 in [22]). This is because $\phi _ { \mathrm { { L } } } = 0$ and $\phi _ { \mathrm { H } } = 1$ were used for a density ratio of 1. If we use $\phi _ { \mathrm { L } } = - 0 . 5$ and $\phi _ { \mathrm { H } } = 0 . 5$ it leads to having exactly the same values for the peak-to-peak lift coefficients. On the other hand, using $\phi _ { \mathrm { { L } } } = 0$ and $\phi _ { \mathrm { H } } = 1$ leads to more stable results when we have a noticeable density ratio. Particularly, at high density ratios, we might encounter numerical instability if $\phi _ { \mathrm { L } } = - 0 . 5$ and $\phi _ { \mathrm { H } } = 0 . 5 .$ . The reason for this is that LBM is weakly compressible, thereby the divergence of velocity is not exactly zero. As such, the advection term $( \nabla \cdot \phi \pmb { u } )$ i n Eq. (1) might not vanish in the bulk of the fluids if $\phi \neq 0 .$ . The compressibility issue is more problematic in the light $( \mathrm { i . e . }$ gas) phase than in the heavy (i.e., liquid) phase. Therefore, using $\phi _ { \mathrm { { L } } } = 0$ for the light phase causes the advection term to vanish in the bulk, and therefore enhances numerical stability. Throughout the simulations presented in this paper we use $\phi _ { \mathrm { { L } } } = 0$ and $\phi _ { \mathrm { H } } = 1$ , which gives $\phi _ { 0 } = 0 . 5$ . 

## A. Two-phase Poiseuille flow

The gravity-driven flow of a two-layer fluid in a rectangular channel is a simple but informative benchmark for assessing multiphase LB models [19,27]. Suppose we have a channel with periodic boundaries in the x direction which is bounded by two walls at the bottom $( y = 0 )$ and top $( y = L )$ . The channel is filled with a light fluid from the bottom wall to the centerline $( y = L / 2 )$ and a heavy fluid from the centerline to the top. The bulk properties of the fluids are $\rho _ { \mathrm { L } }$ and $\mu _ { \mathrm { L } }$ in the lower half of the domain and $\rho _ { \mathrm { H } }$ and $\mu _ { \mathrm { H } }$ in the upper half. A body force $\boldsymbol { F } _ { \mathrm { b } } = \rho g \hat { \boldsymbol { x } }$ , where $g$ is the magnitude of acceleration in the x direction, is applied to the entire domain. In the absence of surface tension, the Navier-Stokes equation simplifies to 

$$
\frac {d}{d y} \left(\mu \frac {d u _ {x}}{d y}\right) + \rho g = 0, \tag {36}
$$

where $u _ { x }$ is the x component of the velocity vector. The density and viscosity of the fluids are given by 

$$
\rho (y) = \frac {\rho_ {\mathrm{H}} + \rho_ {\mathrm{L}}}{2} - \frac {\rho_ {\mathrm{H}} - \rho_ {\mathrm{L}}}{2} \tanh \left(\frac {2 y - L}{\xi}\right), \tag {37a}
$$

$$
\mu (y) = \frac {\mu_ {\mathrm{H}} + \mu_ {\mathrm{L}}}{2} - \frac {\mu_ {\mathrm{H}} - \mu_ {\mathrm{L}}}{2} \tanh \left(\frac {2 y - L}{\xi}\right). \tag {37b}
$$

We can solve Eq. (36) using a second-order, compact FD scheme and consider the result as the diffuse-interface solution. 

First let us evaluate the accuracy of the interpolation scheme for updating the relaxation time in Eqs. (22)–(25) by considering two cases. Denoting the density and viscosity ratios by $\rho ^ { * } = \rho _ { \mathrm { H } } / \rho _ { \mathrm { I } }$ L and $\mu ^ { * } = \mu _ { \mathrm { H } } / \mu _ { \mathrm { L } }$ , we fix the viscosity ratio $\mu ^ { * } = 1 0 0 \mathrm { ~ } ( \tau _ { \mathrm { H } } = 0 . 5$ lu) and consider one case with $\rho ^ { * } = 1$ and another case with $\rho ^ { * } = 1 0$ . Here, the height of the channel is resolved using 64 grid points with $\xi = 4$ l u and $\mathrm { g } = 1 0 ^ { - 6 }$ lu. The results are shown in Fig. 1. The velocity profiles are normalized by the maximum velocity in the channel obtained from the FD solution. When there is no density difference in the system, as is the case in Fig. 1(a), both Eqs. (23) and (25) lead to accurate calculation of the velocity profile in the channel while Eq. (22) overestimates the expected solution. Increasing the density ratio to 10 in Fig. 1(b) reveals that using the local dynamic viscosity to update the relaxation time according to Eq. (25) gives us the most accurate solution. Therefore, we employ Eq. (25) to update the relaxation time in the simulations in the remainder of this section. The reason Eq. (22) overpredicts the velocity is that, as shown in Fig. 2, the harmonic interpolation gives too much weight to the lower viscosity in the system. Similarly for the large density ratio case Eq. (23) does not weight density and viscosity appropriately, while Eq. (25) is the most physically consistent approach. 

Next, we compare the accuracy of three different LB models in calculating this layered Poiseuille flow problem. The first method is the standard, momentum-based phase-field LBM proposed in Refs. [16,22], wherein central differences are employed to calculate the gradient of the phase field. The second model is the standard, momentum-based LBM proposed in Refs. [15,40], wherein mixed differences are employed to calculate the gradient of the phase field. And the third model is the current velocity-based phase-field LBE. 

The steady-state velocity profile obtained using the FD scheme as well as using the aforementioned LB models is shown in Fig. 3 for three different density ratios at $\mu ^ { * } = 1 0 0$ $( \tau _ { \mathrm { L } } = 0 . 5$ lu). As can be seen in Fig. 3, using the momentumbased LBM with central differences [22] deteriorates the accuracy of the results, especially at higher density ratios. The results of the model proposed in Ref. [40] and the current LB model are both in good agreement with the FD results, although the current model performs best in all cases. 

The grid dependence of the results is also shown by conducting a convergence study using different grid resolutions and measuring the $L _ { 2 }$ -norm of the numerical error according to 

$$
\| \delta u \| _ {2} = \sqrt {\frac {\sum_ {y} \left(u _ {x} - u _ {x} ^ {\mathrm{FD}}\right) ^ {2}}{\sum_ {y} \left(u _ {x} ^ {\mathrm{FD}}\right) ^ {2}}}. \tag {38}
$$

Figure 4 shows the $L _ { 2 } \mathrm { - n o r m }$ of the error versus the number of grid points in the y direction for $\rho ^ { * } = 1 0 0 0$ and $\mu ^ { * } = 1 0 0$ $( \tau _ { \mathrm { L } } = 0 . 5$ lu) at a constant Cahn number $\mathrm { C n } = \xi / L = 3 / 3 2$ . As can be seen, the current method produces the lowest error and also has the fastest convergence rate among all three models tested. It is worth noting that we use link bounceback at the bottom and top boundaries, which in the case of single-phase flows would result in a second-order convergence rate [41]. For a two-layer Poiseuille flow, however, the relaxation times at the bottom and top of the domain are different, which leads to a shift in the actual wall location as we refine the mesh. In other words, the effective location of the wall, where the velocity is zero, is not necessarily aligned halfway between the boundary nodes and the adjacent fluid nodes [42]. That might explain why we do not observe a second-order rate of convergence. 

## B. Rayleigh-Taylor instability

The instability created when a heavy fluid layer lies above a lighter fluid within a gravitational field g is a common multiphase flow benchmark problem [12,18,36,43]. Perturbing the interface causes an instability called the Rayleigh-Taylor instability, whereby the heavy fluid penetrates into the lower 

![](images/2ff37427b60874a7e602ab4467deaa9d5a3a1262a232b4e2d2d533524881cfc7.jpg)



(a)


![](images/0113ac2f4e35250480ace3b7fc28f49136c5e393d01d5a86e5878b897fc8f1b4.jpg)



(b)



FIG. 1. Effect of using different interpolation schemes for updating the relaxation time in the layered Poiseuille flow at $\mu ^ { * } = 1 0 0$ $( \tau _ { \mathrm { H } } = 0 . 5$ lu) with (a) $\rho ^ { * } = 1$ and (b) $\rho ^ { * } = 1 0 .$ . The FD solution is shown by the solid back line, the black cross symbols represent the use of the harmonic interpolation in $\operatorname { E q . } \left( 2 2 \right)$ , the blue symbols (I) represent the use of the linear interpolation in Eq. (23), and the red circles represent the use of the dynamic viscosity to update the relaxation time according to Eq. (25).


layer. This problem has been widely studied due to its relevance in numerous natural and engineering phenomena [44]. Our setup consists of a domain $[ 0 , L ] \times [ - 2 L , 2 L ]$ , in which wall boundaries restrict the vertical direction, and periodic conditions are applied horizontally. The top of the domain consists of the heavy fluid $( \rho _ { \mathrm { H } } , \mu _ { \mathrm { H } } )$ , while the light fluid $( \rho _ { \mathrm { L } }$ , $\mu _ { \mathrm { L } } )$ is situated below this. The initial interface position is a flat line at $y = 0$ , which is then perturbed by a cosine function 

![](images/6bc05f784a1b7be5a8d4b1b1298bdd055035ef60ce18d81190db6a0175e2149b.jpg)



FIG. 2. The behavior of the relaxation time using different interpolation schemes $( \rho ^ { * } = 1 0$ and $\mu ^ { * } = 1 0 0 )$ .


$$
\boldsymbol {x} _ {0} = 0. 1 L \times \cos (2 \pi x / L). \tag {39}
$$

The phase field is then initialized according to 

$$
\phi (\boldsymbol {x}) = \phi_ {0} + \frac {\phi_ {\mathrm{H}} - \phi_ {\mathrm{L}}}{2} \tanh \left(\frac {| \boldsymbol {x} - \boldsymbol {x} _ {0} | _ {\perp}}{\xi / 2}\right), \tag {40}
$$

where $\left| x - \pmb { x } _ { 0 } \right| _ { \perp }$ is the signed distance from any grid point to $\scriptstyle { \mathbf { 1 } } _ { 0 } .$ It should be noted that two additional initialization strategies were tested. In the first, a sharp interface was used, while in the second the phase field was defined as 

$$
\phi = \phi_ {0} + \frac {\phi_ {\mathrm{H}} - \phi_ {\mathrm{L}}}{2} \tanh \left(\frac {y - y _ {0}}{\xi / 2}\right), \tag {41}
$$

where $y _ { 0 }$ is the initial height of the interface. In all of the cases considered, we did not observe any significant difference between the final results. 

In order to compare the results of the current model to others existing in the literature, the dimensionless Atwood and Reynolds numbers are defined as 

$$
\mathrm{At} = \frac {\rho_ {\mathrm{H}} - \rho_ {\mathrm{L}}}{\rho_ {\mathrm{H}} + \rho_ {\mathrm{L}}}, \tag {42}
$$

$$
\mathrm{Re} = \frac {\rho_ {\mathrm{H}} U _ {0} L}{\mu_ {\mathrm{H}}}, \tag {43}
$$

where $U _ { 0 } = \sqrt { g L }$ is the reference velocity scale. In order to uniquely define all physical quantities, we need two additional dimensionless parameters, namely the viscosity ratio $\mu ^ { * }$ and the capillary number 

![](images/730d80ea54eeb6ce0b92d6c53fbdd2b5585b00d98b7039d9540e261f1aafd135.jpg)



(a)


![](images/f4114facd3cb90d61347f3725d1194ab53fde2b7f4c1e19b51c8d1e85a59a217.jpg)



(b)


![](images/bba31d6fb4ab97be17607b798421d1ae6e68ad92a70f43c6bf9b9f4f3554521f.jpg)



(c)



FIG. 3. Normalized velocity profile for the two-phase Poiseuille flow at $\mu ^ { * } = 1 0 0 \left( \tau _ { \mathrm { L } } = 0 . 5 \right.$ lu) and density ratio of (a) 10, (b) 100, and (c) 1000. The FD solutions are shown by solid back lines, the LB results of Ref. [22] are labeled “central” and shown by green squares, the LB results of Ref. [40] are labeled “mixed” and shown by blue triangles, and the current results are labeled “current” and shown by red circles.


$$
\mathrm{Ca} = \frac {\mu_ {\mathrm{H}} U _ {0}}{\sigma}. \tag {44}
$$

Additionally, the numerical Péclet number is defined as 

$$
\mathrm{Pe} = \frac {U _ {0} L}{M}. \tag {45}
$$

For verification purposes, the computational parameters are specified consistently with Ref. [36]. A reference length of 256 lu is taken and a reference time is specified as $t _ { 0 } =$ $\sqrt { L / g \mathrm { A t } } = 1 6 0 0 0$ lu, such that $t ^ { * } = t / t _ { 0 }$ is dimensionless time. Other parameters are $\mu ^ { * } = 1 , \mathrm { C a } = 0 . 2 6$ , and $\xi = 5$ lu. In this section, the relaxation time is calculated through a linear update according to Eq. (23). It was found that the very low simulation viscosities used to obtain high Reynolds numbers caused numerical instabilities if the viscosity update in Eq. (25) was used, suggesting that the benefits we identified in the previous section may also come with potential shortcomings. 

![](images/91960e2086c23997edf67e574f6cc9adc9577410d9aeaf8fc0ce0dd47620f078.jpg)



FIG. 4. Convergence study for the layered Poiseuille flow at $\rho ^ { * } =$ 1000 and $\mu ^ { * } = 1 0 0 \left( \tau _ { \mathrm { L } } = 0 . 5 \ln \right)$ ).


The time evolution of the Rayleigh-Taylor instability for $\mathrm { A t } = 0 . 5 0 0$ is depicted in Fig. 5. Here, $ { \mathrm { R e } } = 3 0 0 0$ and $\mathrm { P e } = 1 0 0 0$ are chosen to match the flow regime found in previous studies [18,36,43]. The heavy fluid is observed to symmetrically penetrate the lighter fluid, prior to the generation of counter-rotating vortices. The notable instability of these vortices can be seen as they shed into a wake region behind the heavy liquid front. 

The results of the widely used momentum-based LBM using either isotropic central difference [22] or mixed difference [40] schemes, along with the benchmark data from previous studies [18,36,43], are compared with the currently proposed model. Figure 6 shows the dimensionless positions of the bubble and liquid fronts versus dimensionless time. It is clear that the results obtained using the current model agree well with previously published data. The results obtained using different LB schemes are also in close agreement with each other, suggesting reasonable accuracy for the case where the density ratio is relatively low. 

Currently, there exists few studies which analyze the high density ratio Rayleigh-Taylor instability using phase-field theory. Reference [36] looked to qualitatively assess this problem for a moderate density ratio at a high Reynolds number $( \rho ^ { * } = 9 9 $ and $ { \mathrm { R e } } = 3 0 0 0 )$ , and Ref. [45] presented results for a high density ratio at a moderate Reynolds number $( \rho ^ { * } = 1 0 0 0$ and $ { \mathrm { R e } } { } = 2 0 0 )$ . Here we look to use the model we have proposed to capture both a high density ratio and a high Reynolds number flow $( \rho ^ { * } = 1 0 0 0 \mathrm { a n d R e } = 3 0 0 0 )$ . The viscosity ratio for the simulation is 100 to match a system similar to air-water, and the capillary number is 0.44. Figure 7 shows the time evolution where the model is seen to stably capture the propagation of both the high and low density fronts. This is particularly promising as the model proposed in Ref. [36] with an MRT scheme was reportedly not able to capture the situation investigated here. 

![](images/b5e4d4fd8d82f886dac488726b5dcc26c674ca8814a5b39fefa2a9c212ec74c2.jpg)



FIG. 5. The evolution of a single mode Rayleigh-Taylor instability at At = 0.500 $( \rho ^ { * } = 3 )$ , Re = 3000, μ∗ = 1, Ca = 0.26, and $\mathrm { P e } = 1 0 0 0$ .


![](images/658629c16125e6594182cda4ccc1f37b6e831de214077a812370fb233b184ae6.jpg)



(a)


![](images/081846d0032c41d49d324840153ad4872f127b9dac1303001359880d960abec3.jpg)



FIG. 6. Time evolution of the Rayleigh-Taylor instability at ${ \mathrm { A t } } = 0 . 5 0 0 ~ ( \rho ^ { * } = 3 )$ , Re = 3000, μ∗ = 1, Ca = 0.26, and $\mathrm { P e } = 1 0 0 0$ for (a) bubble front position and the (b) liquid front position. Comparative results were extracted from Refs. [18,36,43]. The value of yi defines the interface position at (a) x = 0 and (b) $x = L / 2$ during the simulation.


![](images/8e6f7d2d0da9ed0d6524ef38d8f4a0dbcb9e50bd9db89e421f43fe51e077c577.jpg)



FIG. 7. The evolution of a single mode Rayleigh-Taylor instability at At = 0.998 $( \rho ^ { * } = 1 0 0 0 )$ , Re = 3000, $\mu ^ { * } = 1 0 0 .$ , Ca = 0.44, and $\mathrm { P e } = 1 0 0 0 .$ .


## C. Planar Taylor bubble

There is significant practical interest in the motion of long bubbles due to their relation to modeling the flow of liquid slugs commonly seen in the oil and gas industry, nuclear reactors, and chemical engineering. The variable rate of gas flow within a confined geometry such as a pipe or channel can lead to a number of characteristic interface topologies, commonly reported as flow regimes. At low gas flow rates a bubbly flow occurs where a large number of small, mostly spherical bubbles rise through the fluid domain. Higher gas rates typically result in an increased rate of bubble coalescence, eventually forming a reduced number of larger bubbles that occupy nearly the entire cross section of the domain. As these bubbles propagate, they form an elongated bullet shape due to the wall confinement and are often referred to as a Taylor bubble. The Taylor bubbles are separated by liquid slugs, within which smaller gas bubbles may still be observed. 

As a single Taylor bubble rises through a dense fluid, the viscous, inertial, and interfacial forces acting on it can have significant influence on both its shape and its rise velocity. The shape of the Taylor bubble can be characterized by a rounded leading edge followed by an almost cylindrical or rectangular body depending on the flow domain. The trailing edge shape depends strongly on the flow condition and liquid properties with flat, rounded, indented, or jagged profiles reported in the literature. Flow separation in the wake can also be expected for Taylor bubbles at moderate Reynolds numbers, with the transition to separation observed at a Reynolds number between 13.4 and 32.6 for tubular flows [46]. An increasing Reynolds number also indicates a transition to an inertial regime, in which viscous and interfacial forces have a lesser, or in some cases negligible, impact on the flow dynamics. 

In this section, the proposed model is used to simulate the rise of a planar Taylor bubble through stagnant fluid in an inertial regime. This case has been studied theoretically [47,48], numerically [49,50], and experimentally [48] by a number of authors, and a summary of these works can be found in [51]. Table I reproduces the findings of Ref. [51] to present the propagation speeds expected for this benchmark case. It is noted here that $V _ { \infty } ^ { * }$ is the dimensionless rise velocity of the bubble, commonly referred to as the Froude number, 


TABLE I. Planar Taylor bubble results for dimensionless rise velocity $( V _ { \infty } ^ { * } )$ with negligible surface tension, re-created from the works of Ref. [51].


<table><tr><td>Authors</td><td>Approach</td><td><eq>V_{\infty}^{*}</eq></td></tr><tr><td>Birkhoff and Carter (1957) [47]</td><td>inviscid theory</td><td>0.23</td></tr><tr><td>Watson (in [47])</td><td>experimental</td><td>0.22–0.23</td></tr><tr><td>Griffith (in [47])</td><td>experimental</td><td>0.23</td></tr><tr><td>Collins (1964) [48]</td><td>inviscid theory</td><td>0.23</td></tr><tr><td>Collins (1964) [48]</td><td>experimental</td><td>0.22–0.23</td></tr><tr><td>Mao and Dukler (1990) [49]</td><td>numerical</td><td>0.22</td></tr><tr><td>Ha-Ngoc and Fabre (2004) [50]</td><td>numerical</td><td>0.22</td></tr></table>

$$
\mathrm{Fr} = V _ {\infty} ^ {*} = \frac {u _ {r}}{U _ {0}}, \tag {46}
$$

where $u _ { r }$ is the rise velocity, $U _ { 0 } = \sqrt { g L }$ is the characteristic velocity, and L is the length of the channel in the $y$ direction. The results presented here are determined under the assumption of small surface tension and a bubble rise Reynolds number, 

$$
\mathrm{Re} _ {r} = \frac {\rho_ {\mathrm{H}} u _ {r} L}{\mu_ {\mathrm{H}}} \geqslant 1 0 0. \tag {47}
$$

In addition to the results provided in Table I, Ha-Ngoc and Fabre [50] provided the numerical results for the bubble Froude number as a function of the Eötvös number, 

$$
\mathrm{Eo} = \frac {(\rho_ {\mathrm{H}} - \rho_ {\mathrm{L}}) g L ^ {2}}{\sigma}. \tag {48}
$$

They were able to conclude that at low surface tensions, the Froude number of the Taylor bubble was independent of the Eötvös number, tending towards $\mathrm { F r } = 0 . 2 2$ . Additionally, the authors managed to predict the Taylor bubble shapes using the boundary element method for ${ \mathrm { E o } } = 1 0 , 1 0 0 { \mathrm { \Omega } }$ , and 1000. In this work, we look to compare stabilized interface profiles, as well as the bubble rise velocity, with $\mathrm { E o } = 1 0 0$ using the proposed LBM. 

Figure 8 indicates the problem construction used to analyze the planar Taylor bubble. Here, a rectangular gas region with a semicircular front is initialized and a gravitational acceleration is applied acting against the direction of curvature. The bubble then propagates along the channel, transported by the liquid movement and gravitational effects. The simulation domain was defined as $[ 1 0 L \times L ]$ , with L being equal to 259 lu and an outer layer of nodes flagged as solid surrounding this with full bounceback applied. The fluid properties are $\mu ^ { * } = 1 0 0$ and $\rho ^ { * } = 1 0 0 0 \ ( \mathrm { A t } = 0 . 9 9 8 )$ , typical of an air-water system, and the gravitational force $\boldsymbol { F } _ { \mathrm { b } } = - \rho \mathrm { g } \hat { x }$ is applied to the entire fluid. A reference time was defined as per Sec. IV B, with $t _ { 0 } =$ 24 000. To match the flow conditions described in Ref. [50], we specify $ { \mathrm { R e } } _ { r } = 2 0 0$ using the expected $\mathrm { F r } = 0 . 2 2$ and Eo 100. 

Figure 9 shows the time evolution of the Taylor bubble with the grey region representing fluid where $\phi = 1$ and the white region where $\phi = 0 .$ . The expansion of the liquid film as it passes the end of the bubble induces a recirculating wake region that causes extension of the trailing edge and is capable of liberating smaller bubbles from the initial gas region. Here the shearing force from the heavy fluid and the recirculation of the falling liquid layer was sufficient to cause a continuous breakup and coalescence-type behavior in the bubble wake. This behavior in the wake region was observed to have no significant impact on the shape of the Taylor bubble front or the rise velocity. 

Six contours of the stabilized shape profile found at the conclusion of the simulation are displayed in Fig. 10. Here we highlight that a diffuse-interface model was used for these simulations and, as such, contours of the phase field are graphed for comparison with the sharp-interface result in Ref. [50]. It is seen that the center of the diffuse interface produces a thinner Taylor bubble, but the curvature of the outer regions of the diffuse layer appear to match quite well with the sharp interface solution. 

The steady rise velocity was found by tracking the position of the bubble front, where $\phi = \phi _ { 0 } = 0 . 5$ , at intervals of $0 . 5 t _ { 0 }$ throughout the simulation. A linear regression was then performed using the final five data points with consistency checked against the remainder. The progression of the bubble front in intervals of $t _ { 0 }$ is displayed in Fig. 11 in comparison to the regression used to determine the velocity. This was additionally verified by assessing the average velocity of the entire gas bubble, as well as the instantaneous velocity at the front of the bubble where $\phi = 0 . 5$ . In the test case, a bubble Froude number of 0.217 was observed, which very closely matches the expected range of 0.22–0.23 from Table I. 

Overall, the results using the proposed LB model were shown to agree well with those based on the sharp-interface model as well as with experimental data in terms of the planar Taylor bubble shape and rise velocity. As was shown in Fig. 1, using a linear interpolation via Eq. (23), although not as accurate as using Eq. (25), is more accurate than using a harmonic interpolation via Eq. (22). Additionally, the work in both Secs. IV B and IV C has indicated that using the linear interpolation improves numerical stability in comparison to the dynamic viscosity update in Eq. (25), particularly when the relaxation time is small. Therefore, for the planar Taylor bubble results presented in this section, a linear interpolation of the relaxation time was used according to Eq. (23). It is noted that instability arises if the relaxation time is updated through the local viscosity via Eq. (25), again highlighting limitations of this approach. 

![](images/48fae486276d31e1c63d185f5572d5c89fbca73b932a6fb0a6728efb2de23b7f.jpg)



FIG. 8. Domain schematic of the slug flow tests for the Taylor bubble rise. The fluid domain size is $1 0 L \times L$ , and the initial bubble size is $3 L \times 4 L / 5$ .


![](images/bb3e2801675174ef52676ff9ca2ed37a043b9a5e65b66beee4c38987787333d8.jpg)



FIG. 9. The time evolution of the planar Taylor bubble with snapshots taken at $t ^ { * } = 0 , 4 , 8 , 1 2 , 1 6 , 2 0$ . The fluid properties are defined by $\rho ^ { * } = 1 0 0 0$ and $\mu ^ { * } = 1 0 0$ , while the flow condition is specified through $ { \mathrm { R e } } _ { r } = 2 0 0$ and $\mathrm { E o } = 1 0 0 \mathrm { . }$ .


## D. Computational efficiency

For many applications of scientific and industrial relevance, the number of lattice sites is often substantial. Hence, efficient parallel performance is essential. The previous model, presented in Ref. [22], uses a stencil consisting of a single layer of neighboring cells, but we have shown in Sec. IV A that its accuracy deteriorates around the interface. In order to capture the interfacial dynamics more accurately, the model in Ref. [40] uses mixed differences, which requires two auxiliary lattice sites (two ghost cells) in each direction. The current model requires only a single stencil and is able to model the flow field at the liquid-gas interface with a high level of accuracy. In the following we aim to probe the computational efficiency of these models. 

To investigate the performance of the aforementioned phase-field models, we implemented a stationary bubble test on a square domain $L \times L$ with a bubble radius of $R = L / 4$ . Taking $L = 5 1 2 0$ lu resulted in a test domain of approximately 26 million cells. With this setup, we analyzed the strong scalability of the methods. The domain was divided into smaller portions, inducing sublinear parallelism. This is in contrast to “weak” scaling, where the size of the mesh is kept proportional to the number of processors. The simulations were completed using the open-source TCLB solver [52] on the Prometheus cluster at Cyfronet, Krakow. This is equipped with CPU nodes fitted with two 12-core Intel Xeon E5-2680 v3 processors and eight additional GPU nodes with two nVidia Tesla K40 cards on each. 

![](images/2191c3f9790da562b9205cefd74623dd941da3a3848dddd32834da003db593cc.jpg)



FIG. 10. Contours of the phase field for a Taylor bubble at $t ^ { * } = 2 0$ with $\mathrm { E o } = 1 0 0$ and $ { \mathrm { R e } } _ { r } = 2 0 0$ . The results from Ref. [50] were supplied by Dr. J. Fabre allowing for the current LBM outputs to be compared with the profile obtained using the boundary element method. The values $x _ { i }$ and yi are used to define the interface location with respect to the bubble nose located at $^ { ( 3 , 0 ) }$ .


![](images/bd2b5bdd62e69e35c74d7fde3f383c03d4b56aafb8053e3e2dd8fbae9bd50f71.jpg)



FIG. 11. Rise of the Taylor bubble front $( \phi = \phi _ { 0 } = 0 . 5 )$ vs time, where Eo 100 and $ { \mathrm { R e } } _ { r } = 2 0 0$ .


Figure 12 shows the performance of the TCLB solver for the various models implemented on a CPU architecture. It is clear that the compared methods have a similar performance, with the speed per node generally decreasing for higher numbers of utilized cores. For the current model, as in all the previous benchmark simulations, we implemented the MRT operator in the calculation of both the stress tensor and the collision step. In the computational results, it is seen that the MRT formulation of the proposed model increases the computation required per node beyond the memory reduction benefit on the CPU architecture for this simulation domain. 

Figure 13 shows the scaling of the TCLB code on a GPU architecture for the compared methods. As expected, the parallel performance can be seen to marginally decrease with core saturation (2880 CUDA cores per nVidia Tesla K40). It is on this parallel architecture that the benefit of the reduced stencil is realized. The difference between the CPU and GPU performances is well explained by the memory access patterns. Here it is clear that the CPU performance is computation bound, whereas on the parallel architecture the performance is bound by memory. As a result of this, the mixed difference approach, which requires a larger computational stencil than the central difference and current model, sees a significant reduction in computational efficiency. In comparison on the CPU architecture, where the memory access speed is higher, there appears no significant distinction between the models for this simulation domain. 

It is often difficult to objectively compare the performance of GPU and CPU codes [53]. In the presented tests, the speed of a single CUDA core is substantially lower than a single CPU core. However, there are 2880 CUDA cores on a single GPU allowing it to vastly outperform a single CPU processor. One technique to compare these computing architectures is to look at the energy efficiency of the computation. The power consumption of a GPU node, consisting of two K40 processors, was measured at 490.5 W. Whereas for a CPU node, with two processors consisting of 24 total cores, the power consumption was 277.0 W. However, this increase in power gives approximately a factor of 4 increase in lattice updates for the proposed model, outweighing the increased energy cost. 

## V. SUMMARY AND CONCLUSION

In this work, we have proposed a multiphase LBM for simulation of immiscible fluids at high density ratios. The conservative phase-field LBE was used to track the interface dynamics while a robust, velocity-based LBE was proposed to capture the hydrodynamics. Since the only nonlocal, macroscopic variable in the current model is the phase field, the proposed LB algorithm is well suited to high-performance computing on massively parallel machines. It was shown that using isotropic central differences, which reduces the computational cost, is adequate for achieving a numerically stable and accurate LBM for multiphase flows at high density and viscosity ratios relieving us from the computational cost and complexity of using biased and mixed finite difference schemes. Additionally, numerous update rules were provided for the relaxation time across the interface. We found that a simple linear interpolation provides an increased stability, while updating the relaxation time via the dynamic viscosity was the most accurate approach. 

![](images/69b7d33774ce6f2a5c29c949c2abffb8d5b8caac60c0a04844e87aae5353240f.jpg)



FIG. 12. The strong scalability of the models on the TCLB solver on CPU.


![](images/cc2e6da55d896dfc8a9dba9b164dca7ae4d6ef9a8258d2f469bc39ad99d942a7.jpg)



FIG. 13. The strong scalability of the models on the TCLB solver on GPU.




[18] Y. Q. Zu and S. He, Phase-field-based lattice Boltzmann model for incompressible binary fluid systems with density and viscosity contrasts, Phys. Rev. E 87, 043301 (2013). 



The proposed LB model was tested against the classical layered Poiseuille flow, where it was able to accurately capture the momentum equation at the phase interface. The velocitybased formulation was able to eliminate the nonphysical velocity oscillations at the interface that are observed when using a momentum-based formulation with central differences. The stability and accuracy of the model in capturing complex interface topologies was assessed through the Rayleigh-Taylor instability. The results from this were shown to closely match the results from the available numerical data in the literature. Furthermore, the model was used to examine the rise of a planar Taylor bubble, and the stabilized interface profile of the bubble front was shown to be in good agreement with previous numerical solutions. The terminal rise velocity was also found to approach the expected range from the available analytical, numerical, and experimental studies. 



[19] Y. Ba, H. Liu, Q. Li, Q. Kang, and J. Sun, Multiple-relaxationtime color-gradient lattice Boltzmann model for simulating twophase flows with high density ratio, Phys. Rev. E 94, 023310 (2016). 





[20] D. Lycett-Brown and K. H. Luo, Cascaded lattice boltzmann method with improved forcing scheme for large-density-ratio multiphase flow at high reynolds and weber numbers, Phys. Rev. E 94, 053313 (2016). 



Overall, the proposed LB formulation allows for accurate and efficient recovery of the hydrodynamics at high density ratios while improving the locality of the LBE method, allowing it to better exploit the inherent computational efficiency of the LBM. Future work is under way to extend the model to three-dimensional lattice structures. Furthermore, a rigorous study of viscosity induced errors, improvements in Galilean invariance, and stability and accuracy analysis of the interface-tracking equation warrants further investigation. 



[21] Q. Li, K. H. Luo, Q. J. Kang, Y. L. He, Q. Chen, and Q. Liu, Lattice Boltzmann methods for multiphase flow and phase-change heat transfer, Prog. Energy Combust. Sci. 52, 62 (2016). 



## ACKNOWLEDGMENTS



[22] A. Fakhari and D. Bolster, Diffuse interface modeling of three-phase contact line dynamics on curved boundaries: A lattice Boltzmann model for large density and viscosity ratios, J. Comput. Phys. 334, 620 (2017). 



We would like to thank Ł. Łaniewski-Wołłk for assisting with the scaling test results of the compared models as well as the authors of the open-source lattice-Boltzmann framework, TCLB (https://github.com/CFD-GO/TCLB). A.F. and D.B. were supported by, or in part by, NSF Grants No. EAR-1351625, No. EAR-1417264, and No. EAR-1446236. The research was supported in part by PL-Grid Infrastructure (Grant No. clb2017). T.M. would also like to acknowledge the support of the Australian Government Research Training Program Scholarship during the development of this work. 



[23] S. Leclaire, A. Parmigiani, O. Malaspinas, B. Chopard, and J. Latt, Generalized three-dimensional lattice Boltzmann colorgradient method for immiscible two-phase pore-scale imbibition and drainage in porous media, Phys. Rev. E 95, 033306 (2017). 





[24] M. Geier, A. Fakhari, and T. Lee, Conservative phase-field lattice Boltzmann model for interface tracking equation, Phys. Rev. E 91, 063309 (2015). 





[1] A. Prosperetti and G. Tryggvason, Computational Methods for Multiphase Flow (Cambridge University Press, Cambridge, England, 2009). 





[25] Z. Guo, C. Zheng, and B. Shi, Force imbalance in lattice Boltzmann equation for two-phase flows, Phys. Rev. E 83, 036707 (2011). 





[2] D. M. Anderson, G. B. McFadden, and A. A. Wheeler, Diffuseinterface methods in fluid mechanics, Annu. Rev. Fluid Mech. 30, 139 (1998). 





[26] P. Lallemand and L.-S. Luo, Theory of the lattice Boltzmann method: Dispersion, dissipation, isotropy, Galilean invariance, and stability, Phys. Rev. E 61, 6546 (2000). 





[3] P.-H. Chiu and Y.-T. Lin, A conservative phase field method for solving incompressible two-phase flows, J. Comput. Phys. 230, 185 (2011). 





[27] S. Leclaire, N. Pellerin, M. Reggio, and J. Trépanier, A multiphase lattice Boltzmann method for simulating immiscible liquid-liquid interface dynamics, Appl. Math. Model. 40, 6376 (2016). 





[4] D. Jacqmin, Calculation of two-phase Navier-Stokes flows using phase-field modeling, J. Comput. Phys. 155, 96 (1999). 





[28] S. M. Allen and J. W. Cahn, Mechanisms of phase transformations within the miscibility gap of Fe-Rich Fe-Al alloys, Acta Metall. 24, 425 (1976). 





[5] S. Chen and G. Doolen, Lattice Boltzmann method for fluid flows, Annu. Rev. Fluid Mech. 30, 329 (1998). 





[29] J. W. Cahn and J. E. Hilliard, Free energy of a nonuniform system. i. interfacial free energy, J. Chem. Phys. 28, 258 (1958). 





[6] M. C. Sukop and D. T. Thorne, Lattice Boltzmann Modeling: An Introduction for Geoscientists and Engineers (Springer, Berlin, 2006). 





[30] R. Folch, J. Casademunt, A. Hernández-Machado, and L. Ramírez-Piscina, Phase-field model for Hele-Shaw flows with arbitrary viscosity contrast. i. theoretical approach, Phys. Rev. E 60, 1724 (1999). 





[7] H. Huang, M. C. Sukop, and X.-Y. Lu, Multiphase Lattice Boltzmann Methods: Theory and Application (Wiley-Blackwell, Chichester, West Sussex, UK, 2015). 





[31] Y. Sun and C. Beckermann, Sharp interface tracking using the phase-field equation, J. Comput. Phys. 220, 626 (2007). 





[8] A. K. Gunstensen, D. H. Rothman, S. Zaleski, and G. Zanetti, Lattice Boltzmann model of immiscible fluids, Phys. Rev. A 43, 4320 (1991). 





[32] X. He and L.-S. Luo, A priori derivation of the lattice Boltzmann equation, Phys. Rev. E 55, R6333 (1997). 





[9] X. Shan and H. Chen, Lattice Boltzmann model for simulating flows with multiple phases and components, Phys. Rev. E 47, 1815 (1993). 





[33] H. L. Wang, Z. H. Chai, B. C. Shi, and H. Liang, Comparative study of the lattice Boltzmann models for Allen-Cahn and Cahn-Hilliard equations, Phys. Rev. E 94, 033304 (2016). 





[34] Z. Chai and T. S. Zhao, Lattice Boltzmann model for the convection-diffusion equation, Phys. Rev. E 87, 063309 (2013). 





[10] X. Shan and G. D. Doolen, Multicomponent lattice-Boltzmann model with interparticle interaction, J. Stat. Phys. 81, 379 (1995). 





[35] A. Fakhari, M. Geier, and D. Bolster, A simple phase-field model for interface tracking in three dimensions, Comput. Math. Appl. (2016), doi:10.1016/j.camwa.2016.08.021. 





[11] M. R. Swift, W. R. Osborn, and J. M. Yeomans, Lattice Boltzmann Simulation of Nonideal Fluids, Phys. Rev. Lett. 75, 830 (1995). 





[12] X. He, S. Chen, and R. Zhang, A lattice Boltzmann scheme for incompressible multiphase flow and its application in simulation of Rayleigh-Taylor instability, J. Comput. Phys. 152, 642 (1999). 





[36] F. Ren, B. Song, M. C. Sukop, and H. Hu, Improved lattice Boltzmann modeling of binary flow based on the conservative Allen-Cahn equation, Phys. Rev. E 94, 023311 (2016). 





[13] X. He and G. D. Doolen, Thermodynamic foundations of kinetic theory and lattice Boltzmann models for multiphase flows, J. Stat. Phys. 107, 309 (2002). 





[37] D. Qian, Bubble motion, deformation, and breakup in stirred tanks, Ph.D. thesis, Clarkson University, 2003. 





[14] T. Inamuro, T. Ogata, S. Tajima, and N. Konishi, A lattice Boltzmann method for incompressible two-phase flows with large density differences, J. Comput. Phys. 198, 628 (2004). 





[38] A. Kumar, Isotropic finite-differences, J. Comput. Phys. 201, 109 (2004). 





[15] T. Lee and L. Liu, Lattice Boltzmann simulations of micronscale drop impact on dry surfaces, J. Comput. Phys. 229, 8045 (2010). 





[39] K. K. Mattila, L. A. H. Júnior, and P. C. Philippi, Highaccuracy approximation of high-rank dderivatives: Isotropic finite differences based on lattice-Boltzmann stencils, Sci. World J. 2014, 1 (2014). 





[16] A. Fakhari and M. H. Rahimian, Phase-field modeling by the method of lattice Boltzmann equations, Phys. Rev. E 81, 036707 (2010). 





[40] A. Fakhari, M. Geier, and T. Lee, A mass-conserving lattice Boltzmann method with dynamic grid refinement for immiscible two-phase flows, J. Comput. Phys. 315, 434 (2016). 





[17] A. Fakhari and T. Lee, Multiple-relaxation-time lattice Boltzmann method for immiscible fluids at high Reynolds numbers, Phys. Rev. E 87, 023304 (2013). 





[41] A. Fakhari and T. Lee, Numerics of the lattice Boltzmann method on nonuniform grids: Standard LBM and finite-difference LBM, Comput. Fluids 107, 205 (2015). 





[42] I. Ginzburg and D. d’Humières, Multireflection boundary conditions for lattice Boltzmann models, Phys. Rev. E 68, 066614 (2003). 





[43] Q. Li, K. H. Luo, Y. J. Gao, and Y. L. He, Additional interfacial force in lattice Boltzmann models for incompressible multiphase flows, Phys. Rev. E 85, 026704 (2012). 





[44] D. Sharp, An overview of Rayleigh-Taylor instability, Physica D 12, 3 (1984). 





[45] J. Shao and C. Shu, A hybrid phase field multiple relaxation time lattice Boltzmann method for the incompressible multipahse flow with large density contrast, Int. J. Numer. Methods Fluids 77, 526 (2015). 





[46] C. W. Kang, S. Quan, and J. Lou, Numerical study of a Taylor bubble rising in stagnant liquids, Phys. Rev. E 81, 066308 (2010). 





[47] G. Birkhoff and D. Carter, Rising plane bubbles, J. Math. Phys. 6, 769 (1957). 





[48] R. Collins, A simple model of plane gas bubble in a finite liquid, J. Fluid Mech. 22, 763 (1964). 





[49] Z. Mao and A. Dukler, The motion of Taylor bubbles in vertical tubes. I. A numerical simulation for the shape and rise velocity of Taylor bubbles in stagnant and flowing liquid, J. Comput. Phys. 91, 132 (1990). 





[50] Ha-Ngoc Hien and J. Fabre, The velocity and shape of 2D long bubbles in inclined channels or in vertical tubes. Part I. In a stagnant liquid, Multiphase Sci. Technol. 16, 177 (2004). 





[51] B. Figueroa-Espinoza and J. Fabre, Taylor bubble moving in a flowing liquid in vertical channel: transition from symmetric to asymmetric shape, J. Fluid Mech. 679, 432 (2011). 





[52] L. Łaniewski-Wołłk and J. Rokicki, Adjoint lattice boltzmann for topology optimization on multi-gpu architecture, Comput. Math. Appl. 71, 833 (2016). 





[53] M. Schönherr, K. Kucher, M. Geier, M. Stiebler, S. Freudiger, and M. Krafczyk, Multi-thread implementations of the lattice Boltzmann method on non-uniform grids for CPUs and GPUs, Comput. Math. Appl. 61, 3730 (2011). 

