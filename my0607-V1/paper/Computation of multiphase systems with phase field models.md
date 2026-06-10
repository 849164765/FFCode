ELSEVIER 
Available online at www.sciencedirect.com 
SCIENCE@DIRECTS 
Journal of Computational Physics 190 (2003) 371–397 
JOURNAL OF COMPUTATIONAL PHYSICS 
www.elsevier.com/locate/jcp 
# Computation of multiphase systems with phase field models
V.E. Badalassi a,*, H.D. Ceniceros b , S. Banerjee a 
a Department of Chemical Engineering, University of California, Santa Barbara, CA 93106, USA 
b Department of Mathematics, University of California, Santa Barbara, CA 93106, USA 
Received 22 July 2002; received in revised form 16 April 2003; accepted 19 May 2003 
# Abstract
Phase field models offer a systematic physical approach for investigating complex multiphase systems behaviors such as near-critical interfacial phenomena, phase separation under shear, and microstructure evolution during solidification. However, because interfaces are replaced by thin transition regions (diffuse interfaces), phase field simulations require resolution of very thin layers to capture the physics of the problems studied. This demands robust numerical methods that can efficiently achieve high resolution and accuracy, especially in three dimensions. We present here an accurate and efficient numerical method to solve the coupled Cahn–Hilliard/Navier–Stokes system, known as Model H, that constitutes a phase field model for density-matched binary fluids with variable mobility and viscosity. The numerical method is a time-split scheme that combines a novel semi-implicit discretization for the convective Cahn– Hilliard equation with an innovative application of high-resolution schemes employed for direct numerical simulations of turbulence. This new semi-implicit discretization is simple but effective since it removes the stability constraint due to the nonlinearity of the Cahn–Hilliard equation at the same cost as that of an explicit scheme. It is derived from a discretization used for diffusive problems that we further enhance to efficiently solve flow problems with variable mobility and viscosity. Moreover, we solve the Navier–Stokes equations with a robust time-discretization of the projection method that guarantees better stability properties than those for Crank–Nicolson-based projection methods. For channel geometries, the method uses a spectral discretization in the streamwise and spanwise directions and a combination of spectral and high order compact finite difference discretizations in the wall normal direction. The capabilities of the method are demonstrated with several examples including phase separation with, and without, shear in two and three dimensions. The method effectively resolves interfacial layers of as few as three mesh points. The numerical examples show agreement with analytical solutions and scaling laws, where available, and the 3D simulations, in the presence of shear, reveal rich and complex structures, including strings. 
 2003 Elsevier Science B.V. All rights reserved. 
Keywords: Cahn–Hilliard equation; Navier–Stokes equations; Phase separation; Model H; Phase separation under shear flow; Interface capturing methods 
* Corresponding author. 
E-mail addresses: badalass@engineering.ucsb.edu (V.E. Badalassi), hdc@math.ucsb.edu (H.D. Ceniceros), banerjee@engineering. ucsb.edu (S. Banerjee). 
0021-9991/$ - see front matter  2003 Elsevier Science B.V. All rights reserved. doi:10.1016/S0021-9991(03)00280-8 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
# 1. Introduction
Phase field-based models replace sharp fluid/material interfaces by thin but nonzero thickness transition regions where the interfacial forces are smoothly distributed. The basic idea is to introduce an order parameter or phase field that varies continuously over thin interfacial layers and is mostly uniform in the bulk phases. Perhaps the best-known example of this type of model is the Cahn–Hilliard equation [1,2] used for modeling phase separation in a binary mixture quenched into the unstable region. The relaxation of the order parameter is driven by local minimization of the free energy subject to phase field conservation and as a result, the interface layers do not deteriorate dynamically. 
One of the applications for which phase field models are particularly well-suited is the complex process of phase separation, structure formation and evolution in flow systems, an area of technological impact in soft materials processing. The hydrodynamics can be introduced in several ways. For density-matched binary liquids, which is the case we focus on this work, this is accomplished with the coupling of the convective Cahn–Hilliard equation with a modified momentum equation that includes a phase field-dependent surface force. This is known as Model H according to the classification of Hohenberg and Halperin [3]. In the case of fluids with different densities a phase field model has been proposed by Lowengrub and Truskinovsky [4]. 
One of the salient points of the phase field description is that the order parameter has a physical meaning and different phenomena can be accounted for by a suitable modification of the free energy. Moreover, complex morphological and topological flow transitions such as coalescence and interface break-up can be captured naturally and in a mass-conservative and energy-dissipative fashion. The main drawback on the other hand is that to properly model relevant physical phenomena the interface layers have to be extremely thin. As a consequence the phase field has large gradients that must be resolved computationally. This is not an easy task. High resolution is required but the Cahn–Hilliard equation and the phase field-dependent surface force have high order derivative components. Fully implicit treatment of these terms yields expensive schemes and explicit discretizations quickly lead to numerical instability or impose impractical time-stepping constraints. 
Here, we propose an efficient and robust numerical method for the coupled Cahn–Hilliard/Navier– Stokes system. The time discretization of the method is a semi-implicit one based on an extraction of constant coefficient leading order terms (at small scales) that are time-step split. The implicit discretization of these constant coefficient terms can be inverted efficiently at optimal cost and relaxes the high order stability constraints. The time splitting allows us to decouple at each time-step the Cahn–Hilliard and the Navier–Stokes solvers. The semi-implicit discretization is combined with an original application of state-ofthe-art high-resolution schemes. We solve the flow using a robust time-discretization of the projection method that is formally second order and has a stronger high modal decay than the popular Crank–Nicolson-based projection methods. For flows confined by walls and with streamwise and spanwise periodicities, we discretize the system in space using a spectral approximation in those directions and a combination of spectral and eighth order compact finite difference approximations [5] in the wall normal direction. We demonstrate the efficacy of the method with examples of pure phase separation and binary shear flow in two and three dimensions. 
Little work has been done on the solution of the coupled Cahn–Hilliard/Navier–Stokes system [6–9] and our three-dimensional simulations for separation under shear flow are, to our knowledge, one of the first ever reported. The overall method proposed here is accurate and robust allowing interface thickness of as few as three mesh points and, as the numerical experiments show, its efficiency makes possible high-resolution 3D simulations even on modest personal computers. The numerical examples show agreement with analytical solutions and scaling laws where available and the 3D simulations in the presence of shear flow reveal rich and complex structures characterized by formation of string-like phases. 
372 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
The rest of the paper is organized as follows. Section 2 gives a brief introduction to the model coupling the phase field and the Navier–Stokes equations. Section 3 discusses our proposed numerical procedure, and in Section 4 the method is validated through numerical examples, and the results of our numerical experiments are presented and discussed. This is followed by some concluding remarks and an Appendix A. 
# 2. The governing equations
# 2.1. The phase field method
Phase field methods are a particular class of diffuse-interface models that have been used successfully in the study of critical phenomena but have not been used much for fluid interfaces. In a phase field method, it is assumed that the state of the system at any given time can be described by an order parameter $\phi$ which is a function of the position vector. For example, in the case of an isothermal binary fluid $\phi$ is the relative concentration of the two components. A free energy can be defined for times when the system is not in equilibrium [10], and this free energy can be written as a functional of $\phi \colon$ 
$$
F [ \phi ] = \int_ {\Omega} \left\{f (\phi (\mathbf {x})) + \frac {1}{2} k | \nabla \phi (\mathbf {x}) | ^ {2} \right\} d \mathbf {x}, \tag {1}
$$
where $\Omega$ is the region of space occupied by the system. The term $\big ( 1 / 2 \big ) k \big | \nabla \phi ( \mathbf { x } ) \big | ^ { 2 }$ accounts for the surface energy, with $k$ a positive constant, and $f ( \phi ( \mathbf { x } ) )$ is the bulk energy density which has two minima corresponding to the two stable phases of the fluid. 
The chemical potential $\mu$ is defined as 
$$
\mu (\phi) = \frac {\delta F [ \phi ]}{\delta \phi (\mathbf {x})} = f ^ {\prime} (\phi (\mathbf {x})) - k \nabla^ {2} \phi (\mathbf {x}). \tag {2}
$$
The equilibrium interface profile can be found by minimizing the functional $F [ \phi ]$ with respect to variations of the function $\phi ,$ i.e., solving $\mu ( \phi ) = 0 .$ . Cahn and Hilliard [1,2] generalized the problem to timedependent situations by approximating interfacial diffusion fluxes as being proportional to chemical potential gradients, enforcing conservation of the field. The convective Cahn–Hilliard equation can be written as 
$$
\frac {\partial \phi}{\partial t} + \mathbf {u} \cdot \nabla \phi = \nabla \cdot (M (\phi) \nabla \mu), \tag {3}
$$
where u is the velocity field and $M ( \phi ) > 0$ is the mobility or Onsager coefficient. Eq. (3) models the creation, evolution, and dissolution of diffusively controlled phase field interfaces [11] (for a review of the Cahn– Hilliard model see for example [12]). At the wall, we adopt the following no-flux boundary conditions: 
$$
\mathbf {n} \cdot \nabla \phi = 0 \quad \text { and } \quad \mathbf {n} \cdot M \nabla \mu = 0, \tag {4}
$$
where n is the unit vector normal to the domain boundary. 
# 2.2. The equations of fluid motion
This work focuses on density-matched binary mixtures with variable viscosity and mobility. The fluid dynamics are described by the Navier–Stokes equations with a phase field-dependent surface force [13] 
$$
\rho \left(\frac {\partial \mathbf {u}}{\partial t} + \mathbf {u} \cdot \nabla \mathbf {u}\right) = - \nabla p + \nabla \cdot \eta (\nabla \mathbf {u} + \nabla \mathbf {u} ^ {\mathrm{T}}) + \mu \nabla \phi , \tag {5}
$$
373 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
$$
\nabla \cdot \mathbf {u} = 0, \tag {6}
$$
where u is the velocity field, $p$ is a scalar related to the pressure that enforces the incompressibility constraint (6), and $\eta$ is the viscosity. The superscript T stands for the transpose operator. At a wall the Dirichlet boundary condition is imposed for the velocity field, i.e., $\mathbf { u } = \mathbf { u } _ { 0 }$ at a fixed domain boundary. 
The coupled Cahn–Hilliard/Navier–Stokes system (3)–(6) is referred to as ‘‘Model $\mathrm { H } '$ according to the nomenclature of Hohenberg and Halperin [3]. 
# 2.3. Interface properties
For the binary fluid we use the following double well potential: 
$$
f (\phi) = \frac {\alpha}{4} \left(\phi - \sqrt {\frac {\beta}{\alpha}}\right) ^ {2} \left(\phi + \sqrt {\frac {\beta}{\alpha}}\right) ^ {2}, \tag {7}
$$
where a and $\beta$ are two positive constants. The equilibrium profile is given by the solutions of the equation 
$$
\mu (\phi) = \frac {\delta F [ \phi ]}{\delta \phi} = \alpha \phi^ {3} - \beta \phi - k \nabla^ {2} \phi = 0. \tag {8}
$$
This leads to two stable uniform solutions $\phi _ { \pm } = \pm \sqrt { \beta / \alpha }$ representing the coexisting bulk phases, and a one-dimensional (say along the z-direction) nonuniform solution 
$$
\phi_ {0} (z) = \phi_ {+} \tanh \left(\frac {z}{\sqrt {2} \xi}\right) \tag {9}
$$
that satisfies the boundary conditions $\phi _ { 0 } ( z  \pm \infty ) = \pm \phi ( \sec { [ 6 , 1 4 ] } )$ . This solution was first found by van der Waals [15] to describe the equilibrium profile for a plane interface normal to the z-direction, of thickness proportional to $\zeta = \sqrt { k / \beta } ,$ that separates the two bulk phases. 
We define the interface thickness to be the distance fromp $0 . 9 \phi .$ - to $0 . 9 \phi _ { + }$ so that the equilibrium interface thickness is $2 \sqrt { 2 } \xi \operatorname { t a n h } ^ { - 1 } ( 0 . 9 ) = 4 . 1 6 4 \xi$ . This width contains 98.5% of the surface tension stress [7]. 
In equilibrium the surface tension r of an interface is equal to the integral of the free energy density along the interface. For a plane interface r is given by [14] 
$$
\sigma = k \int_ {- \infty} ^ {+ \infty} \left(\frac {\mathrm{d} \phi_ {0}}{\mathrm{d} z}\right) ^ {2} \mathrm{d} z = \frac {\sqrt {2}}{3} \frac {k ^ {1 / 2} \beta^ {3 / 2}}{\alpha}. \tag {10}
$$
It is evident from (9) and (10) that we can control the surface tension and interface width through the parameters k, a, and $\beta .$ 
# 2.4. Nondimensionalization
We nondimensionalize the governing equations with the variables 
$$
\mathbf {u} ^ {\prime} = \frac {\mathbf {u}}{U _ {\mathrm{c}}}, \quad t ^ {\prime} = \frac {t}{T _ {\mathrm{c}}}, \quad \mathbf {x} ^ {\prime} = \frac {\mathbf {x}}{L _ {\mathrm{c}}}, \quad p ^ {\prime} = \frac {p L _ {\mathrm{c}}}{\eta_ {\mathrm{c}} U _ {\mathrm{c}}}. \tag {11}
$$
Following Chella and Vinals [6] we choose as characteristic length ~ $L _ { \mathrm { c } }$ the mean-field thickness $\xi$ of the interface, i.e., $L _ { \mathrm { c } } = \xi$ . The characteristic velocity $U _ { \mathrm { c } }$ depends on the problem; for example, it could be the imposed velocity in shear flow. The characteristic time $T _ { \mathrm { c } }$ is the time required for the fluid to be convected a distance of the order of the interface thickness (in the absence of capillarity), $T _ { \mathrm { c } } = \xi / U _ { \mathrm { c } }$ . The order parameter $\phi$ is scaled with its mean-field equilibrium value $\phi _ { + } = \sqrt { \beta / \alpha } .$ Dropping the primes, Eqs. (3)–(6) become 
374 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 

$$
\frac {\partial \phi}{\partial t} + \mathbf {u} \cdot \nabla \phi = \frac {1}{P e} \nabla \cdot \lambda \nabla \mu , \tag {12}
$$
$$
R e \left(\frac {\partial \mathbf {u}}{\partial t} + \mathbf {u} \cdot \nabla \mathbf {u}\right) = - \nabla p + \nabla \cdot \theta (\nabla \mathbf {u} + \nabla \mathbf {u} ^ {\mathrm{T}}) + \frac {1}{C a} \mu \nabla \phi , \tag {13}
$$
$$
\nabla \cdot \mathbf {u} = 0, \tag {14}
$$
where $\theta = \eta / \eta _ { \mathrm { c } }$ and $\lambda = M / M _ { \mathrm { c } }$ are the normalized viscosity and mobility, respectively, and $\mu = \phi ^ { 3 } - \phi -$ $\nabla ^ { 2 } \phi$ is the dimensionless chemical potential. The dimensionless groups used above are the Reynolds number, the Peclet number, and the capillary number given by 
$$
R e = \frac {\rho U _ {\mathrm{c}} \xi}{\eta}, \quad P e = \frac {U _ {\mathrm{c}} \xi}{M _ {\mathrm{c}} \beta}, \quad C a = \frac {\alpha \eta U _ {\mathrm{c}}}{\beta^ {2} \xi} = \frac {2 \eta U _ {\mathrm{c}}}{3 \sigma}, \tag {15}
$$
respectively. Physically, the Peclet number Pe is the ratio between the diffusive time scale $\xi ^ { 2 } / ( M _ { \mathrm { c } } \beta )$ and the convective time scale $\xi / U _ { \mathrm { c } }$ . The Reynolds number $R e$ is the ratio between inertial and viscous forces and the capillary number Ca provides a measure of the relative magnitude of viscous and capillary (or interfacial tension) forces at the interface. Note that with this nondimensionalization the length of the fluid domain is interpreted in units of interface thickness $\xi .$ 
We consider the viscosity g as a linear function of the order parameter $\phi .$ . That is, if $\eta _ { - } \leqslant \eta \leqslant \eta _ { + }$ and $\eta _ { \mathrm { c } } = \eta$ - we get 
$$
\theta = \frac {\theta_ {\max} - 1}{2} \phi + \frac {\theta_ {\max} + 1}{2}, \tag {16}
$$
where $\theta _ { \mathrm { m a x } } = \eta _ { + } / \eta _ { - }$ is the viscosity ratio. In this way g automatically changes across the interface with a profile similar to the tanh function. 
For the mobility M we follow [16] and we consider a profile as $M = M _ { \mathrm { c } } ( 1 - \gamma \phi ^ { 2 } )$ so that we have 
$$
\lambda = (1 - \gamma \phi^ {2}). \tag {17}
$$
where $0 \leqslant \gamma \leqslant 1 . \mathrm { I f } \ \gamma \stackrel { } { \longrightarrow } 0$ we have phase separation dynamics controlled by bulk diffusion, if $\gamma  1$ we have dynamics controlled by interface diffusion. 
# 3. The numerical method
# 3.1. Temporal discretization
We propose a semi-implicit time discretization combined with a time-split strategy. This discretization effectively decouples Cahn–Hilliard and Navier–Stokes solvers and yields an efficient and robust modular scheme. 
The outline of the method is as follows. Given $\phi ^ { n }$ and $ { \mathbf { u } } ^ { n }$ the objective is to solve for $\phi ^ { n + 1 }$ and $ { \mathbf { u } } ^ { n + 1 }$ with the steps: 
(1) Solve the Cahn–Hilliard equation with a second order semi-implicit method and spectral spatial discretization to obtain $\phi ^ { n + 1 }$ . 
375 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
(2) Using $\phi ^ { n + 1 }$ compute the surface force and solve the phase field modified Navier–Stokes equations with a second order SBDF (semi-backward difference formula)-based projection method to obtain $ { \mathbf { u } } ^ { n + 1 }$ . The spatial discretization is spectral in the streamwise and spanwise directions and eighth order compact finite difference in the wall normal direction [5,17,18]. 
Our semi-implicit strategy uses a simple idea that works quite well for diffusion-dominated equations, for example, the variable (even nonlinear) coefficient diffusion equation $u _ { t } = \nabla \cdot ( \chi \nabla u ) , \chi > 0 [ 1 9 , 2 0 ]$ . We rewrite the latter as 
$$
\frac {\partial u}{\partial t} = a \nabla^ {2} u + f (u), \tag {18}
$$
where $f ( u ) = \nabla \cdot ( \chi \nabla u ) - a \nabla ^ { 2 } u$ and a is constant in space (but could be time-dependent). By treating the first term on the right-hand side of (18) implicitly and $f ( u )$ explicitly we can obtain semi-implicit discretizations that can be easily solved. With energy estimates one can show that a first order Euler discretization is unconditionally stable if $a \geqslant ( 1 / 2 )$ max $\chi \ [ 1 9 ]$ . Since the truncation error is dissipative and proportional to $^ { a , }$ we consider $a = ( 1 / 2 )$ max v as an optimal value. Discretizations of this type are of common use in spectral methods [20] as the constant coefficient implicit terms becomes diagonal in Fourier space and thus can be inverted efficiently. However, as noted in [21], these discretizations are less successful for dispersiondominated problems. 
We can apply this idea to deal with variable mobility. However, the application of the same idea to the treatment of the nonlinear term due to the chemical potential is not straightforward. To achieve this, we note that $\nabla ^ { 2 } f ^ { \prime } ( \phi ) = \nabla \cdot ( f ^ { \prime \prime } \nabla \phi )$ where $f ^ { \prime } ( \phi ) = \phi ^ { 3 } - \phi$ and $f ^ { \prime \prime } ( \phi ) = 3 \phi ^ { 2 } - 1$ . Letting $\tau = ( 1 / 2 )$ max $\left( f ^ { \prime } ( \phi ) \right) =$ $( 1 / 2 ) f ^ { \prime \prime } ( \pm 1 ) = 1$ and defining $\lambda _ { \operatorname* { m a x } } = \operatorname* { m a x } \lambda$ as the maximum of the normalized mobility (the mobility ratio if $M _ { \mathrm { c } } = M _ { - } )$ we rewrite (12) as 
$$
\frac {\partial \phi}{\partial t} = \frac {1}{P e} \frac {\lambda_ {\max}}{2} \left[ \tau \nabla^ {2} \phi - \nabla^ {4} \phi \right] + \frac {1}{P e} [ A (\phi) + B (\phi) ] - \mathbf {u} \cdot \nabla \phi , \tag {19}
$$
where $A ( \phi ) = \nabla \cdot \lambda \nabla f ^ { \prime } ( \phi ) - ( \lambda _ { \operatorname* { m a x } } / 2 ) \tau \nabla ^ { 2 } \phi$ and $B ( \phi ) = ( \lambda _ { \operatorname* { m a x } } / 2 ) \nabla ^ { 4 } \phi - \nabla \cdot \lambda \nabla ( \nabla ^ { 2 } \phi )$ : By treating the first term on the right-hand side of (19) implicitly and $A ( \phi ) , B ( \phi )$ and the convective term u  $\nabla \phi$ explicitly we can obtain semi-implicit discretizations that can be solved efficiently at minimal cost. When looking for a second order semi-implicit multi-step method it is fundamental to note that because of the very high frequency content in the Cahn–Hilliard solutions we need a method with high modal damping. The use of weakly damping schemes such as the popular combination of Crank–Nicolson with second or higher order convective terms discretizations is not appropriate [22] since it can lead to extra iterations on the finest grid when using multigrid methods with finite difference spatial discretizations, and to aliasing, when using spectral collocation for spatial discretization as it is in our case. Among the second order multi-step methods the extrapolated Gear (SBDF) scheme has the strongest high modal decay [22]. We experimented numerically with the Crank–Nicolson discretization applied to the modified Cahn–Hilliard equation (19), without convection, and found that unless Crank–Nicolson is used in its dissipative regime $( \Delta t < C h ^ { 2 } )$ it would be unstable. The high modal damping is apparently required to stabilize the high frequency content of the explicitly treated difference between the variable coefficient term and the constant one. The SBDF provides, without the stringent quadratic time-step constraint, the required damping. Applied to (19) this scheme becomes 
$$
\begin{array}{l} \frac {\frac {3}{2} \phi^ {n + 1} - 2 \phi^ {n} + \frac {1}{2} \phi^ {n - 1}}{\Delta t} = \frac {1}{P e} \frac {\lambda_ {\max}}{2} \left[ \tau \nabla^ {2} \phi^ {n + 1} - \nabla^ {4} \phi^ {n + 1} \right] + 2 \left\{\frac {1}{P e} [ A (\phi^ {n}) + B (\phi^ {n}) ] - \mathbf {u} ^ {n} \cdot \nabla \phi^ {n} \right\} \\ - \left\{\frac {1}{P e} \left[ A (\phi^ {n - 1}) + B (\phi^ {n - 1}) \right] - \mathbf {u} ^ {n - 1} \cdot \nabla \phi^ {n - 1} \right\}. \tag {20} \\ \end{array}
$$
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
In the absence of convection, this new discretization appears in our numerical experiments to be unconditionally stable. Eyre [23] considers a discretization of this type for the one-dimensional Cahn–Hilliard equation with constant mobility, Smereka [24] uses it in the context of interface motion by surface diffusion while Zhu et al. [25] use it for the mobility term but not for the nonlinear one, resulting in a conditionally stable method. 
For the Navier–Stokes equations (13) and (14) we use the Gear scheme combined with the above semiimplicit discretization applied to the variable viscosity term. This discretization will provide the necessary damping for the high mode components due to the near discontinuities in the derivatives of the velocity and the presence of the almost singular surface-tension source term. The discretized Navier–Stokes equations are 
$$
\begin{array}{l} \frac {\frac {3}{2} \mathbf {u} ^ {n + 1} - 2 \mathbf {u} ^ {n} + \frac {1}{2} \mathbf {u} ^ {n - 1}}{\Delta t} = - \frac {\nabla p ^ {n + 1}}{R e} + \frac {\theta_ {\max}}{2 R e} \nabla^ {2} \mathbf {u} ^ {n + 1} + \frac {1}{R e C a} \mu (\phi^ {n + 1}) \nabla \phi^ {n + 1} + 2 \left[ \frac {C (\mathbf {u} ^ {n})}{R e} - \mathbf {u} ^ {n} \cdot \nabla \mathbf {u} ^ {n} \right] \\ - \left[ \frac {C (\mathbf {u} ^ {n - 1})}{R e} - \mathbf {u} ^ {n - 1} \cdot \nabla \mathbf {u} ^ {n - 1} \right], \tag {21} \\ \end{array}
$$
where $C ( \mathbf { u } ^ { m } ) = \nabla \cdot \theta ^ { m + 1 } ( \nabla \mathbf { u } ^ { m } + ( \nabla \mathbf { u } ^ { m } ) ^ { \mathrm { T } } ) - ( \theta _ { \operatorname* { m a x } } / 2 ) \nabla ^ { 2 } \mathbf { u } ^ { m }$ with $m = n , n - 1$ . Karniadakis et al. [26] employ this scheme in the context of single phase flow, i.e., without the source term and the scheme was called ‘‘stiffly’’ stable. We use the same splitting with the addition of the source, i.e., surface tension, term. The method can be summarized as follows: 
Step 1: 
$$
\frac {\mathbf {u} ^ {*} - 2 \mathbf {u} ^ {n} + \frac {1}{2} \mathbf {u} ^ {n - 1}}{\Delta t} = 2 \left[ \frac {C (\mathbf {u} ^ {n})}{R e} - \mathbf {u} ^ {n} \cdot \nabla \mathbf {u} ^ {n} \right] - \left[ \frac {C (\mathbf {u} ^ {n - 1})}{R e} - \mathbf {u} ^ {n - 1} \cdot \nabla \mathbf {u} ^ {n - 1} \right] + \frac {1}{R e C a} \mu (\phi^ {n + 1}) \nabla \phi^ {n + 1}. \tag {22}
$$
Step 2: 
$$
\frac {\mathbf {u} ^ {* *} - \mathbf {u} ^ {*}}{\Delta t} = - \frac {\nabla p ^ {n + 1}}{R e}. \tag {23}
$$
Step 3 (Helmholtz equation): 
$$
\frac {\frac {3}{2} \mathbf {u} ^ {n + 1} - \mathbf {u} ^ {* *}}{\Delta t} = \frac {\theta_ {\max}}{2 R e} \nabla^ {2} \mathbf {u} ^ {n + 1} \tag {24}
$$
with Dirichlet boundary conditions 
$$
\mathbf {u} ^ {n + 1} = \mathbf {u} _ {0}. \tag {25}
$$
We need to introduce two further assumptions for the intermediate velocity fields $\mathbf { u } ^ { * } , \mathbf { u } ^ { * * }$ . First the incompressibility constraint 
$$
\nabla \cdot \mathbf {u} ^ {* *} = 0 \tag {26}
$$
and second that the same field u also satisfies the prescribed Dirichlet condition in the direction normal to the boundary 
$$
\mathbf {u} ^ {* *} \cdot \mathbf {n} = \mathbf {u} _ {0} \cdot \mathbf {n}. \tag {27}
$$
Incorporating these assumptions into Eq. (23) we finally derive a separately solvable equation for the pressure (Poisson equation) 
$$
\nabla^ {2} p ^ {n + 1} = \frac {R e}{\Delta t} \nabla \cdot \mathbf {u} ^ {*}. \tag {28}
$$
377 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
Karniadakis et al. [26] derive the Neumann boundary conditions that allows second order accuracy in the velocity and pressure in the context of single phase flow and constant viscosity. We follow the same procedure and we evaluate the normal component of (21) at the boundary and let the term $\nabla ^ { 2 } \mathbf { u } = - \nabla \times ( \nabla \times \mathbf { u } )$ (due to the incompressibility constraint (6)) to yield 
$$
\begin{array}{l} \left. \frac {\partial p ^ {n + 1}}{\partial n} \right| _ {\Gamma} = \mathbf {n} \cdot \left[ 2 (R e \mathbf {u} ^ {n} \cdot \nabla \mathbf {u} ^ {n} + \nabla \theta^ {n + 1} \nabla \mathbf {u} ^ {n} - \theta^ {n + 1} \nabla \times (\nabla \times \mathbf {u} ^ {n}) + \nabla \cdot \theta^ {n + 1} (\nabla \mathbf {u} ^ {n}) ^ {\mathrm{T}}) \right. \\ \left. - \left(R e \mathbf {u} ^ {n - 1} \cdot \nabla \mathbf {u} ^ {n - 1} + \nabla \theta^ {n} \nabla \mathbf {u} ^ {n - 1} - \theta^ {n} \nabla \times (\nabla \times \mathbf {u} ^ {n - 1}) + \nabla \cdot \theta^ {n} (\nabla \mathbf {u} ^ {n - 1}) ^ {\mathrm{T}}\right) \right]. \tag {29} \\ \end{array}
$$
Note that we calculate the term $\nabla ^ { 2 } \mathbf { u } ^ { n + 1 }$ with an extrapolation from the time levels n and $n - 1$ . 
# 3.2. Stability
A rigorous stability analysis for the overall scheme is quite difficult. Nevertheless one can obtain valuable information about the stability and robustness of the scheme through numerical tests (see Section 4). In particular, through numerical experiments we find that the semi-implicit discretization (20) for the Cahn– Hilliard equation appears to be unconditionally stable when $\mathbf { u } \equiv \mathbf { 0 } _ { : }$ , regardless of the interface thickness. Moreover, the unconditional stability seems to hold for $a \geqslant f ^ { \prime \prime } ( \pm 1 ) / 2 = 1$ just as for the corresponding discretization of the variable diffusion equation $\boldsymbol { u } _ { t } = \nabla \cdot \left( \chi \nabla \boldsymbol { u } \right)$ . Thus, for a given nonzero u, the scheme for the convective Cahn–Hilliard equation has only a CFL stability condition 
$$
\Delta t _ {\mathrm{cfl}} \leqslant \left(\frac {| u | _ {\max}}{\Delta x} + \frac {| v | _ {\max}}{\Delta y} + \frac {| w | _ {\max}}{\Delta z}\right) ^ {- 1}, \tag {30}
$$
where $( u , v , w )$ are the components of the velocity field. 
When coupled with the time discretization of the modified Navier–Stokes equations (22)–(24), in addition to the natural CFL condition, we have to consider time step restrictions due to surface tension and viscosity. For the surface tension we observe a mild stability constraint of the form 
$$
\Delta t _ {s} \leqslant C _ {1} \sqrt {\text { Re }   C a} (\min \{\Delta x, \Delta y, \Delta z \}) ^ {3 / 2}, \tag {31}
$$
where $C _ { 1 }$ is a constant. $C _ { 1 } = 1 0$ , works well for all our numerical examples. Note that spatial mesh sizes are nondimensional so that min $( \Delta x , \Delta y , \Delta z ) = \mathbf { O } ( 1 )$ . The same type of condition is found for capturing $\mathrm { ( ^ { 6 6 } c o l o r ^ { 9 } ) }$ 号 methods (with the appropriate nondimensionalization) such as the level set method [27] and the continuum surface force method (CSF) [28] that rely both on less stiff evolution equations for the ‘‘color’’ function. 
We now derive the stability constraint associated with the variable viscosity term. Using the incompressibility condition, the Navier–Stokes equation (13) in indicial notation (repeated index summation implied) becomes 
$$
\left(\frac {\partial u _ {i}}{\partial t} + u _ {k} \frac {\partial u _ {i}}{\partial x _ {k}}\right) = - \frac {1}{R e} \frac {\partial p}{\partial x _ {i}} + \frac {1}{R e} \left\{\frac {\partial}{\partial x _ {k}} \left(\theta \frac {\partial u _ {i}}{\partial x _ {k}}\right) + \frac {\partial \theta}{\partial x _ {k}} \frac {\partial u _ {k}}{\partial x _ {i}} \right\} + \frac {1}{C a R e} \mu \frac {\partial \phi}{\partial x _ {i}}. \tag {32}
$$
The semi-implicit discretization removes the severe stability constraint due to the term $( \hat { \ v O } / \hat { \ v O } x _ { k } ) ( \theta ( \hat { \ v O } u _ { i } / \hat { \ v O } x _ { k } ) )$ but has limited effect on the term $\widehat { \sf { O } } \theta / \widehat { \sf { O } } x _ { k } \widehat { \sf { O } } u _ { k } / \widehat { \sf { O } } x _ { i }$ . This term gives rise to a CFL-like stability constraint that can be determined by estimating max $| \nabla \theta |$ . In the limit of gently curved interfaces, and when the motion of the interface is slow compared with the local relaxation times of $\phi ,$ we can approximate $\phi$ by the one-dimensional stationary solution $\phi _ { 0 }$ in (9) along the direction normal to the interface, i.e., $\nabla \phi \simeq \nabla \phi _ { 0 }$ . From (9) and (16) we have that $\nabla \theta \propto ( \theta _ { \mathrm { m a x } } - 1 ) \sec h ^ { 2 } x$ , then max $| \nabla \theta | \propto$ $( \theta _ { \mathrm { m a x } } - 1 )$ . Thus, the variable viscosity time-step constraint has the form 
378 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
$$
\Delta t _ {\mathrm{vr}} \leqslant C _ {2} \frac {R e}{\theta_ {\max} - 1} (\min \{\Delta x, \Delta y, \Delta z \}), \tag {33}
$$
where $C _ { 2 }$ is a constant. For $\theta _ { \mathrm { m a x } } = 1$ the discretization is unconditionally stable since it reduces to an implicitly treated constant viscosity case. For $\theta _ { \mathrm { { m a x } } } > 1$ we could use successfully $C _ { 2 } = 1 0$ for all our simulations. Note that if we were treating the viscous term purely explicitly we would have the more restrictive constraint $\Delta t \leqslant ( R e / \theta _ { \mathrm { m a x } } ) [ ( \Delta x ) ^ { - 2 } + ( \bar { \Delta } y ) ^ { - 2 } + ( \Delta z ) ^ { - 2 } ] ^ { - 1 }$ . 
We can now express our adaptive time-stepping strategy as 
$$
\Delta t ^ {n + 1} = \min (\Delta t _ {\mathrm{cfl}}, \Delta t _ {\mathrm{s}}, \Delta t _ {\mathrm{vr}}). \tag {34}
$$
The discretization (20) effectively removes the high order stability constraints associated with the Cahn– Hilliard equation and makes the phase field-based method computationally competitive and robust. To relax more the viscous stability constraint in the case of very small Re one can use a predictor–corrector iteration strategy. Increasing the constant leading order term $\theta _ { \mathrm { m a x } } / 2$ in (21) also relaxes the constraint by allowing a larger constant $C _ { 2 }$ , albeit at the cost of increasing the truncation error. For example if $\theta _ { \mathrm { m a x } }$ is used instead of $( 1 / 2 ) \theta _ { \mathrm { m a x } }$ we find that one can use $C _ { 2 } = 1 8 0$ giving a significant saving in time stepping. 
# 3.3. Spatial discretization
We employ high-resolution spatial discretizations to be able to accurately resolve thin interfaces. The Cahn–Hilliard equation is discretized in space (pseudo) spectrally (via FFT for periodic boundary conditions or Cosine transform for the no-flux conditions). For the Navier–Stokes equations we use spectral derivatives in the streamwise and spanwise (periodic) directions and an eighth order finite difference compact scheme [5] for the wall normal derivatives of the velocity and pressure. Note that compact finite difference approximations are used only for the wall normal derivatives of the velocity in (22)–(24) and for the first order wall normal derivative of / in (20). We compute the other derivatives spectrally in the x- and y-directions with the fast Fourier transform (FFT). The details of the spatial discretization are given in Appendix A. 
# 4. Numerical experiments and validation
We present three types of numerical experiments to validate the proposed method and test its capabilities. The experiments are simulations of drop deformation, pure phase separation (spinodal decomposition) and phase separation under shear flow. A resolution study is also performed to check the accuracy and the stability of the method. This is briefly described next. 
# 4.1. Drop deformation in a shear flow
We consider an initially 2D spherical drop in a shear flow. This is a classical problem that was solved analytically for sharp interfaces and small deformations in the creeping flow approximation for unbounded domain by Taylor [29] and in the presence of two walls by Shapira and Haber [30]. The drop will assume the shape of an ellipsoid with a deformation that depends on the capillary number and the viscosity ratio. Taylor [29] found that for equal viscosity blends at steady state, i.e., when deformation due to the externally imposed shear flow and interfacial relaxation balance one another, the deformation parameter $D = ( l - s ) / ( l + s )$ is related to the capillary number Ca as 
$$
D = \frac {3 5}{3 2} C a, \tag {35}
$$
379 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
where l; s denote, respectively, the longest and shortest axes of the ellipsoid in the shear gradient plane. This relation is valid in the limit of vanishing deformations and holds in good approximation for $D < 0 . 3$ . We use this problem to demonstrate the convergence and accuracy of the numerical results under grid refinement. At the same time, we validate the calculation of surface tension and viscosity ratio. 
As initial condition we start with a 2D circular drop in the center of the domain with a ‘‘tanh’’ profile of the interface and we solve the Cahn–Hilliard equation without convection to reach a steady state that leads to a completely saturated mixture. Then we impose a shear flow with the top and bottom lid moving in opposite directions and with the dimensionless velocity equal to plus or minus one, respectively. We consider three capillary numbers $C a = 0 . 6 , 0 . 9$ , and 1.2. The fluids have the same viscosity and $R e = 0 . 1 $ , $P e = 1 0$ . We employ two resolutions $1 2 8 \times 1 2 8$ and $2 5 6 \times 2 5 6 ,$ , and domain sizes of $L = 1 7 8$ and $L = 3 5 5$ , respectively. This combination of parameters determines an interface thickness of three mesh points. Recall 
/ae04e2bfed3d5546c639e6b9c32b30b13aa3c186df6ee0f4e01521ce00164c97.jpg
/b75412df7bcd2572349765523bd6d2db77f24d07a36f3265adbfc418b54c90fb.jpg
/c883c1734bdfa659f00f1b7b447fe15c69ef243ded6f5ebf3e1c47cd5e553134.jpg
/a6403428fa97b8936611b1e1a33dd2c0a0506526133bb316bb410ecf1defed6e.jpg

Fig. 1. 2D deformation of an initially spherical drop described by /, $P e = 1 0 ,$ $R e = 0 . 1$ for: (a) $C a = 0 . 6 ,$ (b) $C a = 0 . 9 ,$ and (c) $C a = 1 . 2 .$ First column $N = 1 2 8$ and $L = 1 7 8 ,$ , second column $N = 2 5 6$ and $L = 3 5 5 .$ .

380 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
/40d1af6dd1eee931fd3d45f00cd9ed85fec7d7c0cb3c7c4d15334fb75f6a6f55.jpg
/c43ca0a6d74f4f8e579865799f193be45bd239c8a85cb5736445e90be0299698.jpg

Fig. 1. (continued)

that, based on our nondimensionalization, the length of the domain L is to be interpreted in units of the interface thickness n. In Fig. 1 we plot the final equilibrium stage. The convergence of the results under grid refinement is evident. The drop shape is ellipsoidal with the major axis converging to an angle of $4 5 ^ { \circ }$ as Ca decreases, just as predicted by the analytic results [29,30]. Moreover, the contours )0.9 and 0.9 describing the interface are well behaved since they stay parallel throughout the computation. Increasing the capillary number results in increased deformation and the angle diminishes in the direction of the major axis to the undisturbed (horizontal) streamlines. This result matches the numerical result of Rallison [31] for deformations where the analytical solution is not available. 
Now we perform 3D simulations with a grid size $1 2 8 \times 1 2 8 \times 1 2 8$ and $L = 1 7 8$ which correspond to a three mesh-point thick interface. We choose the droplet radius to be 35 grid points which is large enough to avoid effects due to the finite interfacial width and the presence of the walls [30]. We see from Fig. 2 that the deformation parameters obtained from simulations with different capillary numbers Ca and $R e = 0 . 0 1$ , $P e = 1 0 0$ (error bars) correspond well with the theoretical predictions of (35). The error bars in Fig. 2 result from errors due to the use of diffuse interface (i.e., errors in estimating l and s) and the use of a finite Re instead of a $R e = 0$ . These numerical results are analogous to the ones reported in [32]. 
To test the accuracy of the time discretization we perform a sequence of simulations with $5 1 2 \times 5 1 2$ mesh points keeping the spatial resolution fixed and halving the time step. As a parameter we use interface thicknesses of three, four and five mesh points (respectively L ¼ 711, L ¼ 533, L ¼ 426). Again $C a = 1 . 5 .$ , $P e = 1 0$ and $R e = 1$ and the mean initial drop radius is 128 mesh points. Denoting by $\mathbf { V } _ { i , j } ( \Delta t ) = ( \phi _ { i , j } , \mathbf { u } _ { i , j } )$ the approximation obtained using a step-size $\Delta t ,$ and defining the error ratios 
$$
R (\Delta t) = \frac {\sum_ {i , j = 1} ^ {N _ {x} , N _ {y}} | \mathbf {V} _ {i , j} (\Delta t) - \mathbf {V} _ {i , j} (\Delta t / 2) |}{\sum_ {i , j = 1} ^ {N _ {x} , N _ {y}} | \mathbf {V} _ {i , j} (\Delta t / 2) - \mathbf {V} _ {i , j} (\Delta t / 4) |}, \tag {36}
$$
we calculate the order of convergence in time as 
$$
\mathrm{O} (V) = \frac {\log R (\Delta t)}{\log 2}. \tag {37}
$$
In Table 1 we show the results for $\phi$ and for the w component of $\mathbf { u } _ { i , j } = ( u _ { i , j } , w _ { i , j } )$ . 
381 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
/33972b62cd3ae23633000babbe2ab6c8cc036ba9487b70a8c66db751d94d5e3c.jpg

Fig. 2. Deformation of a droplet subject to shear flow. The continuous line represent the theoretical prediction of Taylor [29] while the error bars represent the steady-state deformation parameter D obtained from 3D simulations as a function of the capillary number Ca.


Table 1 Order of convergence in time

<table><tr><td>Interfacial thickness</td><td>O(φ)</td><td>O(w)</td></tr><tr><td>3</td><td>1.01</td><td>0.98</td></tr><tr><td>4</td><td>1.20</td><td>1.18</td></tr><tr><td>5</td><td>1.42</td><td>1.39</td></tr></table>
To test the accuracy of the space discretization we compute a sequence of 2D simulations with $C a = 1 . 5 ,$ $P e = 1 0$ and $R e = 1$ to time 2.0 using a fixed time step that satisfies the stability requirement on three different grids, $1 2 8 \times 1 2 8 , 2 5 6 \times 2 5 6$ and $5 1 2 \times 5 1 2$ with L ¼ 133, L ¼ 266 and $L = 5 3 3$ , respectively (four mesh-points thick interface). Further, we set the initial drop mean radius equal to 1/4 the domain length. We proceed as for the time accuracy check (Eqs. (36) and (37)) and by comparing the difference in the numerical solutions for adjacent resolutions we estimate the maximum error point-wise. We then use these estimates of the error to compute a numerical convergence rate. For points close to the periodic boundaries and away from the interface we found a convergence rate of 2.7 for the velocity field and 2.8 for $\phi .$ In points close to the walls and to the interface the convergence rate deteriorates with, respectively, 1.9 and 0.91 for the velocity and 2.1 and 0.95 for $\phi .$ These results compare favorably with the level set [27] and volume of fluid [33] methods. Even though we cannot preserve spectral accuracy due to the presence of the interface, the high accuracy discretization is important as interface layers of only a few mesh points need to be resolved and numerical diffusion has to be limited to avoid unphysical coalescence of interfaces. 
Finally, we examine drop deformation in the case of variable viscosity. Shown in Fig. 3 are the results for $P e = 1 0 , R e = 0 . 1 , C a = 0 . 8$ . Three viscosity ratios are considered: $\theta _ { \mathrm { m a x } } = 2 , 5$ , and 10. We plot the contour of $\phi = 0$ only. The observed deformation increases in accordance with the predictions in [29,30] without any appreciable change in the orientation. 
382 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
/651e1544e2a5b5329355cf828bfea6a68ad28dc7a7e9c3d7b5e2d2c9cdb2891c.jpg

Fig. 3. 2D deformation of an initially spherical drop described by /, Pe ¼ 10, Re ¼ 0:1, Ca ¼ 0:8, hmax ¼ 2, 5 and 10. $N = 2 5 6$ and $L = 3 5 5$ .


# 4.2. 2D and 3D phase separation
We begin the numerical experiments with an example of pure spinodal phase separation of a binary mixture. An initially homogeneous disordered phase separates into ordered structures when quenched into a metastable region. The Cahn–Hilliard equation (without convection) models this process. For pure phase separation it is convenient to nondimensionalize (3), with $\mathbf u = 0$ , using variables (11) with $L _ { \mathrm { c } }$ as the domain size and $T _ { \mathrm { c } } = M _ { \mathrm { c } } \beta$ . Dropping the primes, Eq. (3) becomes 
$$
\frac {\partial \phi}{\partial t} = \nabla \cdot (1 - \gamma \phi^ {2}) \nabla (f ^ {\prime} (\phi) - C ^ {2} \nabla^ {2} \phi), \tag {38}
$$
where $C = \xi / L _ { \mathrm { c } }$ is the Cahn number and $f ^ { \prime } ( \phi ) = \phi ^ { 3 } - \phi$ . The Cahn number represents the ratio between the interface thickness and the domain size. Characteristic properties of (38) are the conservation of the order parameter [12] 
383 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
$$
\frac {\mathrm{d}}{\mathrm{d} t} \int_ {\Omega} \phi (t, \mathbf {x}) \mathrm{d} \mathbf {x} = 0, \tag {39}
$$
and a monotonic decrease in the total energy 
$$
\frac {\mathrm{d}}{\mathrm{d} t} F [ \phi ] = \int_ {\Omega} \left\{f (\phi) + \frac {C ^ {2}}{2} | \nabla \phi | ^ {2} \right\} \mathrm{d} \mathbf {x} \leqslant 0. \tag {40}
$$
![image](/16b96cbfb148f80d3df9565e13704bc728c3f9cccfd85b5ddbb7d3b543a7cc19.jpg)

/0e13cd8daa9290e986b790d669800ab5dcd01461532219772a10ca5b8b30a1f3.jpg
/fef2291a4d05c3abe32104469fd4819469ba675d2e56e9042f094c2814ffe4de.jpg
/e805358de209593b24bfa0441d0322d128cbd4cf82fe005e4589adaff834c59e.jpg

Fig. 4. Evolution of /, represented in flooded contours, at different times: (a) $t = 6 . 9 3 \times 1 0 ^ { - 4 }$ , (b) $t = 0 . 3 6 ,$ (c) t ¼ 1:16, (d) t ¼ 2:76. N ¼ 1024, c ¼ 0, /m ¼ 0, C ¼ 7:03  10-4 and $\Delta t = 4 . 9 5 \times 1 0 ^ { - 7 }$ for (a), Dt ¼ 0:0002 for (b)–(d).

384 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
We take as initial condition a random perturbation of a uniform mixture as follows: 
$$
\phi (0, \mathbf {x}) = \phi_ {\mathrm{m}} + C r (\mathbf {x}), \tag {41}
$$
where the random $r ( \mathbf { x } )$ is in ½-1; 1 and has zero mean. $\phi _ { \mathrm { m } }$ is the constant concentration of the uniform mixture. The domain is the unit square. 
In our first example we consider periodic boundary conditions, constant mobility $( \gamma = 0 )$ and $\phi _ { \mathrm { m } } = 0$ which corresponds to the well known case of spinodal phase separation controlled by bulk diffusion. We take $C = 7 . 0 3 \times 1 0 ^ { - 4 }$ and using a spatial mesh of $1 0 2 4 \times 1 0 2 4$ points we have interfacial thickness of three points. According to linear analysis (see e.g. [34]) the fastest growth rate is $1 / ( 4 C ^ { 2 } )$ . The solution quickly develops two spatial length scales, one associated with the wavelength k of the fastest growing mode and the other, the shortest one, with the transitions between phases. For $\phi _ { \mathrm { m } } = 0 , \lambda = 2 \pi \sqrt { 2 } C$ , while the phase transition layers are approximately of size C and thus a mesh size of $\mathbf { O } ( C )$ is needed. After the fast initial stage the dynamics are very slow and it takes a long time to reach a quasi-stationary state. With the second order semi-implicit scheme we can compute stably the solution and resolve both the fast initial dynamics and the slow long-time behavior, varying Dt to adjust to the dynamics, while at the same time retaining the required high spatial resolution. 
Fig. 4 shows snapshots of the solution plotted as flooded contours. The lightest phase corresponds to $\phi = 1$ and the darkest one to $\phi = - 1$ . The initially homogeneous mixture undergoes a fast separation followed by slow coarsening where typical spinodal structures can be observed. Due to the very small Cahn number and the high resolution the interfaces separating the structures appear fairly sharp. We start the computation with $\Delta t = 4 . 9 5 \times 1 0 ^ { - 7 }$ to resolve the early fast growth of solution, but we only compute with this time-step up to $t = C = 7 . 0 3 \times 1 0 ^ { - 4 }$ . For the longer time computation we use $\Delta t = 0 . 0 0 0 2$ . This timestep selection is based on accuracy as the method appears to be unconditionally stable, and any choice of $\Delta t$ produces a stable computation. 
/f4f00f679c89b88725d0fa6c5274179df428cdcaf7aaedace96fcd04fddf55d2.jpg

Fig. 5. Structure function as a function of k at five different time steps for bulk-diffusion-controlled coarsening. Time increases in the direction of vertical axes.

385 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 

As shown in [35] a two-phase morphology undergoing coarsening can be characterized by the timedependent structure function 
$$
S (\mathbf {k}, t) = \frac {1}{N} \left\langle \sum_ {\mathbf {r}} \sum_ {\mathbf {r} ^ {\prime}} \mathrm{e} ^ {- \mathrm{i} \mathbf {k} \cdot \mathbf {r}} \left[ \phi (\mathbf {r} + \mathbf {r} ^ {\prime}, t) \phi (\mathbf {r} ^ {\prime}, t) - \langle \phi \rangle^ {2} \right] \right\rangle , \tag {42}
$$
where both sums run over the lattice, N is the total number of points in the lattice, and h i stands for average over all lattice points. The normalized structure function $s ( k , t )$ is given by 
$$
s (k, t) = \frac {S (k , t)}{N \left[ \langle \phi^ {2} (\mathbf {r}) \rangle - \langle \phi \rangle^ {2} \right]} \tag {43}
$$
and we can characterize the typical length scale $R ( t )$ with the first moment of $s ( k , t )$ [25], 
$$
k _ {1} (t) = \frac {\sum k s (k , t)}{\sum s (k , t)}. \tag {44}
$$
In Fig. 5 we plot the normalized and circularly averaged structure function at five different time steps. The lines are spline fits to the simulation data. As time increases, the maximum value of the structure function increases and shifts to lower k, indicating an increase in the real-space average length scale. This is consistent with the results reported in [25]. In Fig. 6 we plot the cubic of the average domain size vs time. The straight line behavior confirms the expected cubic growth law [25]. 
/968a1ddc52c11d4a15f9a6a1396356cd1b340b100003f05231121badd9ec467e.jpg

Fig. 6. The cubic of the average domain size vs time at the late stage of bulk-diffusion-controlled coarsening (domain size characterized by $( 1 / k _ { 1 } ) )$ .

386 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 

We now consider a case of variable mobility by setting c ¼ 0:9. Fig. 7 shows the morphological evolution of the mixture for $\phi _ { \mathrm { m } } = 0$ and Cahn number $C = 0 . 0 0 1$ using a 1024  1024 resolution. This is the case of interface-diffusion-controlled coarsening that is characterized by much slower dynamics but with similar morphological patterns. These results are analogous to the ones reported in [25]. But here, with the unconditionally stable scheme, we are able use a large time step $( \Delta t = 0 . 0 1 )$ to follow the very slow coarsening dynamics. Moreover, we can resolve a thinner interface of only three mesh points, with second order time integration. 
![image](/66f3a05a39bbf4420f980807abd13404b23025ae01d5978dd027f71a457fcb60.jpg)

![image](/36fc545ebade0ec29b97e84f0d9fcb342b11d829093519d5f9c895a918b0d7eb.jpg)

![image](/5043ae1d8569f3e43bc8bd212dd1103bf2b224a96f800a08097062d205c6613d.jpg)

![image](/4a25667cbfadf4cdf43d8f7f6dbaeed1056c449d32cd9a436990de76a6432099.jpg)


Fig. 7. Variable mobility: evolution of $\phi ,$ represented in flooded contours, at different times: (a) t ¼ 0:1, (b) t ¼ 0:47, (c) t ¼ 1:4, (d) $t = 3 . 0 . ~ N = 1 0 2 4 , \gamma = 0 . 9 , \phi _ { \mathrm { m } } = 0 . 0 , C = 0 . 0 0 1 \mathrm { ~ a n d ~ } \Delta t = 0 . 0 1$ for (a)–(d).

387 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 

We turn now to two 3D simulations of pure phase separation with constant mobility $( \gamma = 0 )$ and no-flux boundary conditions, i.e., n  $\nabla \phi = 0$ and $\mathbf { n } \cdot \nabla ( f ^ { \prime } ( \phi ) - C ^ { 2 } \nabla ^ { 2 } \phi ) = 0$ (Figs. 8 and 9). We take first $\phi _ { \mathrm { m } } = 0$ and $C = 0 . 0 1$ and we render the iso-surface of separation of the two fluids at $\phi = 0$ . Fig. 8 depicts representative snapshots of the iso-surface. Notice the complexity of the patterns that cannot be extrapolated from the 2D counterpart. The simulation begins with $\Delta t = 2 . 5 \times 1 0 ^ { - 4 }$ up to $t = C = 0 . 0 1$ and for longer times $\Delta t = 0 . 0 1$ is used. Fig. 9 presents very different separation morphology. For this simulation we take $\phi _ { \mathrm { m } } = - 0 . 5$ and we render the initial stages using $\Delta t = 2 . 5 \times 1 0 ^ { - 4 }$ . The initial uniform mixture evolves into a system consisting of a large array of round particles at $t = 0 . 0 1$ . The coarsening takes place and the spherical drops grow until they coalesce. 
![image](/9fc3f578094baf562bc8cd7b7e636055b0b0b5c8cd15cc91dff53b6396ebdb08.jpg)

![image](/0d55dc2c15a35570cee5030136db278d4d80c222a16c0bb57a781cb450b56e4c.jpg)

![image](/e57cdc9b936e9612f86d98e76ff91c53e1c2ee85b6e15a17ba55d8b6d59f8a4f.jpg)

![image](/c5c5b3a1eb2c22b0a0d64027b0268fec860a9602e79bb0bf4ef44000e45e5b18.jpg)


Fig. 8. Evolution of /, represented by the iso-surfaces of separation of the two fluids at $\phi = 0 . 0 $ at different times: $( \mathrm { a } ) \ t = 0 . 0 3 7 5 ,$ (b) $t = 6 . 5 , \mathrm { ~ ( c ) ~ } t = 1 4 . 0 , \mathrm { ~ ( d ) ~ } t = 9 4 . 0 , N = 2 5 6 , \gamma = 0 . 0 , \phi _ { \mathrm { m } } = 0 . 0 \mathrm { ~ a n d ~ } \Delta t = 2 . 5 \times 1 0 ^ { - 4 } \mathrm { ~ f o r ~ ( a ) } , \Delta t = 0 . 0 1 \mathrm { ~ f o r ~ ( b ) - ( d ) ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t = 0 . 0 0 \mathrm { ~ a n d ~ } t o r i a = 0 .$ ).

388 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 

![image](/cf1bf8948564e553cfd2daaab2d0a1e5b16e6984ca857908ed7669761f30d7e6.jpg)

![image](/6de919d512ed01065d44fa7748e88ca8edd428f5d36e173674a66940f2becfe3.jpg)

![image](/8878d0b2f526ebb249c25aa7fe1ed1769c59632b9ba98e468634cbd383a060a2.jpg)

![image](/61c15e3836d0dceb244262016ed33bf1fbbe2c3a7928410973efe38a2f877753.jpg)


Fig. 9. Evolution of $\phi ,$ represented by the iso-surfaces of separation of the two fluids at $\phi = 0 . 0$ , at different times: $( \mathrm { a } ) \ t = 0 . 1 0 3 7 5 ,$ (b) $t = 0 . 1 0 4 7 5 , ( \mathrm { c } ) \ t = 0 . 1 0 6 2 5 , ( \mathrm { d } ) \ t = 0 . 1 5 . N = 2 5 6 , \gamma = 0 . 0 , \phi _ { \mathrm { m } } = - 0 . 5 \ \mathrm { a n d } \ \Delta t = 2 . 5 \times 1 0 ^ { - 4 } \ \mathrm { f m }$ or (a)–(d).

389 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
/57f1441e9e75f9f9df7315c19e016f71c3df3598ea6913ae516fcdd07ba18c9d.jpg
/7e81537ac15747bbe63b5ea2ee4072c05781bc3ba517438ea043066532b5538d.jpg

Fig. 10. Behavior of the mean $\phi _ { \mathrm { m } }$ and of the energy F ð/Þ in time for the semi-implicit scheme for $\phi _ { \mathrm { m } } ( t = 0 ) = 0$ .

Fig. 10 shows the time behavior of the phase field mean and the energy for a spinodal decomposition with a resolution of $1 2 8 \times 1 2 8 \times 1 2 8$ . We find that the mean is preserved within 3–4 digits and the energy decreases monotonically (and smoothly) throughout the entire computation as required by (40). 
# 4.3. 2D and 3D phase separation and pattern formation in a channel under shear
We consider phase separation (spinodal decomposition) of a density-matched binary fluid mixtures in a channel under shear. As we will see, linear shear plays a crucial role in the morphology and evolution of the patterns. The initial conditions are a random perturbation around the uniform concentration $\phi = 0 .$ . Figs. 11(a) and (b) show 2D results at two shear rates with the top lid and bottom lid moving horizontally in opposite directions. The shear rate is defined as $s r = U _ { \mathrm { c } } / h$ where h is the distance between the two plates. 
Since we impose a fixed geometry, $s r \propto U _ { \mathrm { c } }$ . Thus, to change the shear rate we need to change Pe, Re and Ca as they all contain $U _ { \mathrm { c } }$ . The flow in Fig. 11(b) has five times the shear rate as that in Fig. 11(a). We notice that after a transient stage characterized by the formation of patterns in the mixture under the influence of the Cahn–Hilliard term (spinodal decomposition), the domains get elongated into long layers against their intrinsic surface tension instabilities. Moreover, the patterns formed in the early stage are quite different when the shear velocity increases, and the number of layers in the late stage increases when the shear rate is higher. This behavior is in accordance with experiments reported in [36] and simulations in [9]. 
390 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
![image](/a73b2ba8d6a13475b3c567c5d0ffb60d2229147ff7c594a6a964a06dc5a6cb4a.jpg)


Fig. 11. 2D spinodal decomposition in a channel under shear: (a) (first column) $P e = 7 . 5 , R e = 0 . 1$ , and Ca ¼ 0:5; (b) (second column) $P e = 3 7 . 5 , R e = 0 . 5$ , and Ca ¼ 2:5. N ¼ 256 and L ¼ 355.


391 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 

In Fig. 12 we consider a 3D simulation in the presence of shear. Here the structures are much more complex with strings forming. String-like structures have been observed in polymer blends which are thermodynamically near a phase transition point [37,38] and in immiscible viscoelastic systems in complex flow fields [39] and in dispersed droplets [39]. There is great current interest in micro- and nanolengthscale technologies in which polymer blends could play an important role. For example, if we create strings with a conductive material in an insulating matrix with good mechanical properties, then one could produce wires. In other ways [39] it might be possible to manufacture ultrathin materials of high one-dimensional strength or scaffolds. A detailed study of the string process formation with our numerical procedure is under way and it will be reported elsewhere. The methodology presented here appears quite promising for the design and analysis of multiphase and complex fluid formulations. 
![image](/d9b64fbbf4c950e258fb96c03fe20bbc3f7d606f8a2809bcaf1cfa053b7430ed.jpg)

![image](/62ca26232e5f3434b9d9fdab68147584543ecc2160df3149c8e3be56c040f3da.jpg)

![image](/16939d64eaa28d996182e903aaf977e0cee098accabeb58fbc53b6dead65be72.jpg)

![image](/5adaac96ba5765f758ef80c0adaa3fe20041f54f335a91e1ab63244361507df9.jpg)


Fig. 12. 3D spinodal decomposition in a channel under shear, Pe ¼ 10, Re ¼ 0:5, Ca ¼ 5, N ¼ 128 and L ¼ 178: (a) t ¼ 800, (b) t ¼ 2500, (c) t ¼ 5000, (d) t ¼ 57,100.

392 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 

# 5. Concluding remarks
An accurate and efficient numerical method for computing phase ordering kinetics coupled with fluid dynamics was presented. The numerical method is a time-split scheme that combines a novel semi-implicit discretization for the convective Cahn–Hilliard equation with a ‘‘stiffly stable’’ time-discretization of the projection method for the Navier–Stokes equations. The numerical method is robust and has minimal cost. Some of the capabilities of the method were illustrated with numerical examples in two and three dimensions, including the technologically important problem of phase separation under shear flow. In particular, the 3D simulations in the presence of shear flow reveal rich and complex structures, including strings. The method can be extended to general geometries through the use of other spatial high order discretizations such as in spectral element methods, while retaining the same characteristic of stability and efficiency. The type of discretizations presented here also offer great promise for the computation of complex fluid systems such as polymeric flows. 
# Acknowledgements
We thank G. Fredrickson, C. Garcia-Cervera, D. Jacqmin, G. Leal, J. Lowengrub, and T.Y. Hou for helpful discussions and M. Frigo for providing us with the latest FFTW libraries [40]. This work was partially supported by the National Aeronautics and Space Administration, Microgravity Research Division under the Contract No. NAG3-2414. H.D.C. acknowledges partial support from the Academic Senate Junior Faculty Research Award. 
# Appendix A. The spatial discretization
# A.1. Helmholtz equation
We rewrite the third step of the projection method (24) (Helmholtz equation) as 
$$
\begin{array}{l} \frac {2 R e}{\theta_ {\max}} \frac {\mathbf {u} ^ {*}}{\Delta t} - \nabla^ {2} \mathbf {u} ^ {*} = \frac {2 R e}{\theta_ {\max}} \left(\frac {\mathbf {u} ^ {n}}{\Delta t} - \mathbf {u} ^ {n} \cdot \nabla \mathbf {u} ^ {n}\right) + \frac {2}{\theta_ {\max} C a} \mu (\phi^ {n + 1}) \nabla \phi^ {n + 1} \\ + \frac {2}{\theta_ {\max}} \left[ \nabla \cdot \theta^ {n + 1} (\nabla \mathbf {u} ^ {n} + (\nabla \mathbf {u} ^ {n}) ^ {\mathrm{T}}) - \frac {\theta_ {\max}}{2} \nabla^ {2} \mathbf {u} ^ {n} \right]. \tag {A.1} \\ \end{array}
$$
Since we have periodic boundary conditions in the horizontal direction we can Fourier transform to obtain (dropping the asterisk) 
$$
\left(\frac {2 R e}{\theta_ {\max} \Delta t} + k _ {x} ^ {2} + k _ {y} ^ {2}\right) \hat {\mathbf {u}} - \hat {\mathbf {u}} ^ {\prime \prime} = \widehat {\boldsymbol {\Omega}} (k _ {x}, k _ {y}, z), \tag {A.2}
$$
393 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
where the prime denotes derivative with respect to $z , \Omega ( x , y , z )$ is the right-hand side of (A.1) and the caret stands for the 2D Fourier transform in the streamwise direction. Thus (i is the z-index of the grid) 
$$
\hat {\mathbf {u}} _ {i} ^ {\prime \prime} = k ^ {2} \hat {\mathbf {u}} _ {i} - \widehat {\boldsymbol {\Omega}} _ {i}, \tag {A.3}
$$
where $k ^ { 2 } = ( R e / \theta _ { \mathrm { m a x } } \Delta t ) + k _ { x } ^ { 2 } + k _ { \nu } ^ { 2 }$ . An eighth order finite difference compact scheme discretization of (A.6), as we will see in Section A.3, yields the pentadiagonal system 
$$
C \hat {\boldsymbol {u}} _ {i - 2} + B \hat {\boldsymbol {u}} _ {i - 1} + A \hat {\boldsymbol {u}} _ {i} + B \hat {\boldsymbol {u}} _ {i + 1} + C \hat {\boldsymbol {u}} _ {i + 2} = \beta \widehat {\boldsymbol {\Omega}} _ {i - 2} + \alpha \widehat {\boldsymbol {\Omega}} _ {i - 1} + \widehat {\boldsymbol {\Omega}} _ {i} + \alpha \widehat {\boldsymbol {\Omega}} _ {i + 1} + \beta \widehat {\boldsymbol {\Omega}} _ {i + 2}, \tag {A.4}
$$
where $A = k ^ { 2 } + ( b / 2 + 2 a ) / ( \Delta z ) ^ { 2 } , B = \alpha k ^ { 2 } - a / ( \Delta z ) ^ { 2 } , C = \beta k ^ { 2 } - b / ( 2 \Delta z ) ^ { 2 }$ . The parameters a, b, a and $\beta$ (given in Section A.3) are chosen to achieve formal eighth order accuracy [5]. 
# A.2. Poisson equation
Since we have periodic boundary conditions in the horizontal direction we can Fourier transform the Poisson equation (28) to obtain 
$$
(k _ {x} ^ {2} + k _ {y} ^ {2}) \hat {\boldsymbol {p}} - \hat {\boldsymbol {p}} ^ {\prime \prime} = \widehat {\Xi} (k _ {x}, k _ {y}, z), \tag {A.5}
$$
where the prime denotes derivative with respect to $z , \Xi ( x , y , z )$ is the right-hand side of (28) and the caret stands for the 2D Fourier transform in the streamwise direction. Thus (i is the z-index of the grid) 
$$
\hat {p} _ {i} ^ {\prime \prime} = k ^ {2} \hat {p} _ {i} - \widehat {\Xi} _ {i}, \tag {A.6}
$$
where $k ^ { 2 } = k _ { x } ^ { 2 } + k _ { \nu } ^ { 2 }$ . An eighth order finite difference compact scheme discretization of (A.6) yields the pentadiagonal system: 
$$
C \hat {p} _ {i - 2} + B \hat {p} _ {i - 1} + A \hat {p} _ {i} + B \hat {p} _ {i + 1} + C \hat {p} _ {i + 2} = \beta \widehat {\Xi} _ {i - 2} + \alpha \widehat {\Xi} _ {i - 1} + \widehat {\Xi} _ {i} + \alpha \widehat {\Xi} _ {i + 1} + \beta \widehat {\Xi} _ {i + 2}, \tag {A.7}
$$
where $A = k ^ { 2 } + ( b / 2 + 2 a ) / ( \Delta z ) ^ { 2 } , B = \alpha k ^ { 2 } - a / ( \Delta z ) ^ { 2 } , C = \beta k ^ { 2 } - b / ( 2 \Delta z ) ^ { 2 }$ and the parameters $a , b ,$ a and $\beta$ as provided in the next section. The Neumann boundary condition (29), applied at i ¼ 1 and $i = N _ { z }$ , is implemented via second order approximations: 
$$
\frac {3}{2 \Delta z} \hat {p} _ {1} - \frac {2}{\Delta z} \hat {p} _ {2} + \frac {1}{2 \Delta z} \hat {p} _ {3} = - \hat {p} _ {1} ^ {\prime}, \tag {A.8}
$$
$$
- \frac {1}{2 \Delta z} \hat {p} _ {N _ {z} - 2} + \frac {2}{\Delta z} \hat {p} _ {N _ {z} - 1} - \frac {3}{2 \Delta z} \hat {p} _ {N _ {z}} = - \hat {p} _ {N _ {z}} ^ {\prime} \tag {A.9}
$$
where $\hat { p } _ { 1 } ^ { \prime }$ and $\hat { p } _ { N _ { z } } ^ { \prime }$ are the ðx; yÞ-transforms of the z-derivatives at the walls. The pentadiagonal matrix for this linear system is well-conditioned, except for $k = 0 .$ , in which case it is singular. This situation arises because, with Neumann conditions at both ends of the domain, the solution for the pressure is nonunique since pressure is only defined within a constant. Rewriting the momentum equation (13) at the wall with the use of the incompressibility condition (14) 
$$
\frac {\partial p}{\partial z} = \theta \frac {\partial^ {2} w}{\partial z ^ {2}}, \tag {A.10}
$$
where w is the wall normal component of $\mathbf { u } = ( u , v , w )$ . Fourier transforming (A.10) in x- an y-directions and using the incompressibility condition (14) we get 
394 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
$$
\frac {\partial \hat {p}}{\partial z} = \mathrm{i} k _ {x} \theta \frac {\partial \hat {u}}{\partial z} + \mathrm{i} k _ {y} \theta \frac {\partial \hat {v}}{\partial z}, \tag {A.11}
$$
i.e., for the singular case $k _ { x } = k _ { y } = 0$ the two Neumann conditions at both ends of the domain (29) reduce to 
$$
\frac {\partial \hat {p}}{\partial z} = 0. \tag {A.12}
$$
To solve for the case $k = k _ { x } = k _ { y } = 0$ we use the cosine transform since it automatically satisfy (A.12); we then deal with a third wave number $k _ { z }$ and for the case $k _ { x } = k _ { y } = k _ { z } = 0$ we set the solution as a constant. This is inconsequential since as noted before pressure is only defined within a constant. 
# A.3. Finite difference compact schemes
For the first derivative in the z (wall normal direction) we use the compact approximation scheme [5] 
$$
\beta \hat {\boldsymbol {u}} _ {i - 2} ^ {\prime} + \alpha \hat {\boldsymbol {u}} _ {i - 1} ^ {\prime} + \hat {\boldsymbol {u}} _ {i} ^ {\prime} + \alpha \hat {\boldsymbol {u}} _ {i + 1} ^ {\prime} + \beta \hat {\boldsymbol {u}} _ {i + 2} ^ {\prime} = b \frac {\hat {\boldsymbol {u}} _ {i + 2} - \hat {\boldsymbol {u}} _ {i - 2}}{4 \Delta z} + a \frac {\hat {\boldsymbol {u}} _ {i + 1} - \hat {\boldsymbol {u}} _ {i - 1}}{2 \Delta z}, \tag {A.13}
$$
where the prime denotes derivative with respect z and the caret stands for the 2D Fourier transform in the streamwise direction. The optimized coefficients for an eighth order compact stencil are $\alpha = 4 / 9 , \beta = 1 / 3 6$ , $a = 4 0 / 2 7 , b = 2 5 / 5 4$ . For the points neighboring boundaries i ¼ 2 and $i = N _ { z } - 1$ we use a fourth order scheme with $\alpha = 1 / 4 , \beta = 0 , a = 3 / 2 , b = 0$ . 
The compact approximation schemes for the boundaries i ¼ 1 and $N _ { z }$ are 
$$
\hat {u} _ {1} ^ {\prime} + 2 \hat {u} _ {2} ^ {\prime} = \frac {1}{\Delta z} \left(- \frac {5}{2} \hat {u} _ {1} + 2 \hat {u} _ {2} + \frac {1}{2} \hat {u} _ {3}\right), \tag {A.14}
$$
$$
\hat {\boldsymbol {u}} _ {N _ {z}} ^ {\prime} + 2 \hat {\boldsymbol {u}} _ {N _ {z} - 1} ^ {\prime} = \frac {1}{\Delta z} \left(\frac {5}{2} \hat {\boldsymbol {u}} _ {N _ {z}} - 2 \hat {\boldsymbol {u}} _ {N _ {z} - 1} - \frac {1}{2} \hat {\boldsymbol {u}} _ {N _ {z} - 2}\right). \tag {A.15}
$$
These are third order schemes [5,18]. 
To approximate the second derivative we use 
$$
\beta \hat {u} _ {i - 2} ^ {\prime \prime} + \alpha \hat {u} _ {i - 1} ^ {\prime \prime} + \hat {u} _ {i} ^ {\prime \prime} + \alpha \hat {u} _ {i + 1} ^ {\prime \prime} + \beta \hat {u} _ {i + 2} ^ {\prime \prime} = b \frac {\hat {u} _ {i + 2} - 2 \hat {u} _ {i} + \hat {u} _ {i - 2}}{4 (\Delta z) ^ {2}} + a \frac {\hat {u} _ {i + 1} - 2 \hat {u} _ {i} + \hat {u} _ {i - 1}}{(\Delta z) ^ {2}} \tag {A.16}
$$
and the optimized coefficients for an eighth-order compact stencil are $\alpha = 3 4 4 / 1 1 7 9 , \beta = 2 3 / 2 3 5 8$ ; $a = 3 2 0 / 3 9 3 , b = 3 1 0 / 3 9 3$ ; while for the points i ¼ 2 and $i = N _ { z } - 1$ we use a fourth order scheme with $\alpha = 1 / 1 0 , \beta = 0 , a = 6 / 5 . b = 0 $ ;. For the boundaries we choose 
$$
\hat {u} _ {1} ^ {\prime \prime} + 1 1 \hat {u} _ {2} ^ {\prime \prime} = \frac {1}{(\Delta z) ^ {2}} (1 3 \hat {u} _ {1} - 2 7 \hat {u} _ {2} + 1 5 \hat {u} _ {3} - \hat {u} _ {4}), \tag {A.17}
$$
$$
\hat {\boldsymbol {u}} _ {N _ {z}} ^ {\prime \prime} + 1 1 \hat {\boldsymbol {u}} _ {N _ {z} - 1} ^ {\prime \prime} = \frac {1}{(\Delta z) ^ {2}} (1 3 \hat {\boldsymbol {u}} _ {N _ {z}} - 2 7 \hat {\boldsymbol {u}} _ {N _ {z} - 1} + 1 5 \hat {\boldsymbol {u}} _ {N _ {z} - 2} - \hat {\boldsymbol {u}} _ {N _ {z} - 3}). \tag {A.18}
$$
These are third order accurate schemes with a truncation error 10 times smaller than that of the analogous explicit one (see [5,18]). 
395 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 
# A.4. Properties of cosine transforms
We define the Fourier cosine transforms of a function f ðxÞ as 
$$
C [ f (x) ] = \frac {2}{\pi} \int_ {0} ^ {\infty} f (x) \cos \omega x d x. \tag {A.19}
$$
For the second derivatives we have 
$$
C \left[ \frac {\mathrm{d} ^ {2} f}{\mathrm{d} x ^ {2}} \right] = - \frac {2}{\pi} \frac {\mathrm{d} f}{\mathrm{d} x} (0) - \omega^ {2} C [ f ] \tag {A.20}
$$
under the hypothesis of compact support for $f ( x )$ and $f ^ { \prime } ( x )$ . 
For the fourth derivative we have 
$$
C \left[ \frac {\mathrm{d} ^ {4} f}{\mathrm{d} x ^ {4}} \right] = - \frac {2}{\pi} \left(\frac {\mathrm{d} ^ {3} f}{\mathrm{d} x ^ {3}} (0) - \omega^ {2} \frac {\mathrm{d} f}{\mathrm{d} x} (0)\right) + \omega^ {4} C [ f ] \tag {A.21}
$$
under the hypothesis of compact support for $f ( x ) , f ^ { \prime } ( x )$ and $f ^ { \prime \prime \prime } ( x )$ . We used these properties to solve the Cahn–Hilliard and the Poisson equation for the case $k _ { x } = k _ { y } = 0 ;$ for the Cahn–Hilliard equation $f ^ { \prime } ( 0 ) = f ^ { \prime \prime \prime } ( 0 ) = 0$ due to the boundary conditions (4) and for the Poisson equation $f ^ { \prime } ( 0 ) = 0$ due to the boundary condition (A.12). Note that these properties hold for the discrete transforms as well. 
# References


[1] J.W. Cahn, J.E. Hilliard, Free energy of a nonuniform system I, J. Chem. Phys. 28 (1958) 258. 




[2] J.W. Cahn, J.E. Hilliard, Free energy of a nonuniform system III, J. Chem. Phys. 31 (1959) 688. 




[3] P.C. Hohenberg, B.I. Halperin, Theory of dynamic critical phenomena, Rev. Mod. Phys. 49 (3) (1977) 435. 




[4] J. Lowengrub, L. Truskinovsky, Quasi-incompressible Cahn–Hilliard fluids and topological transitions, Proc. R. Soc. Lond. A 454 (1998) 2617. 




[5] S.K. Lele, Compact finite difference schemes with spectral-like resolution, J. Comput. Phys. 103 (1992) 16. 




[6] R. Chella, V. Vinals, Mixing of a two-phase fluid by a cavity flow, Phys. Rev. E 53 (1996) 3832. ~ 




[7] D. Jacqmin, Calculation of two phase Navier Stokes flows using phase-field modeling, J. Comput. Phys. 115 (1999) 96. 




[8] V.M. Kendon, M.E. Cates, I.P. Barraga, J.-C. Desplat, P. Blandon, Inertial effects in three-dimensional spinodal decomposition of a symmetric binary fluid mixture: a lattice Boltzmann study, J. Fluid Mech. 440 (2001) 147. 




[9] Y. Wu, H. Skrdla, T. Lookman, S. Chen, Spinodal decomposition in binary fluids under shear flow, Physica A 239 (1997) 428– 436. 




[10] O. Penrose, P. Fife, Thermodynamically consistent models of phase-field type for the kinetics of phase transitions, Physica D 43 (1990) 44. 




[11] P.W. Bates, P.C. Fife, The dynamics of nucleation for the Cahn–Hilliard equation, SIAM J. Appl. Math. 53 (1993) 990. 




[12] C.M. Elliot, The Cahn–Hilliard model for the kinetics of phase separation, in: J.F. Rodrigues (Ed.), Mathematical Models for Phase Change Problems, International Series of Numerical Mathematics, vol. 88, Bikh€auser, Basel, 1989, pp. 35–72. 




[13] M.E. Gurtin, D. Polignone, J. Vinals, Two-phase binary fluids and immiscible fluids described by an order parameter, Math. ~ Models Meth. Appl. Sci. 6 (6) (1996) 815. 




[14] A.J. Bray, Theory of phase-ordering kinetics, Adv. Phys. 43 (3) (1994) 357–459. 




[15] J.D. van der Waals, The thermodynamic theory of capillarity flow under the hypothesis of a continuous variation of density (in Dutch), Verhandel/Konink. Akad. Weten. 1 (1879) 8. 




[16] J.S. Langer, M. Baron, H. Miller, New computational method in theory of spinodal decomposition, Phys. Rev. A 11 (4) (1975) 1417. 




[17] A.W. Cook, P.E. Dimotakis, Transition stages of Rayleigh–Taylor instability between miscible fluids, J. Fluid Mech. 443 (2001) 69. 




[18] J.C. Buell, A hybrid numerical-method for 3-dimensional spatially-developing free shear flows, J. Comput. Phys. 95 (2) (1991) 313–338. 


396 
V.E. Badalassi et al. / Journal of Computational Physics 190 (2003) 371–397 


[19] J.J. Douglas, T. Dupont, Alternating-direction Galerkin methods on rectangles, in: B. Hubbard (Ed.), SYNSPADE-1970, Numerical Solution of Partial Differential Equations-II, Academic Press, New York, 1971, pp. 133–213. 




[20] D. Gottlieb, S.A. Orszag, in: Numerical Analysis of Spectral Methods: Theory and Applications, CBMSNSF Regional Conference Series in Applied Mathematics, SIAM Press, 1977. 




[21] H.D. Ceniceros, A semi-implicit moving mesh method for the focusing Schrodinger equation, Commun. Pure Appl. Anal. 1 (2002) € 1. 




[22] U.M. Ascher, S.J. Ruuth, B.T.R. Wetton, Implicit–explicit methods for time dependent partial differential equations, SIAM J. Numer. Anal. 32 (1995) 797. 




[23] D.J. Eyre, An unconditionally stable one-step scheme for gradient systems (preprint). 




[24] P. Smereka, Semi-implicit level set methods for motion by mean curvature and surface diffusion (preprint). 




[25] J. Zhu, L.-Q. Chen, J. Shen, V. Tikare, Coarsening kinetics from a variable-mobility Cahn–Hilliard equation: application of a semi-implicit Fourier spectral method, Phys. Rev. E 60 (4) (1999) 3564–3572. 




[26] G.E. Karniadakis, M. Israeli, S.A. Orszag, High-order splitting methods for the incompressible Navier–Stokes equations, J. Comput. Phys. 97 (1991) 414–443. 




[27] M. Sussman, P. Smereka, S. Osher, A level set approach for computing solutions to incompressible two-phase flow, J. Comput. Phys. 114 (1994) 146–159. 




[28] J.U. Brackbill, D.B. Kothe, C. Zemach, A continuum method for modeling surface tension, J. Comput. Phys. 100 (1992) 335. 




[29] G.I. Taylor, The formation of emulsions in definable fields of flows, Proc. R. Soc. Lond. A 146 (1934) 501–523. 




[30] M. Shapira, S. Haber, Low Reynolds number motion of a droplet in shear flow including wall effects, Int. J. Multiphase Flow 16 (2) (1990) 305–321. 




[31] J.M. Rallison, The deformation of small viscous drops and bubbles in shear flow, Annu. Rev. Fluid Mech. 16 (1984) 45–66. 




[32] T. Roths, C. Friedrich, M. Marth, J. Honerkamp, Dynamics and rheology of the morphology of immiscible polymer blends—on modeling and simulation, Rheologica Acta 41 (3) (2002) 211–222. 




[33] E.G. Puckett, J.B.B.A.S. Almgren, D.L. Marcus, W.J. Rider, A high-order projection method for tracking fluid interfaces in variable density incompressible flows, J. Comput. Phys. 130 (2) (1997) 269–282. 




[34] M. Copetti, C. Elliot, Kinetics of phase decomposition process: numerical solutions to the Cahn–Hilliard equation, Mater. Sci. Technol. 6 (1990) 273. 




[35] A. Chakrabarti, R. Toral, J.D. Gunton, Late-stage coarsening for off-critical quenches: scaling functions and the growth law, Phys. Rev. E 47 (1993) 3025–3038. 




[36] A. Onuki, Phase transitions of fluids in shear flow, J. Phys.: Condes. Matter 9 (1997) 6119. 




[37] A. Frischknecht, Stability of cylindrical domains in phase-separating binary fluids in shear flow, Phys. Rev. E 58 (3) (1998) 3495– 3514. 




[38] T. Hashimoto, K. Matsuzaka, E. Moses, A. Onuki, String phase in phase-separating fluids under shear flow, Phys. Rev. Lett. 74 (1) (1995) 126–129. 




[39] K.B. Migler, String formation in sheared polymer blends: coalescence, breakup, and finite size effects, Phys. Rev. Lett. 86 (6) (2001) 1023–1026. 




[40] M. Frigo, S.G. Johnson, FFTW: an adaptive software architecture for the fft, ICASSP Conf. Proc. 3 (1998) 1381–1384. 


397 