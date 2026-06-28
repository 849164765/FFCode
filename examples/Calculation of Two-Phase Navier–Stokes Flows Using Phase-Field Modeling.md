# Calculation of Two-Phase Navier–Stokes Flows Using Phase-Field Modeling

David Jacqmin 

NASA Lewis Research Center, Cleveland, Ohio 44135 

Received February 15, 1999; revised July 7, 1999 

Phase-field models provide a way to model fluid interfaces as having finite thickness. This can allow the computation of interface movement and deformation on fixed grids. This paper applies phase-field modeling to the computation of two-phase incompressible Navier–Stokes flows. The Navier–Stokes equations are modified by the addition of the continuum forcing $- C \vec { \nabla } \phi$ , where $C$ is the composition variable and $\phi$ is $C \mathrm { { : } }$ chemical potential. The equation for interface advection is replaced by a continuum advective-diffusion equation, with diffusion driven by $C \mathrm { { : } }$ chemical potential gradients. The paper discusses how solutions to these equations approach those of the original sharp-interface Navier–Stokes equations as the interface thickness ² and the diffusivity both go to zero. The basic flow-physics of phase-field interfaces is discussed. Straining flows can thin or thicken an interface and this must be resisted by a high enough diffusion. On the other hand, too large a diffusion will overly damp the flow. These two constraints result in an upper bound for the diffusivity of $O ( \epsilon )$ and a lower bound of $O ( \epsilon ^ { 2 } )$ . Within these two bounds, the phasefield Navier–Stokes equations appear to generate an $O ( \epsilon )$ error relative to the exact sharp-interface equations. An ${ \cal O } ( h ^ { 2 } / \epsilon ^ { 2 } )$ numerical method is introduced that is energy conserving in the sense that creation of interface energy by convection is always balanced by an equal decrease in kinetic energy caused by surface tension forcing. An ${ \cal O } ( h ^ { 4 } / \epsilon ^ { 4 } )$ compact scheme is introduced that takes advantage of the asymptotic, comparatively smooth, behavior of the chemical potential. For $O ( \epsilon )$ accurate phasefield models the optimum path to convergence for this scheme appears to be $\epsilon \propto h ^ { 4 / 5 }$ . The asymptotic rate of convergence corresponding to this is $O ( h ^ { 4 / 5 } )$ but results at practical resolutions show that the practical convergence of the method is generally considerably faster than linear. Extensive analysis and computations show that this scheme is very effective and accurate. It allows the accurate calculation of two-phase flows with interfaces only two cells wide. Computational results are given for linear capillary waves and for Rayleigh–Taylor instabilities. The first set of computations is compared to exact solutions of the diffuse-interface equations and of the original sharp-interface equations. The Rayleigh–Taylor computations test the ability of the method to compute highly deforming flows. These flows include near-singular phenomena such as interface coalescences and breakups, contact line movement, and the formation and breakup of thin wall-films. Grid-refinement studies are made and rapid convergence is found for macroscopic flow features such as instability growth rate and propagation speed, wavelength, and the general physical characteristics of the instability and mass transfer rates. 

# 1. INTRODUCTION

Diffuse interface models provide a way of modeling interfacial forces as continuum forces, the effect being that delta-function forces and discontinuities at an interface are smoothed by distributing them over thin but numerically resolvable layers. Such models have attracted much interest recently because of their advantages for making numerical calculations. Diffuse-interface models for multiphase Navier–Stokes flow are much easier to solve than the exact equations because calculations can be done on fixed grids—diffuse interfaces simply propagate through the grids—while calculations of the exact sharp-interface equations generally require adaptive, interface fitting grids. Interface fitting grids are impractical for flows involving coalescing or splitting phases or, in general, for 3-D flows. Diffuse-interface flow models can be used to calculate these flows and many others that are currently impossible for sharp-interface solvers. Their ease of use compensates for their relatively low accuracy. 

There are currently three main types of diffuse-interface models, a tracking/distributed force model introduced by Unverdi and Tryggvason [23], the continuum surface force method (Brackbill et al. [4], Lafaurie et al. [14], Kothe et al. [13]), and phase-field (or mean-field) based models (Anderson and McFadden [1], Antonovskii [2], Chella and Vi˜nals [6], Jacqmin [8, 9], Jasnow and Vi˜nals [11], and Nadiga and Zaleski [17]). 

The method of Unverdi and Tryggvason tracks interfaces by following the advection of control points. These points mark the smeared interface’s center. The interfaces are further defined by connecting the control points by curves or line segments (in 2-D) or triangular surfaces (in 3-D). Surface tension forces are calculated from the control point positions and distributed to the fixed grids. Changes in fluid properties across the interface are smoothed so as to take place over several grid cells. 

The continuum surface force model of Brackbill et al. uses a continuum variable, such as a color function or density, to mark each phase. The local surface tension forcing is set equal to the local gradient of the continuum variable times its field curvature times the surface tension. The total forcing on the fluid through an interface is thus proportional to the interface’s gradient-weighted curvature. The model has been applied using volume-of-fluid (Lafaurie et al. [14], Rider et al. [20], Kothe et al. [13]), TVD (Jacqmin [7], and level-set (Sussman et al. [21, 22]) methods. 

The above methods are based on models of surface tension forces. Phase-field methods are based on models of fluid free energy. The simplest model of free energy density that gives two phases is 

$$
f = \frac {1}{2} \alpha | \vec {\nabla} C | ^ {2} + \beta \Psi (C) \tag {1.1}
$$

a formulation that goes back to van der Waals [24]. The first term is gradient energy, the second bulk energy. Two phases are possible if 9 has two minima. Interfaces separating two phases are $O ( { \sqrt { \alpha / \beta } } )$ in width and have a surface tension proportional to $\sqrt { \alpha \beta }$ . The surface tension forcing on the fluid is derived variationally from its energy density field. Numerical implementations of phase-field models are able to use conventional advection routines because interface profiles can be maintained against distortion by the use of highorder energy-downgradient anti-diffusion. 

Each of the above models has its advantages and disadvantages. Unverdi and Tryggvason’s method has so far met with the most success—a number of significant research results have been generated using it. Its chief drawbacks are that it requires intervention to handle topological changes, that it doesn’t conserve mass or volume, and that it can be difficult to use for three-dimensional calculations because of the need then to utilize adaptive surface grids. The CSF model handles topological changes well and it can be implemented so as to conserve mass or volume. VOF-CSF methods can be difficult to implement in three dimensions. These methods also have some instability problems and convergence issues that are not yet understood theoretically. 

Phase-field methods appear to have several potential advantages over the VOF-CSF approach. Because phase-field models allow the use of standard advection techniques they are relatively easy to implement in three dimensions, with unstructured grids, or using finite element techniques. It is easy to generate phase-field numerical implementations that are dissipative of energy, and that therefore are free of parasitic flows. So far, however, the phase-field method has fallen short in a very important respect. Phase-field interface structure is important in determining interface energy and thus surface tension. Because of the need to calculate this structure, numerical phase-field interfaces have usually been made wide, typically four to eight cells. Wide interfaces exacerbate other problems of the phase-field method. For example, many phase-field models exhibit curvaturedependent solubilities that are proportional to interface thickness. Also, wider interfaces require stronger anti-diffusion in order to keep them from being distorted by advective straining. 

The main purpose of this paper is to introduce a method that allows the use of much thinner interfaces. The asympotics of convected phase-field interfaces are outlined and it is shown how to take advantage of this asymptotics to derive simple, high-order, compact reconstruction and convection schemes. In many practical cases these schemes allow the accurate and useful calculation of phase-field convection with interfaces that are less than two cells wide. In order to lay the groundwork for this approach, the paper first discusses the convergence of phase-field modeling. This requires a discussion of both physics and numerics. The convergence of phase-field numerical calculations is dependent both on the accuracy of the phase-field model and on the accuracy of the numerical methods used to calculate the model. 

The paper proceeds as follows. The next section gives a brief introduction to the two-phase phase-field Navier–Stokes equations. The third and fourth sections look quickly at phasefield interfaces and at the effects of convection and model diffusion on those interfaces. The sixth section discusses a very simple second-order implementation of the equations that conserves energy. The seventh section discusses some fourth order methods that take advantage of interface asymptotics. The eighth, ninth, and tenth sections discuss convergence issues. The eighth section points out grid effects on interface energy calculations, the ninth discusses one-dimensional convection, and the tenth discusses linear capillary waves. Section 11 discusses a “real-world” fully nonlinear problem: the computation of Rayleigh–Taylor instabilities. 

# 2. CONTINUUM INTERFACE ENERGETICS AND EQUATIONS

A general model for an isothermal two-phase fluid’s free energy density is 

$$
f = \frac {1}{q} \alpha | \vec {\nabla} C | ^ {q} + \beta \Psi (C), \tag {2.1}
$$

C is a “measure” of phase. The free energy density is made up of two components. The first is the gradient energy $\frac { 1 } { q } \alpha | \vec { \nabla } C | ^ { q }$ and the second is the bulk energy density $\beta \Psi ( C )$ . $\Psi ( C )$ models the fluid components’ immiscibility. It has two minima corresponding to the fluids’ two stable phases. The case $q = 1 , \alpha = \sigma , \sigma$ being the surface tension, and $\beta = 0$ gives the free energy density for the CSF model. With phase-field methods, $q$ is set to $2 ,$ , $\alpha$ is set to $O ( \epsilon )$ , and $\beta$ is set to $O ( 1 / \epsilon )$ . This choice of parameters produces phase-field interfaces with $O ( \epsilon )$ thickness and $O ( 1 )$ surface tension. The most frequently used and simplest $\Psi ( C )$ is $( C + 1 / 2 ) ^ { 2 } ( C - 1 / 2 ) ^ { 2 }$ , which has a peak of high energy at $C = 0$ and minima at $C _ { \mathrm { b u l k p h a s e } } = \pm 1 / 2$ . This is the $\Psi ( C )$ used for the calculations in this paper. 

The potential, $\phi$ , is the rate of change of the free energy $\textstyle { \mathcal { F } } = \int f d V$ with respect to $C$ , 

$$
\phi = \frac {\delta \mathcal {F}}{\delta C} = \beta \Psi^ {\prime} (C) - \alpha \vec {\nabla} \cdot | \vec {\nabla} C | ^ {q - 2} \vec {\nabla} C. \tag {2.2}
$$

For the CSF method the potential is equal to the surface tension times the field curvature. For the phase-field method the potential is $\beta \Psi ^ { \prime } ( C ) - \alpha \nabla ^ { 2 } C$ . 

Van der Waals [24] hypothesized that equilibrium interface profiles are those that minimize the integral of $f$ . From the calculus of variations, these profiles satisfy $\beta \Psi ^ { \prime } ( C ) -$ $\alpha \vec { \nabla } \cdot | \vec { \nabla } C | ^ { q - 2 } \vec { \nabla } C \equiv \phi = \mathrm { c o n s t a n t }$ . Cahn and Hilliard [5] extended van der Waals’ hypothesis to time-dependent situations by approximating interfacial diffusion fluxes as being proportional to chemical potential gradients. The Cahn–Hilliard equation 

$$
\frac {\partial C}{\partial t} = \kappa \nabla^ {2} \phi = - \kappa \nabla^ {2} (\alpha \nabla^ {2} C - \beta \Psi^ {\prime} (C)) \tag {2.3}
$$

models the creation, evolution, and dissolution of diffusively controlled phase-field interfaces (Bates and Fife [3]). 

The further extension to diffuse-interface fluid-dynamics is discussed by, among others, Antanovskii, Jasnow and Vi˜nals, and Joseph [12]. The derivation of the diffuse-interface fluid-dynamical forcing is fairly simple, especially for compressible flow. The key ideas are (1) convection can change the amount of free energy by either lengthening or thickening/thinning interfaces, (2) there must be a diffuse-interface force exerted by the fluid such that the change in kinetic energy is always opposite to the change in free energy, (3) this must be true for arbitrary interface configurations and (compressible) velocity fields. The rate of change of free energy due to convection is $\begin{array} { r } { \int \phi ( \partial C / \partial t ) | _ { \mathrm { c o n v e c t i o n } } d V = } \end{array}$ $\begin{array} { r } { - \int \phi \sum _ { j } ( \partial u _ { j } C / \partial x _ { j } ) d V = \int \sum _ { j } u _ { j } C ( \partial \phi / \partial x _ { j } ) } \end{array}$ dV . The rate of change of kinetic energy due to surface tension forcing $\vec { F }$ is $\textstyle { \int \sum _ { j } F _ { j } u _ { j } } ~ d V$ . For the two to be equal and opposite for arbitrary $C$ and $\vec { u } .$ , it must be true that $F _ { j } = - C ( \partial \phi / \partial x _ { j } )$ . The argument is essentially the same for incompressible flows. An additional force, the gradient of a potential, must be introduced to enforce the incompressibility constraint. 

The incompressible Navier–Stokes equations with either CSF or phase-field surface tension forcing are 

$$
\vec {\nabla} \cdot \vec {u} = 0 \tag {2.4}
$$

$$
\rho \frac {D u _ {i}}{D t} = \rho \frac {\partial u _ {i}}{\partial t} + \rho \sum_ {j} u _ {j} \frac {\partial u _ {i}}{\partial x _ {j}} = - \vec {\nabla} S + \vec {\nabla} \cdot \vec {\tau} _ {\mu , i} - C \vec {\nabla} \phi + g _ {i} \rho . \tag {2.5}
$$

S enforces the incompressibility condition. $\vec { \tau } _ { \mu , i }$ is the viscous stress tensor and $\vec { g }$ is the gravity vector. The phase-field advection-diffusion equation for $C$ is 

$$
\frac {D C}{D t} = \frac {\partial C}{\partial t} + \sum_ {j} u _ {j} \frac {\partial C}{\partial x _ {j}} = \vec {\nabla} \cdot \kappa (C) \vec {\nabla} \phi = \vec {\nabla} \cdot \kappa (C) \vec {\nabla} (\beta \Psi^ {\prime} (C) - \alpha \nabla^ {2} C). \tag {2.6}
$$

This is the Cahn–Hilliard equation plus advection. $\kappa$ is the diffusion parameter, called the mobility. (The diffusivity in the bulk phases is $\kappa \Psi ^ { \prime \prime } ( C _ { \mathrm { b u l k p h a s e } } ) . )$ The evolution of the fluid’s total energy is found by multiplying (2.5) by uE, (2.6) by $\phi$ , adding, and integrating. The result, neglecting some very small effects due to density diffusion, is 

$$
\frac {\partial \mathcal {E}}{d t} = - \int (\kappa | \vec {\nabla} \phi | ^ {2} + \mu (C) \vec {\nabla} \vec {u} \cdot \vec {\nabla} \vec {u}) d V, \tag {2.7}
$$

where $\mu ( C )$ is the dynamic viscosity. 

In order to isolate interface and surface-tension effects the calculations in this paper will be restricted to a Boussinesq fluid with the two phases having the same viscosity and mobility. Equations (2.4)–(2.6) then simplify to 

$$
\vec {\nabla} \cdot \vec {u} = 0 \tag {2.8}
$$

$$
\rho_ {0} \frac {D u _ {i}}{D t} = - \vec {\nabla} S + \mu \nabla^ {2} u _ {i} - C \vec {\nabla} \phi + g _ {i} \tilde {\rho} (C) \tag {2.9}
$$

$$
\frac {D C}{D t} = \kappa \nabla^ {2} \phi = \kappa \nabla^ {2} (\beta \Psi^ {\prime} (C) - \alpha \nabla^ {2} C). \tag {2.10}
$$

These are the equations that will be solved in this paper. $\rho _ { 0 }$ is the mean density and $\tilde { \rho } ( C )$ is the perturbation from $\rho _ { 0 }$ . 

$- C \vec { \nabla } \phi$ is the continuum surface tension forcing in its potential form. This forcing can be manipulated into a stress form, which for general $q$ is 

$$
\tau_ {j j} = \alpha \sum_ {i \neq j} \left(\frac {\partial C}{\partial x _ {i}}\right) ^ {2} / | \vec {\nabla} C | ^ {2 - q} \tag {2.11}
$$

$$
\tau_ {i j, i \neq j} = - \alpha \left(\frac {\partial C}{\partial x _ {i}} \frac {\partial C}{\partial x _ {j}}\right) / | \vec {\nabla} C | ^ {2 - q}. \tag {2.12}
$$

The principle axes of the tensor are directed in and perpendicular to the tangent plane of the interface. The normal stress perpendicular to the plane is zero and the two tangent normal stresses are equal. When this form is used S becomes the true pressure. The phase-field relationship between $S _ { \mathrm { s t r e s s } }$ and $S _ { \mathrm { p o t } }$ is $\begin{array} { r } { S _ { \mathrm { s t r e s s } } = p = S _ { \mathrm { p o t } } + C \phi - \beta \Psi + \frac { 1 } { 2 } \alpha | \vec { \nabla } C | ^ { 2 } } \end{array}$ . 

Because of incompressibility, the potential form of the surface tension forcing can also be written as $\phi \vec { \nabla } C$ . The new S for this is equal to the old S plus $\phi C$ . The actual motion causing component of the surface tension—as versus the pressure forcing component—can be seen in the vorticity equation. In two dimensions the surface tension vortical forcing is $C _ { y } \phi _ { x } - C _ { x } \phi _ { y }$ . Both $\vec { \nabla } C$ and $\vec { \nabla } \phi$ must simultaneously be nonzero for there to be non-trivial velocity forcing. Asymptotically, the forcing occurs only at interfaces. The component of the potential that actually causes motion is that which varies parallel to the interface. 

The viscous wall boundary condition used for this paper is no-slip. Two boundary conditions are needed for C. The no-flux condition is 

$$
\frac {\partial \phi}{\partial x _ {n}} = 0. \tag {2.13}
$$

$x _ { n }$ denotes the direction normal to the wall. The second boundary condition depends on the interface at the wall being at or near local equilibrium. Postulating that the wall free energy is of the form 

$$
\mathcal {F} _ {w} = \int \gamma g (C) d A, \tag {2.14}
$$

that is, that the wall-fluid interfacial energy is a function only of the fluid composition right against the wall, then the resulting phase-field natural boundary condition, which corresponds to a diffusively controlled local equilibrium at the wall, is 

$$
\alpha \frac {\partial C}{\partial x _ {n}} + \gamma g ^ {\prime} (C) = 0. \tag {2.15}
$$

This condition is analogous to the classical contact angle condition, in which the dynamic contact angle right at the wall (the microscopic contact angle) is taken to be the same as the static equilibrium angle. A more general condition that allows nonequilibrium is 

$$
\frac {D C}{D t} = D _ {w} \left(\alpha \frac {\partial C}{\partial x _ {n}} + \gamma g ^ {\prime} (C)\right). \tag {2.16}
$$

This results in the microscopic contact angle being a function of wall velocity. C approaches the equilibrium condition as $D _ { w }$ , the “wall diffusion,” increases to infinity. All the computations presented in this paper use boundary condition (2.15) with g(C) equal to zero, so both the equilibrium and the dynamic contact angle will be 90◦. 

# 3. PHASE-FIELD SURFACE TENSION AND INTERFACE WIDTH

A calculation’s surface tension and interface width are controlled through 9, α, and $\beta$ . In an isothermal fluid system in equilibrium the surface tension of an interface is equal to the integral of the free energy density through the interface. The equilibrium interface profile is the profile that minimizes $\mathcal { F }$ and it can thus be found from the free energy functional via the calculus of variations. 

The plane interface gives the simplest case. From the stress form of the phase-field equations the surface tension of a phase-field plane interface is 

$$
\sigma = \alpha \int_ {- \infty} ^ {+ \infty} \left(\frac {d C}{d x}\right) ^ {2} d x. \tag {3.1}
$$

The interface profile that minimizes $\mathcal { F }$ obeys 

$$
\alpha \frac {d ^ {2} C}{d x ^ {2}} - \beta \Psi^ {\prime} C) = - \mu = 0. \tag {3.2}
$$

Multiplying by $d C / d x$ and integrating, this becomes 

$$
\frac {\alpha}{2} \left(\frac {d C}{d x}\right) ^ {2} = \beta \Psi (C). \tag {3.3}
$$

The $\Psi ( C )$ used for all the calculations in this paper is $( C + 1 / 2 ) ^ { 2 } ( C - 1 / 2 ) ^ { 2 }$ . This is the simplest non-singular $\Psi ( C )$ that has two equal energy minima. The equilibrium profile of $C$ for this is $\begin{array} { r } { C = \frac { 1 } { 2 } \operatorname { t a n h } ( \zeta ) } \end{array}$ , where $\zeta = \sqrt { \beta / 2 \alpha } x$ . The surface tension stresses $\tau _ { y y }$ and $\tau _ { z z }$ are $\beta \ \mathrm { s e c h } ^ { 4 } ( \zeta ) / 8$ . The surface tension is $\sqrt { \alpha \beta / 1 8 } . \epsilon$ for this interface will be defined in this paper to be the distance from $C = - . 4 5 \mathrm { t o } C = + . 4 5 ( 9 0 \%$ of the variation of C ). This is given by $2 \operatorname { t a n h } ^ { - 1 } ( . 9 ) { \sqrt { 2 \alpha / \beta } } \simeq 4 . 1 6 4 { \sqrt { \alpha / \beta } }$ . This width contains 98.5% of the surface tension stress. In general, for general 9, the surface tension of an interface is proportional to $\sqrt { \alpha \beta }$ while its thickness is proportional to $\sqrt { \alpha / \beta }$ . 

A class of 9 of interest is 9 that have singular behavior at $C _ { \mathrm { b u l k p h a s e } }$ . One example is $\Psi = | C - 1 / 2 | ^ { 3 / 2 } | C + 1 / 2 | ^ { 3 / 2 }$ , for which $\Psi _ { \mathrm { b u l k p h a s e } } ^ { \prime }$ has square root behavior. Another is the double obstacle energy recently used by Oono and Puri [19] and Nochetto [18]. This $\Psi ( C ) \mathrm { i s } - { \textstyle \frac { 1 } { 2 } } ( C + 1 / 2 ) ( C - 1 / 2 )$ for $\begin{array} { r } { - \frac { 1 } { 2 } < C < \frac { 1 } { 2 } } \end{array}$ and +∞ for $\begin{array} { r } { | C | > \frac { 1 } { 2 } } \end{array}$ . These models have some advantages over non-singular models, but they also necessarily raise some difficult numerical issues. This paper will therefore discuss some of their characteristics but no numerical work will be done with them. One advantage of the double obstacle model is that it gives a sharply defined interface. Its equilibrium plane interface profile is $\begin{array} { r } { C = \frac { 1 } { 2 } \sin ( \sqrt { \beta / \alpha } x ) } \end{array}$ for $\textstyle | x | \leq { \frac { \pi } { 2 } } { \sqrt { \alpha / \beta } } ; C = \pm { \frac { 1 } { 2 } }$ otherwise. The interface width is $\pi \sqrt { \alpha / \beta }$ and its surface tension is $\textstyle { \frac { \pi } { 8 } } { \sqrt { \alpha \beta } }$ . 

Interface curvature changes phase-field surface tension. The incurred error in both surface tension and pressure jump is a quadratic function of interface thickness times curvature. The coefficient of this error is small. For example, the error for the double obstacle energy is less than 0.25% for a thickness to Gibbs-radius ratio of $1 / 3$ . For a ratio of 1 the error is about 4%. For $\Psi ( C ) = ( C + 1 / 2 ) ^ { 2 } ( C - 1 / 2 ) ^ { 2 }$ the error is less than 0.13% for thickness times curvature equal to .2. 

# 4. PHASE-FIELD FLOW PHYSICS

This section gives an overview of how solutions to the phase-field Navier–Stokes equations behave as $\epsilon  0$ . The emphasis is on behaviors that are not seen in the original sharp-interface equations and on how these behaviors can, in the limit, be suppressed. It is desired, of course, that the diffuse-interface solutions converge to solutions of the sharpinterface Navier–Stokes equations. This, it will be shown, places constraints on $\kappa ( \epsilon )$ . κ must go to zero along with ², otherwise there is a formally $O ( 1 )$ error due to diffusive transport. But it must approach zero slowly enough so that interface profiles can be maintained against convective distortion. Unlike a sharp interface, a diffuse interface can be subject to thickening/thinning modes. These can create chemical potential boundary layers that can lead to incorrect interface behavior and that can also be a major source of numerical error. These modes are suppressed as $\epsilon \to 0 \mathrm { i f } \kappa ( \epsilon )$ is given the right behavior. 

The chemical potential is the phase-field analogue to surface tension times curvature and as such it is a very important variable. It can sometimes provide the key to understanding a particular physical or numerical issue. $\mathbf { A s } \ \epsilon \to 0$ , for appropriate $\kappa ( \epsilon )$ , interfaces tend more and more to take on their equilibrium profile and the chemical potential tends to definite $O ( 1 )$ values. The variation in interface profile corresponding to this $O ( 1 )$ interfacial chemical potential is $O ( \epsilon )$ and the variation in the interface energy and surface tension is $O ( \epsilon ^ { 2 } )$ . Exceptions to this occur during near-singular events such as interface creation or disappearance or during interface coalescence or at moving contact lines. These, however, are all instances when the Navier–Stokes equations themselves fail and are invalid. 

Solutions of the purely diffusive Cahn–Hilliard equation have potential fields that are generally smooth. The Cahn–Hilliard equation has curvature-dependent solubility (the Gibbs– Thomson effect), which is what makes it useful for modeling nucleation, evaporation, and coarsening. Regions of high curvature are generally also regions of high potential and high solubility—material from these regions fluxes into the surrounding lower-potential medium. The extent of this solubility depends on the model of $\Psi ,$ . The solubility for the 9 used for this paper is $O ( \epsilon )$ . 9 singular at $C _ { \mathrm { b u l k p h a s e } }$ have lower orders of solubility. For $\Psi = | C - 1 / 2 | ^ { 3 / 2 } | C + 1 / 2 | ^ { 3 / 2 }$ the solubility is $O ( \epsilon ^ { 2 } )$ and for the double obstacle energy it is zero. 

When fluid convection is introduced, the chemical potential is no longer necessarily smooth. This has two causes. The first is that convective straining can tend to thicken or thin an interface. This strain is $O ( 1 )$ and will be opposed by an $O ( 1 )$ divergence of the diffusion flux. Since this divergence occurs over an $O ( \epsilon )$ thick interface the strain induced perturbation to the chemical potential is $O ( \epsilon ^ { 2 } / \kappa )$ . The second cause is related to curvature dependent solubility. Oscillation of an interface changes its curvature and thus the local solubility on an $O ( 1 )$ time scale. Alternatively, an interface may be convected through a surrounding inhomogenous fluid. In either case, C becomes out of balance across the interface by $O ( \epsilon )$ (the magnitude of the solubility) and the chemical potential gains a jump of $O ( 1 )$ . The diffusive divergence that can correct this imbalance is $O ( \epsilon )$ . Denoting the flux magnitude as $O ( \gamma )$ , the corresponding boundary layer thickness (the distance from the interface that is put into balance in an $O ( 1 )$ timeframe) is $O ( \gamma / \epsilon )$ . This thickness implies a flux of $O ( \kappa \epsilon / \gamma )$ . Equating this to $O ( \gamma )$ gives a boundary layer thickness of $O ( \sqrt { \kappa / \epsilon } )$ . 

The first boundary layer is important for numerical reasons. It sets limits on the accuracy of $\phi$ interpolation and differentiation for the fourth order accurate method discussed in Section 7. Dissipation due to this type of boundary layer is $O ( \epsilon ^ { 3 } / \kappa )$ . Momentum forcing by it is comparatively negligible. A model problem for this is the steady state straining $u = - x , v = y$ , with the interface parallel to the x-axis. Equation (2.10) becomes 

$$
\kappa \nabla^ {2} \phi = y C _ {y}. \tag {4.1}
$$

Both C and the resulting $\phi$ are antisymmetric. The force on the fluid, which for this discussion is most conveniently written as $\phi C _ { y }$ , therefore integrates to zero. These symmetries are approximately maintained with curved interfaces. The integrated force is then approximately $O ( \epsilon ^ { 4 } / \kappa )$ . Note that this dominance of dissipation over momentum forcing is desirable. The interface maintains its form though diffusion rather than by distorting the velocity field. 

In addition to curvature induced solubility, the phase-field Navier–Stokes equation can exhibit (1) the generation of wall layers, (2) disjoining pressures, and (3) overshoots or undershoots of C past $C _ { \mathrm { b u l k p h a s e } }$ . The wall layers are due to boundary conditions (2.15)– (2.16). They can be largely eliminated by choosing a form for $g ( C )$ so that $g ^ { \prime } ( C _ { \mathrm { b u l k p h a s e } } ) = 0$ . A disjoining pressure builds when two interfaces become very close. It can cause instabilities that hasten coalescence. This is imitative of what really occurs at microscale lengths and has no noticeable negative effect on simulations. For the 9 used in this paper, the overshoots are $O ( \epsilon )$ . the singular 9 that reduce or eliminate curvature induced solubility also reduce or eliminate these overshoots. 

C is O(1) with nth interfacial derivative $O ( \epsilon ^ { - n } )$ . Velocities are smooth across the interfaces. In general, nth interfacial derivatives of the velocities are $O ( \epsilon ^ { 1 - n } )$ . The pressure jumps across interfaces but S, which will be used in all the numerics, is smooth. 

The boundary layer analyses above imply constraints on κ relative to ². The desired asymptotic behavior of the interfacial chemical potential is that it be constant across an interface. The component of the potential induced by straining must therefore be asymptotically small compared to the potential’s smoother components. This requires 

$$
\kappa = O (\epsilon^ {\delta}), \quad \delta <   2. \tag {4.2}
$$

Also, solubility-related boundary layers must be thick compared to the interface. This yields the same constraint. A general physical argument that once again yields this constraint is that as $\epsilon , \kappa \to 0$ the interface should stay closer and closer to its equilibrium profile, so that its tension remains closer and closer to its desired value. For diffusion to dominate convective distortion (4.2) must be true. 

There are also constraints on the minimum value of δ. Assuming 9 is such that results can be $O ( \epsilon )$ accurate (phase-field results $O ( \epsilon )$ different from sharp-interface results for variables of interest), then diffusive fluxes of C across $O ( 1 )$ length scales should be allowed to be no greater than $O ( \epsilon )$ . With the potential having an $O ( 1 )$ variation across an $O ( 1 )$ domain this implies that 

$$
\delta \geq 1. \tag {4.3}
$$

# 5. MODEL AND NUMERICAL CONVERGENCE

Analysis of the accuracy of phase-field computations is complicated by the fact that convergence is governed by three factors, not just mesh spacing but also the interface thickness and the mobility. The mobility affects the thickness and perturbation magnitude of the chemical potential boundary layers. The rate of convergence of a set of calculations is given by a double limit which is a combination of the asymptotics of the approach of the phase-field model to the physics of the exact sharp interface and of the convergence of the numerical methods to the exact solution of the phase-field model. In this limit, mesh size h, interface thickness ², and mobility κ must all be reduced to zero. The interface thickness must be reduced at a slower rate than the mesh size, so as to obtain a more and more accurate estimation of the interfacial forces. The relative local truncation error of numerical differentiation of C is proportional to $( h / \epsilon ) ^ { n }$ , where $h / \epsilon$ is the mesh size scaled by interface thickness, and n is the order of accuracy of the numerical approximation. Truncation error order is not the same as solution error order; the solution error order can be the same, worse, or, since the truncation error is restricted to interface regions, even better. But assuming for the moment it is the same then equating this truncation error to the $O ( \epsilon )$ error incurred by this paper’s phase-field model indicates an optimal convergence rate of $O ( h ^ { n / n + 1 } )$ . The corresponding optimal relationship between interface thickness and mesh size is then ² ∝ hn/n+1. $\epsilon \propto h ^ { n / n + 1 }$ 

For more rapid convergence a more accurate phase-field model must be used. $O ( \epsilon ^ { 2 } )$ accurate results are sometimes possible with the model used in this paper. An example is linear plane waves, which will be discussed in Section 10. Then $O ( h ^ { 2 n / n + 2 } )$ convergence can be achieved. 

# 6. SECOND ORDER CENTRAL DIFFERENCE METHODS AND ENERGY CONSERVATION

This section discusses some simple central-differenced staggered-grid methods and shows how they conserve energy. The equations being discretized are (2.8)–(2.10) with no-slip and $\partial \phi / \partial x _ { n } = \partial C / \partial x _ { n } = 0$ at walls. The discussion here and throughout the paper will be for a uniform, square, Cartesian, or axisymmetric grid. The discrete S, C, and $\phi$ are located at cell centers and the velocities are at cell faces. The discretizations of the viscous, convective, and gravitational terms in the Navier–Stokes equations are all made using standard second order centered differences. Standard 5-point discrete Laplacians are used to calculate the chemical potential and chemical potential driven diffusion. The chemical potential forcing of the momentums is discretized as 

$$
F _ {x, i + 1 / 2, j} = - \frac {h}{2} (C _ {i, j} + C _ {i + 1, j}) (\phi (i + 1, j) - \phi (i, j)) \tag {6.1a}
$$

$$
F _ {y, i, j + 1 / 2} = - \frac {h}{2} (C _ {i, j} + C _ {i, j + 1}) (\phi (i, j + 1) - \phi (i, j)). \tag {6.1b}
$$

Convective fluxes across cell faces are approximated by 

$$
\frac {h}{2} u _ {i + 1 / 2, j} (C _ {i, j} + C _ {i + 1, j}), \quad \frac {h}{2} v _ {i, j + 1 / 2} (C _ {i, j} + C _ {i, j + 1}). \tag {6.2}
$$

Both (6.1) and (6.2) use second order central differencing. The discretization has an ${ \cal O } ( h ^ { 2 } / \epsilon ^ { 2 } )$ relative truncation error in interfaces and an $O ( h ^ { 2 } )$ truncation error elsewhere. From Section 5, its optimum convergence, for $O ( \epsilon )$ accurate phase-field models, is hypothesized to be $O ( h ^ { 2 / 3 } )$ . 

Ignoring wall energies, the free energy of the discretized system is 

$$
\mathcal {F} = \frac {1}{2} \alpha \sum_ {i, j} (C _ {i + 1, j} - C _ {i, j}) ^ {2} + \frac {1}{2} \alpha \sum_ {i, j} (C _ {i, j + 1} - C _ {i, j}) ^ {2} + h ^ {2} \beta \sum_ {i, j} \Psi_ {i, j}. \tag {6.3}
$$

The kinetic energy is 

$$
\mathcal {K} = \frac {1}{2} h ^ {2} \rho_ {0} \sum_ {i, j} u _ {i + \frac {1}{2}, j} ^ {2} + \frac {1}{2} h ^ {2} \rho_ {0} \sum_ {i, j} v _ {i, j + \frac {1}{2}} ^ {2}. \tag {6.4}
$$

The discretization is energy conserving, as follows. The rate of change in free energy due to convection is 

$$
\begin{array}{l} \frac {d \mathcal {F}}{d t} = \frac {h}{2} \sum_ {i, j} \phi_ {i, j} u _ {i - \frac {1}{2}, j} (C _ {i - 1, j} + C _ {i, j}) - \frac {h}{2} \sum_ {i, j} \phi_ {i, j} u _ {i + \frac {1}{2}, j} (C _ {i, j} + C _ {i + 1, j}) \\ + \frac {h}{2} \sum_ {i, j} \phi_ {i, j} v _ {i, j - \frac {1}{2}} (C _ {i, j - 1} + C _ {i, j}) - \frac {h}{2} \sum_ {i, j} \phi_ {i, j} v _ {i, j + \frac {1}{2}} (C _ {i, j} + C _ {i, j + 1}). \tag {6.5} \\ \end{array}
$$

This is found by multiplying the discretized advection-diffusion equations for $C _ { i , j }$ by $h ^ { 2 } \phi _ { i , j }$ and summing. The rate of change in kinetic energy due to surface tension forcing is 

$$
\begin{array}{l} \frac {d \mathcal {K}}{d t} = - \frac {h}{2} \sum_ {i, j} u _ {i + \frac {1}{2}, j} (C _ {i + 1, j} + C _ {i, j}) (\phi_ {i + 1, j} - \phi_ {i, j}) \\ - \frac {h}{2} \sum_ {i, j} v _ {i, j + \frac {1}{2}} (C _ {i, j + 1} + C _ {i, j}) (\phi_ {i + 1, j} - \phi_ {i, j}). \tag {6.6} \\ \end{array}
$$

Equations (6.5) and (6.6) sum to zero, as can be found by reindexing the first $( i  i + 1 )$ and third $( j  j + 1 )$ sums of (6.5). 

The method is very easy to implement. Its major problem is that it requires fairly wide interfaces. As will be discussed in Sections 8 and 9, for acceptable results α and $\beta$ must be such that the interfaces are at least $3 \textstyle { \frac { 1 } { 2 } }$ cells wide. For $\Psi = ( C - 1 / 2 ) ^ { 2 } ( C + 1 / 2 ) ^ { 2 }$ , and using the definition of interface width given in Section 3, this constrains $\alpha / \beta$ to be $\geq ( 3 . 5 / 4 . 1 6 4 ) ^ { 2 } h ^ { 2 } \simeq . 7 1 h ^ { 2 }$ . 

One way to reduce interface width is to use a finer grid for the color function and potential. This is acceptable costwise because, compared to the velocities and pressure field, C and $\phi$ are fairly inexpensive to calculate. The method works well when the potential form of the momentum forcing is used. It is then easy to manage energy transfers between the coarse velocity grid and the fine C grid so that energy is conserved. 

For a uniform grid, the simplest approach is to divide each pressure cell into $n \times n \left\{ C , \phi \right\}$ cells. u is approximated in each pressure cell as varying linearly in the x direction and as constant in the y direction, and vice versa for v. This yields the discrete velocities needed for solving (2.10). Equation (2.10) is discretized on the fine grid, using all centered differences. Numbering the $\{ C , \phi \}$ cells in each $p _ { i , j }$ cell from $k = n ( i - 1 ) + 1$ to ni and $l = n ( j - 1 )$ to $n j$ , the effect of convection on the fluid free energy can be written as 

$$
\begin{array}{l} \frac {d \mathcal {F}}{d t} = \frac {h}{2} \sum_ {i, j} u _ {i + \frac {1}{2}, j} \sum_ {l = n (j - 1) + 1} ^ {l = n j} \sum_ {k = n (i - 1)} ^ {k = n (i + 1)} \frac {n - | k - n i |}{n} (\phi_ {k + 1, l} - \phi_ {k, l}) (C _ {k, l} + C _ {k + 1, l}) \\ + \frac {h}{2} \sum_ {i, j} v _ {i, j + \frac {1}{2}} \sum_ {k = n (i - 1) + 1} ^ {k = n i} \sum_ {l = n (j - 1)} ^ {l = n (j + 1)} \frac {n - | l - n j |}{n} (\phi_ {k, l + 1} - \phi_ {k, l}) (C _ {k, l} + C _ {k, l + 1}). \tag {6.7} \\ \end{array}
$$

h is the fine mesh size. The second order accurate u forcing that conserves energy is 

$$
F _ {x, i + \frac {1}{2}, j} = - \frac {h}{2} \sum_ {l = n (j - 1) + 1} ^ {l = n j} \sum_ {k = n (i - 1)} ^ {k = n (i + 1)} \frac {n - | k - n i |}{n} (\phi_ {k + 1, l} - \phi_ {k, l}) (C _ {k, l} + C _ {k + 1, l}) \tag {6.8}
$$

and similarly for v. In general, for energy conservation the x-direction force on each discrete $u _ { n }$ and the y direction force on each discrete $v _ { n }$ must obey 

$$
F _ {x, n} = - \frac {\delta}{\delta u _ {n}} \frac {d \mathcal {F}}{d t}; \quad F _ {y, n} = - \frac {\delta}{\delta v _ {n}} \frac {d \mathcal {F}}{d t}. \tag {6.9}
$$

Relations (6.9) can be stated in words as: The rate in change of kinetic energy of a particular velocity component is opposite to the rate of change of free energy caused by that velocity component. This holds for arbitrary grid systems and interpolations. 

# 7. FOURTH ORDER COMPACT METHODS

Section 4 discussed how the chemical potential forms two types of interfacial boundary layers. One is $O ( \epsilon )$ in thickness but with a perturbation to $\phi$ of only $O ( \epsilon ^ { 2 } / \kappa )$ . The second has a thickness of $O ( \sqrt { \kappa / \epsilon } )$ with a perturbation of $O ( 1 )$ . Away from coalescence, interface appearance and disappearance, and moving contact line near-singularities, the chemical potential is otherwise generally $O ( 1 )$ and smooth. nth derivatives of $\phi$ in interfacial regions are $O ( \epsilon ^ { 2 - n } / \kappa )$ or $O ( \epsilon ^ { n / 2 } / \kappa ^ { n / 2 } )$ . In either case nth derivatives of $\phi$ are always of smaller magnitude than nth derivatives of C provided that $\kappa = O ( \epsilon ^ { \delta } ) , \delta < 2$ . 

There are various ways to take advantage of both $\phi \ ' { } s { } O ( 1 )$ magnitude (the following discussion will not consider near-singularities) and its comparative smoothness. This section discusses a fairly easy way to construct a compact ${ \cal O } ( h ^ { 4 } / \epsilon ^ { 4 } )$ finite volume discretization on a square Cartesian or axisymmetric grid. It also very briefly discusses a formulation for more general grids. As discussed in Section 5, the optimum choice of $\epsilon$ in terms of h, for fourth order discretizations for $O ( \epsilon )$ accurate phase-field models, appears to be $\epsilon \propto h ^ { 4 / 5 }$ . This gives a convergence rate of $O ( h ^ { 4 / 5 } )$ ). 

To calculate fluxes, finite volume methods need to find interfacial values of C from given cell averages of $C$ . The second order discretizations discussed in Section $^ 6$ are equivalent to finite volume methods that take $\bar { \bar { C } } _ { i , j }$ , the cell average, to be identical to $C _ { i , j }$ , the cell midpoint value. This approximation has an error of ${ \cal O } ( h ^ { 2 } / \epsilon ^ { 2 } )$ . The linear interpolation to then find C interface values from C midcell values has the same order of error. 

However, the cell midpoint value can easily be found more accurately. With error terms included, 

$$
C _ {i, j} = \bar {\bar {C}} _ {i, j} - \frac {h ^ {2}}{2 4} (\nabla^ {2} C) _ {i, j} + O (h ^ {4} / \epsilon^ {4}) + O (h ^ {6} / \epsilon^ {6}). \tag {7.1}
$$

The equation for the chemical potential is 

$$
\alpha \nabla^ {2} C - \beta \Psi^ {\prime} (C) = - \phi = O (1) \tag {7.2}
$$

from which 

$$
(\nabla^ {2} C) _ {i, j} = \frac {\beta}{\alpha} \Psi^ {\prime} (C _ {i, j}) - \phi_ {i, j} / \alpha . \tag {7.3}
$$

The first term on the right hand side is $O ( \epsilon ^ { - 2 } )$ , the second is $O ( \epsilon ^ { - 1 } )$ . Substituting into (7.1) gives 

$$
C _ {i, j} + \frac {h ^ {2}}{2 4} \frac {\beta}{\alpha} \Psi^ {\prime} (C _ {i, j}) = \bar {\bar {C}} _ {i, j} + O (h ^ {4} / \epsilon^ {4}) + O (h ^ {2} / \epsilon) + O (h ^ {6} / \epsilon^ {6}) \tag {7.4}
$$

a local nonliner equation for $C _ { i , j }$ . For $\epsilon = O ( h ^ { 4 / 5 } )$ the dominant error term is ${ \cal O } ( h ^ { 4 } / \epsilon ^ { 4 } )$ . $O ( h ^ { 2 } / \epsilon )$ then is ${ \cal O } ( h ^ { 6 } / \epsilon ^ { 6 } )$ so there is no reason to evaluate the $\phi _ { i , j } / \alpha$ term in (7.3). 

C must now be found on the cell interfaces. This can be done indirectly but compactly by approximating $\phi$ between cell midpoints to be bilinear and then solving the potential equation (7.2) using high order differencings for $C$ . This can be viewed as generating a system of equations subject to constraints. The constraints are the previously calculated midpoint values of C. The unknowns are the midpoint values of $\phi$ and the values of C needed on the cell interfaces. There is one unknown discrete $\phi$ per cell but there can be an arbitrary number of unknown Cs. In strained interfaces, from Section 4, a linear approximation to $\phi$ between cells incurs an $O ( h ^ { 2 } / \kappa )$ error. For the interpolation to be everywhere $O ( h ^ { 4 / 5 } )$ accurate, κ must therefore be at least equal in magnitude to $O ( h ^ { 6 / 5 } ) = O ( \epsilon ^ { 3 / 2 } )$ . 

The simplest approach is to find C at cell corners and cell interface midpoints. Once found, these then allow an ${ \cal O } ( h ^ { 4 } / \epsilon ^ { 4 } )$ accurate Simpson’s rule integration to find interface fluxes. Simpson’s rule gives an average interface C according to 

$$
\bar {C} _ {i + 1 / 2, j} = \frac {1}{6} (C _ {i + 1 / 2, j - 1 / 2} + 4 C _ {i + 1 / 2, j} + C _ {i + 1 / 2, j + 1 / 2}) \tag {7.5a}
$$

$$
\bar {C} _ {i, j + 1 / 2} = \frac {1}{6} (C _ {i - 1 / 2, j + 1 / 2} + 4 C _ {i, j + 1 / 2} + C _ {i + 1 / 2, j + 1 / 2}). \tag {7.5b}
$$

u and v convective fluxes are then 

$$
h u _ {i + 1 / 2, j} \bar {C} _ {i + 1 / 2, j}, \quad h v _ {i, j + 1 / 2} \bar {C} _ {i, j + 1 / 2}. \tag {7.6}
$$

The present method, unlike the second-order methods, has no explicit discrete energy. The approximate rate of change of free energy by convection is given by 

$$
\begin{array}{l} h \sum_ {i, j} \phi_ {i, j} u _ {i - \frac {1}{2}, j} \bar {C} _ {i - 1 / 2, j} - h \sum_ {i, j} \phi_ {i, j} u _ {i + \frac {1}{2}, j} \bar {C} _ {i + 1 / 2, j} \\ + h \sum_ {i, j} \phi_ {i, j} v _ {i, j - \frac {1}{2}} \bar {C} _ {i, j - 1 / 2} - h \sum_ {i, j} \phi_ {i, j} v _ {i, j + \frac {1}{2}} \bar {C} _ {i, j + 1 / 2}. \tag {7.7} \\ \end{array}
$$

The surface tension forcing corresponding to this is 

$$
F _ {x, i + 1 / 2, j} = - h \bar {C} _ {i + 1 / 2, j} (\phi (i + 1, j) - \phi (i, j)) \tag {7.8a}
$$

$$
F _ {y, i, j + 1 / 2} = - h \bar {C} _ {i, j + 1 / 2} (\phi (i, j + 1) - \phi (i, j)). \tag {7.8b}
$$

Its relative error is the maximum of ${ \cal O } ( h ^ { 4 } / \epsilon ^ { 4 } )$ and $O ( h ^ { 2 } / \kappa )$ . The first is the error in the approximation of $C _ { i }$ the second is the error due to the second-order-accurate differentiation of $\phi$ . Diffusive fluxes of $C$ can be approximated by 

$$
f _ {x, i + 1 / 2, j} = - \kappa (\phi (i + 1, j) - \phi (i, j)) \tag {7.9a}
$$

$$
f _ {x, i + 1 / 2, j} = - \kappa (\phi (i + 1, j) - \phi (i, j)). \tag {7.9b}
$$

It remains to discuss how to discretize the chemical potential equation. This is solved for C at the cell corners and cell interface midpoints. At cell corners $\phi$ is approximated as being the average of the $\phi _ { i , j }$ in the four surrounding cells; at interface midpoints it is the average taken from the two cell neighbors. The discretization must be at least ${ \cal O } ( h ^ { 4 } / \epsilon ^ { 4 } )$ accurate. A good method that takes advantage of interface asymptotics is to use a variation of the compact mehrstellungen scheme. Expressed in stencil form, this, together with leading order error terms, is 

$$
\frac {\alpha}{6 h ^ {2}} \left[ \begin{array}{c c c} 1 & 4 & 1 \\ 4 & - 2 0 & 4 \\ 1 & 4 & 1 \end{array} \right] C - \frac {\beta}{1 2} \left[ \begin{array}{c c c} 0 & 1 & 0 \\ 1 & 8 & 1 \\ 0 & 1 & 0 \end{array} \right] \Psi^ {\prime} (C) = - \left[ \begin{array}{c c c} 0 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{array} \right] \phi + \frac {h ^ {2}}{1 2} \nabla^ {2} \phi + O (h ^ {4} / \epsilon^ {5}). \tag {7.10}
$$

Because $\nabla ^ { 2 } \phi$ is $O ( 1 / \kappa )$ the dominant relative truncation error (for C found from $\phi )$ is the maximum of ${ \cal O } ( h ^ { 4 } / \epsilon ^ { 4 } )$ and $O ( \epsilon h ^ { 2 } / \kappa )$ . The equation is solved iteratively. Because of the local constraints it converges very quickly, usually in 4 to 6 iterations, and at a rate that appears to be independent of the number of mesh points. 

The above discretizations can be fairly easily adapted for axisymmetric coordinates. Equation (7.4) remains valid. Its $O ( h ^ { 2 } / \epsilon )$ error now includes the term $- ( h ^ { 2 } / 2 4 ) ( 1 / r ) ( \partial C /$ $\partial r )$ . The transformations needed for the various Simpson rule integrations are obvious. The chemical potential equation becomes 

$$
\frac {\alpha}{6 h ^ {2}} \left[ \begin{array}{c c c} M _ {-} & 4 & M _ {+} \\ 4 M _ {-} & - 2 0 & 4 M _ {+} \\ M _ {-} & 4 & M _ {+} \end{array} \right] C - \frac {\beta}{1 2} \left[ \begin{array}{c c c} 0 & 1 & 0 \\ M _ {-} & 8 & M _ {+} \\ 0 & 1 & 0 \end{array} \right] \Psi^ {\prime} (C) = - \phi + \frac {h ^ {2}}{1 2} \nabla^ {2} \phi + O (h ^ {4} / \epsilon^ {5}). \tag {7.11}
$$

$M _ { \pm } = ( r _ { \pm } + r _ { 0 } ) / 2 r _ { 0 }$ , where $r _ { 0 }$ is $r$ at the stencil’s central point, and $r _ { - }$ and $r _ { + }$ are the r at the stencil’s inner and outer points. Errors now also include $O ( h ^ { 2 } / \epsilon )$ terms that contain the first and second radial derivatives of $C$ . 

With more general grids it becomes hard to apply (7.3). A methodology applicable to a general grid of K finite volumes is to use the known $\bar { \bar { C } } _ { k }$ directly as constraints. Each cell is assigned a midpoint or centroid and a $\phi _ { k }$ located at it. $\phi$ is then linearly interpolated between these points. A fine mesh of L points is used to solve the chemical potential equation so as to find C on the interfaces. The unknowns are the K $\phi _ { k }$ plus the $L C _ { l }$ . The L chemical potential equations are supplemented by the K integral constraints. The constraints are expressed by high-order numerical integrations over the $C _ { l }$ . 

This approach has been implemented for the square mesh by replacing (7.4) with the constraint equation 

$$
\begin{array}{l} \frac {1}{3 6} (C _ {i - 1 / 2, j - 1 / 2} + C _ {i - 1 / 2, j + 1 / 2} + C _ {i + 1 / 2, j - 1 / 2} + C _ {i + 1 / 2, j + 1 / 2}) \\ + \frac {1}{9} (C _ {i - 1 / 2, j} + C _ {i, j + 1 / 2} + C _ {i, j - 1 / 2} + C _ {i + 1 / 2, j}) + \frac {4}{9} C _ {i, j} = \bar {\bar {C}} _ {i, j}. \tag {7.12} \\ \end{array}
$$

This is again Simpson’s rule in two dimensions. Equations (7.5), (7.6), (7.8), (7.9) are unchanged. Using Eq. (7.12) instead of (7.4) results in a somewhat slower iterative solution of the chemical potential equation. There is no significant difference in solution results. 

To repeat a major point contained in the preceding: for optimum convergence κ has to be $O ( \epsilon ^ { \delta } ) , \delta \leq 3 / 2$ . Assuming the phase-field model is $O ( \epsilon )$ accurate, the error of the method is the maximum of $O ( \epsilon ) , O ( h ^ { 4 } / \epsilon ^ { 4 } )$ , and $O ( h ^ { 2 } / \kappa )$ . 

# 8. GRID ROUGHNESS AND ORIENTATION EFFECTS

The movement of a drop over a solid can be noticeably affected by the chemical inhomogeneities and roughness of the solid surface. For example, a small drop moving down a window pane tends to move intermittently and erratically. In numerical simulations the grid imposes the equivalent of a spatial roughness. Interfaces moving through a grid exhibit small structural and energy oscillations as they move from being cell centered to cell-interface centered and then back to being cell centered. This is especially manifest with very narrow interfaces. 


TABLE I Interface Energy Using Second Order Numerics


<table><tr><td>Interface cellwidth ε/h</td><td>Cell centered energy error ×100</td><td>Cell-interface centered energy error ×100</td><td>Energy difference ×100</td></tr><tr><td>3.0</td><td>-2.38954</td><td>-5.37324</td><td>2.98370</td></tr><tr><td>3.5</td><td>-2.19554</td><td>-3.16595</td><td>0.97041</td></tr><tr><td>4.0</td><td>-1.82815</td><td>-2.12010</td><td>0.29195</td></tr><tr><td>4.5</td><td>-1.48281</td><td>-1.56581</td><td>0.08300</td></tr><tr><td>5.0</td><td>-1.20554</td><td>-1.22811</td><td>0.02257</td></tr><tr><td>6.0</td><td>-0.82987</td><td>-0.83137</td><td>0.00150</td></tr><tr><td>7.0</td><td>-0.60436</td><td>-0.60445</td><td>0.00009</td></tr><tr><td>8.0</td><td>-0.45998</td><td>-0.45998</td><td>0.00000</td></tr><tr><td>10.0</td><td>-0.29238</td><td>-0.29238</td><td>0.00000</td></tr></table>

Since the numerical methods used in this analysis are energy conserving or almost energy conserving, a grid-dependent oscillation in interface energy means there must also be an oscillation in kinetic energy. If the kinetic energy at its maximum is less than the difference between the interface’s maximum free energy and minimum then the interface cannot move through the grid. With very narrow interfaces this can occur at fairly high velocities. 

A simple way to estimate the energy roughness of a grid is to calculate the difference between the energy of a static one-dimensional interface when cell centered and when cellinterface centered. To find this, (3.2) is solved for $x > 0$ with C antisymmetric about $x = 0$ . For discrete $C _ { n } , n = 0 , 1 , 2 , \ldots$ , this antisymmetry condition becomes $C _ { 0 } = 0$ for the cell centered case and $C _ { 0 } = - C _ { 1 }$ for the cell-interface centered case. 

Table I gives results for standard second-order differencing. The percentage relative error of the numerical surface tension, $1 0 0 \times ( \sigma _ { \mathrm { e x a c t } } - \sigma _ { \mathrm { n u m } } ) / \sigma _ { \mathrm { e x a c t } } ,$ , is given for both the cell centered and cell-interface centered cases and results are given as a function of the interface cellwidth $\epsilon / h$ . The surface tension error for both cases is ${ \cal O } ( h ^ { 2 } / \epsilon ^ { 2 } )$ but the difference between the two cases decreases exponentially. The difference is acceptable for rough calculations beginning at about $\epsilon / h = 3 . 5$ (a relative difference of 1%), of little effect (difference of .2%) at $\epsilon / h = 4 . 2$ , and negligible (difference of .02%) at $\epsilon / h = 5$ . 

Table II gives results for the fourth order mehrstellungen approximation. The numerical surface tension is fourth order accurate. As discussed in Section 7, the mehrstellungen discretization doesn’t have an explicit discrete energy. For the special case of equilibrium plane interfaces, however, it is possible to calculate the discretization’s energy to high accuracy. From Eqs. (3.1), (3.3) the energy of a plane interface in equilibrium is given by $\sigma = 2 \beta \int _ { - \infty } ^ { + \infty }$ σ = 2β R +∞−∞ 9 dx. A trapezoidal numerical integration of this has exponential accuracy, so errors due to the discretization can also be found to exponential accuracy. The results show that the error in the energy decays like $h ^ { 4 } / \epsilon ^ { 4 }$ and, like the second order approximation, that the interfacial energy roughness decays exponentially. The roughness is acceptable down to about $\epsilon / h = 1 . 6$ and negligible at $\epsilon / h \simeq 2 . 3$ . 


TABLE II Interface Energy Using Fourth Order Numerics


<table><tr><td>Interface cellwidth ε/h</td><td>Cell centered energy error ×100</td><td>Cell-interface centered energy error ×100</td><td>Energy difference ×100</td></tr><tr><td>1.4</td><td>0.38126</td><td>-2.05892</td><td>2.44018</td></tr><tr><td>1.6</td><td>0.01335</td><td>-0.91346</td><td>0.92681</td></tr><tr><td>1.8</td><td>-0.09191</td><td>-0.43001</td><td>0.33810</td></tr><tr><td>2.0</td><td>-0.10300</td><td>-0.22259</td><td>0.11959</td></tr><tr><td>2.2</td><td>-0.08694</td><td>-0.12817</td><td>0.04123</td></tr><tr><td>2.5</td><td>-0.05849</td><td>-0.06651</td><td>0.00802</td></tr><tr><td>3.0</td><td>-0.02903</td><td>-0.02951</td><td>0.00048</td></tr><tr><td>3.5</td><td>-0.01553</td><td>-0.01555</td><td>0.00002</td></tr><tr><td>4.0</td><td>-0.00902</td><td>-0.00902</td><td>0.00000</td></tr><tr><td>5.0</td><td>-0.00365</td><td>-0.00365</td><td>0.00000</td></tr></table>

Another important type of error stems from grid/interface orientation. Interfaces at angles to the grid are better resolved and so more accurately approximated. A $4 5 ^ { \circ }$ orientation gives the best resolution. For the second order method, the energy error for an interface with thickness $\epsilon$ at a $4 5 ^ { \circ }$ orientation is exactly equal to the energy error of a $9 0 ^ { \circ }$ oriented interface with thickness $\sqrt { 2 } \epsilon$ . Both the surface tension error and energy roughness of a $4 5 ^ { \circ }$ oriented interface can therefore be found from Table I by substituting $\sqrt { 2 } \epsilon / h$ for $\epsilon / h$ . From the table, the numerical surface tension of an interface is greater at the $4 5 ^ { \circ }$ orientation than it is when parallel. In general, the variation in numerical surface tension with change in angle is $O ( h ^ { n } / \epsilon ^ { n } )$ . This variation can affect capillary vibrations but it does not significantly affect transport. Given the much greater energy effects of changing interface length, and given volumetric constraints, energy reduction by rotation is unlikely to play any role in forcing interface evolution. 

To summarize, at small interface widths and at low velocities interfaces can become trapped at energy minima in the grid. Just above the trapping threshold the kinetic energy, free energy, and velocity can all exhibit large grid-related oscillations. Fortunately, these effects fade exponentially with increasing interface width. Grid anisotropy effects are $O ( h ^ { n } / \epsilon ^ { n } )$ . 

# 9. ONE-DIMENSIONAL CONVECTION

One-dimensional convection through a grid provides one of the few discrete numerical systems that can be considered analytically. In continuum one-dimensional convection an interface remains at its equilibrium profile, the chemical potential is identically zero, and there is no dissipation. In convection through a grid the interface becomes distorted from its equilibrium profile and an interfacial chemical potential boundary layer builds up. For a given discretization the magnitude of $\phi$ in this layer is a function of the three parameters h, ², and κ . 

If the numerical method is energy conserving, this boundary layer slows the fluid down. Any distortion caused by the numerical convection necessarily increases the interface’s free energy. The interfacial force $- C d \phi / d x$ compensates for this by exerting an opposing force on the fluid that reduces the kinetic energy. The chemical potential gradients that build up act to dissipate the excess free energy and to restore the interface to equilibrium. When these various effects are in approximate balance and the interface is traveling through the grid with a quasi-periodic profile and energy the time-averaged decrease in kinetic energy caused by $C d \phi / d x$ must be equal to the time-averaged dissipation of free energy. 

The potential is forced by the error in the discretization of the interface convection. For the second order method the truncation error is $O ( U h ^ { 2 } / \epsilon ^ { 3 } )$ . The resulting chemical potential should be this order multiplied by $\epsilon ^ { 2 } / \kappa$ , or $O ( U h ^ { 2 } / \kappa \epsilon )$ . The free energy dissipation, from Eq. (2.7) applied over the $O ( \epsilon )$ interfacial region where $\phi$ is large, should be $O ( U ^ { 2 } h ^ { 4 } / \kappa \epsilon ^ { 3 } )$ and the integral of the interfacial fluid forcing should necessarily be the same. These estimates have been checked by carrying out calculations of interface convection with increasingly refined grids. The convection was given a fixed velocity so a steady-periodic state could be achieved. The integral of $- C d \phi / d x$ was calculated at each time step but this force was not then applied against the flow. Table III shows maximum chemical potential, the free energy dissipation, the interfacial force exerted on the fluid, and the relative surface tension error for the particular case of U = 1 cm/s and $\sigma = 3 0$ dynes/cm. $\kappa = O ( \epsilon )$ , with $\kappa$ for the narrowest cellwidth being $6 \times 1 0 ^ { - 6 } \mathrm { c m s ^ { 5 } / e r g - s }$ . Interface width at the narrowest cellwidth is 0.04 cm, so in one second this interface traverses 25 times its width and 100 mesh cells. Dissipation per unit interface area is given in $\mathrm { e r g s } / \mathrm { s } { \mathrm { - } } \mathrm { c m } ^ { 2 }$ and flow resistance per area in dynes/cm2. Since $U = 1$ the flow resistance and the rate of decrease in kinetic energy are numerically the same. 


TABLE III Instantaneous Maximum Chemical Potential, Dissipation Error and Relative Interface Energy Error as a Function of Grid Resolution, for a Convected Interface Discretized Using Second Order Methods


<table><tr><td>Interface cellwidth ε/h</td><td>|φ|max(t)</td><td>Dissipation</td><td>Flow resistance</td><td>Energy error ×100</td></tr><tr><td>4.00</td><td>61.5/86.2</td><td>2.8/11.8</td><td>-14.2/33.2</td><td>-2.07/-1.77</td></tr><tr><td>5.04</td><td>40.5/43.5</td><td>1.93/2.61</td><td>-1.5/6.0</td><td>-1.19/-1.17</td></tr><tr><td>6.35</td><td>26.17/26.32</td><td>0.83/0.85</td><td>0.63/1.05</td><td>-.736/-.735</td></tr><tr><td>8.00</td><td>16.66/16.67</td><td>.323</td><td>.319/.327</td><td>-.459</td></tr><tr><td>10.08</td><td>10.56</td><td>.12572</td><td>.12570/.12574</td><td>-.287548</td></tr></table>

The increase in interface cellwidth from one line of the table to the next lower is by a factor of $2 ^ { 1 / 3 }$ . The predicted decrease in $| \phi | _ { \mathrm { m a x } }$ is by $2 ^ { 2 / 3 }$ (1.587) and in the dissipation and flow resistance by $2 ^ { 4 / 3 }$ (2.520). This behavior is seen for cellwidths of 6 and greater. For example, between cellwidths 6.35 and 8.00, $| \phi | _ { \mathrm { m a x } }$ decreases by a factor of 1.579 and the dissipation by 2.632. $| \phi | _ { \mathrm { m a x } }$ and the other variables are all functions of the interface position relative to the grid, hence they are (quasi) periodic. The table gives the variables’ minima and maxima. The range of variation decreases steadily and rapidly as cellwidth increases. This decrease is related to the exponential decay in the grid energy roughness. Though it is not too evident from the results, the deviation of the surface tension from its averaged static value should be approximately equal to the second variation of the free energy integral over the interface, or $( O ( \epsilon \beta ( C - C _ { \mathrm { e q } } ) ^ { 2 } ) = O ( \epsilon ^ { 2 } | \phi | ^ { 2 } ) = O ( U ^ { 2 } h ^ { 4 } / \kappa ^ { 2 } )$ . For the present case this is ${ \cal O } ( h ^ { 8 / 3 } )$ . The energy error, as it must be, is always such that the moving interface has more energy than the static equilibrium interface (the energy errors in Table III are less negative than the static energy errors in Table I). 

The fourth order mehrstellungen method has no explicit free energy and thus no exact energy conservation. Free energy dissipation is not directly measurable. The measurable quantity of significance that relates to dissipation is the interfacial flow resistance. The interfacial convective truncation error for the method is ${ \cal O } ( U h ^ { 4 } / \epsilon ^ { 5 } )$ . The indicated maximum chemical potential is $O ( U h ^ { 4 } / \kappa \epsilon ^ { 3 } )$ . The chemical potential is concentrated about the interface. There is, so far, no theory for the magnitude of the flow resistance; its behavior will be found from numerical experiment. 


TABLE IV Maximum Chemical Potential as a Function of Grid Resolution, for a Convected Interface Discretized Using Fourth Order Methods


<table><tr><td>Interface cellwidth ε/h</td><td>max{|φ|max(t)} κ = O(ε)</td><td>max{|φ|max(t)} κ = O(ε3/2)</td></tr><tr><td>1.636</td><td>188.2</td><td>188.2</td></tr><tr><td>1.879</td><td>103.2</td><td>107.9</td></tr><tr><td>2.158</td><td>45.6</td><td>49.6</td></tr><tr><td>2.479</td><td>15.7</td><td>18.2</td></tr><tr><td>2.848</td><td>4.32</td><td>5.90</td></tr><tr><td>3.272</td><td>1.05</td><td>2.32</td></tr><tr><td>3.758</td><td>0.301</td><td>1.29</td></tr></table>

As a practical matter, the convergence of $\phi$ and the flow resistance to zero is much faster than asymptotic. Very good results are obtained at interface thicknesses of only two to three cellwidths. At this thickness, convergence appears to be almost exponential. The asymptotic regime is entered beginning at about a cellwidth of 5. Table IV gives max $\{ | \phi | _ { \operatorname* { m a x } } ( t ) \}$ as a function of cellwidth in the narrow interface, non-asymptotic regime. Two sequences are given, one with $\kappa = O ( \epsilon )$ and one with $\kappa = O ( \epsilon ^ { 3 / 2 } )$ . For both, κ at the narrowest width is $5 \times 1 0 ^ { - 6 }$ $\mathrm { c m s } ^ { 5 } / \mathrm { e r g - s } .$ U and $\sigma$ are the same as for the second-order calculations. The increase in interface cellwidth from one line of the table to the next lower is by $2 ^ { 1 / 5 }$ . The expected rate of decrease in max $\{ | \phi | _ { \operatorname* { m a x } } ( t ) \}$ is by a factor of $2 ^ { 4 / 5 } \left( 1 . 7 4 1 \right)$ for $\kappa = O ( \epsilon )$ and $2 ^ { 2 / 5 } \left( 1 . 3 1 9 \right)$ for $\kappa = O ( \epsilon ^ { 3 / 2 } )$ . The observed decrease in the neighborhood of $\epsilon / h \simeq 3$ is by a factor greater than four for the first case, by a factor of two for the second. Table V gives results that extend into the asymptotic regime. In order to hasten convergence κ was set to twice that of Table IV. $| \phi | _ { \operatorname* { m a x } } ( t )$ and the flow resistance are shown. $| \phi | _ { \operatorname* { m a x } } ( t )$ decreases as expected. For example, between $\epsilon / h$ of 7.516 and 8.634, $| \phi | _ { \operatorname* { m a x } } ( t )$ for $\kappa = O ( \epsilon )$ decreases by a factor of 1.74, for $\kappa = O ( \epsilon ^ { 3 / 2 } )$ it decreases by 1.32. The flow resistance decreases like $h ^ { 8 / 5 }$ for $\kappa = O ( \epsilon )$ , like $h ^ { 6 / 5 }$ for $\kappa = O ( \epsilon ^ { 3 / 2 } )$ . Note that the flow resistance is several orders of magnitude lower than for the second order method. 

The error in the potential due to convective truncation error does not have a directly deleterious effect on computations. This is because it is not as important as the error in the dissipation, the flow forcing, or in the surface energy, which are all much smaller. These latter errors are truer gauges of how the flow and interface are being distorted by the grid. Thus, the $O ( h ^ { 2 / 5 } )$ error in the potential for the fourth order method with $\kappa = O ( \epsilon ^ { 3 / 2 } )$ is acceptable, because the flow resistance error is only $O ( h ^ { 6 / 5 } )$ . 


TABLE V Instantaneous Maximum Chemical Potential and the Flow Resistance as a Function of Grid Resolution, for a Convected Interface Discretized Using Fourth Order Methods


<table><tr><td>Interface cellwidth <eq>\epsilon/h</eq></td><td><eq>|\phi|_{\text{max}}(t)</eq><eq>\kappa = O(\epsilon)</eq></td><td>Flow resistance<eq>\kappa = O(\epsilon)</eq></td><td><eq>|\phi|_{\text{max}}(t)</eq><eq>\kappa = O(\epsilon^{3/2})</eq></td><td>Flow resistance<eq>\kappa = O(\epsilon^{3/2})</eq></td></tr><tr><td>3.272</td><td>.380/.858</td><td>-.0184/.0213</td><td></td><td></td></tr><tr><td>3.758</td><td>.257/.473</td><td>-.00197/.00284</td><td></td><td></td></tr><tr><td>4.316</td><td>.156/.252</td><td><eq>-0.56/3.31 \times 10^{-4}</eq></td><td></td><td></td></tr><tr><td>4.958</td><td>.085/.134</td><td><eq>3.23/5.54 \times 10^{-5}</eq></td><td></td><td></td></tr><tr><td>5.696</td><td>.0488/.0719</td><td><eq>1.34/1.48 \times 10^{-5}</eq></td><td>.148/.217</td><td><eq>4.08/4.47 \times 10^{-5}</eq></td></tr><tr><td>6.543</td><td>.0293/.0393</td><td><eq>4.52/4.61 \times 10^{-6}</eq></td><td>.117/.157</td><td><eq>1.81/1.84 \times 10^{-5}</eq></td></tr><tr><td>7.516</td><td>.0172/.0218</td><td><eq>1.48 \times 10^{-6}</eq></td><td>.091/.115</td><td><eq>7.82 \times 10^{-6}</eq></td></tr><tr><td>8.634</td><td>.0102/.0122</td><td><eq>4.84 \times 10^{-7}</eq></td><td>.0718/.0854</td><td><eq>3.38 \times 10^{-6}</eq></td></tr></table>

# 10. CAPILLARY WAVES

Capillary wave computations provide a test of the numerics of the surface tension momentum forcing. This section considers small-amplitude capillary waves on a plane interface. Sections 8 and 9 have already shown that the fourth order method is far superior to the second order. Accordingly, only the fourth-order method will be discussed from now on. The problem is to calculate linear capillary wave frequencies as a function of wavelength, fluid viscosities, densities, and other fluid parameters. Numerical results can be compared to analytic results that are available for the original sharp-interface flow and for semi-analytic results for the linearized diffuse-interface model flow. An analytical expression is available for sharp-interface viscous capillary-wave frequencies and the frequency eigenvalue problem for capillary waves on plane diffuse interfaces is easily solvable numerically as a one-dimensional boundary value problem. 

The previous two sections have given an indication of how difficult it is even in one dimension to reach regimes in which numerical error is decaying in true asymptotic fashion. In two dimensions these regimes, as a practical matter, are completely unreachable. In this section this problem is partially got around by comparing numerical diffuseinterface frequencies to exact diffuse-interface frequencies and then comparing these exact diffuse-interface frequencies to those of the actual sharp interface. There are, in a way, two regimes of convergence, the “practical” and the asymptotic. The two have very different behaviors with the practical regime showing much more rapid convergence. It is easy to solve the frequency eigenvalue problem in the practical regime and sometimes difficult but possible to solve it in the asymptotic regime. In the latter, one finds the asymptotic rate of convergence of diffuse-interface to sharp-interface frequencies. It is also possible to show the “practical” convergence of two-dimensional numerical frequencies to exact diffuse interface frequencies in the “practical” regime. This practical convergence is much faster than the asymptotic numerical convergence hypothesized in Section 5. This hypothesis has been supported by the results reported in Sections 8 and 9 but these results have also shown that the coefficients for the asymptotic error, as indicated by Tables II and V, are extremely small. Essentially, as will be shown below, these asymptotic error terms are so small that they are invisible in capillary wave simulations except at impossibly fine grids. 

The exact physical system under consideration is capillary waves on an infinite plane interface. The interface separates two fluids with identical viscosities and densities and it runs along the x axis. To make numerical computations easier the system is bounded by nostress walls at $y = \pm y _ { s }$ . The boundary conditions at the walls are that the vertical velocity v is zero and the horizontal velocity u and C and $\phi$ are symmetric. The sharp-interface equation for the frequency ω is 

$$
\omega^ {2} = \frac {1}{2} \frac {\sigma k ^ {3}}{\rho} \left(\tanh k y _ {s} - \frac {k}{l} \tanh l y _ {s}\right), \quad l = \sqrt {k ^ {2} - i \omega / v}. \tag {10.1}
$$

k is the (real) wavenumber. ν is the kinematic viscosity. i is the square root of −1. ω is complex, with its real part being the frequency and its imaginary part the damping rate. This is an implicit equation as ω also appears in the viscous-related term l. 

The diffuse-interface eigenvalue equations are derived by linearizing (2.8)–(2.10) around a motionless plane interface and then assuming solutions of the form $\{ u , v , S , C , \phi \} $ $\exp ^ { - i \omega t }$ {sin kx u, cos kx v, cos kx S, cos kx C, cos kx φ}. The result can then be rearranged as a system of two fourth-order ordinary differential equations, one for v and one for φ: 

$$
\begin{array}{l} \mu \left(\frac {d ^ {2}}{d y ^ {2}} - k ^ {2}\right) ^ {2} v + i \omega \rho \left(\frac {d ^ {2}}{d y ^ {2}} - k ^ {2}\right) v = k ^ {2} \bar {C} _ {y} \phi (10.2a) \\ \alpha \kappa \left(\frac {d ^ {2}}{d y ^ {2}} - k ^ {2}\right) ^ {2} \phi - \beta \kappa \Psi^ {\prime \prime} (\bar {C}) \left(\frac {d ^ {2}}{d y ^ {2}} - k ^ {2}\right) \phi - i \omega \phi \\ = - \alpha k ^ {2} \bar {C} _ {y} v + 2 \alpha \bar {C} _ {y y} v _ {y} + \alpha \bar {C} _ {y} v _ {y y}. (10.2b) \\ \end{array}
$$

$\bar { C } ( y )$ is the $C ( y )$ of the unperturbed diffuse interface. Equations (10.2a) and (10.2b) were discretized using fourth-order differences. $\omega ( k , \sigma , \rho , \nu , \epsilon , \kappa )$ was found to eight-digit accuracy via shooting techniques coupled with Newton–Raphson iterations. 

Two-dimensional calculations of (2.8–2.10) were made for the particular case of $\sigma =$ 30 dynes/cm, ν = 1 cP, ρ = 1 gram/cm3, k = π, and $y _ { s } = 1 / 2 \ \mathrm { c m }$ . The calculations took advantage of horizontal symmetry and were of half a wavelength. The domain of the calculation was thus the square $0 \leq x \leq 1 , - 1 / 2 \leq y \leq + 1 / 2$ . For these parameters the exact sharp-interface frequency is $\omega = 2 0 . 1 0 3 1 3 - . 5 7 9 8 6 i$ . Calculations were time-accurate and were made on $1 6 \times 1 6$ up to $2 5 6 \times 2 5 6$ grids. The initial condition was a finite amplitude disturbance (velocity zero but the interface perturbed from planar). This sets off a capillary wave that gradually decays to being linear. Frequencies and decay rates were estimated by calculating and storing the kinetic energy at each time step and then computing times between kinetic energy peaks and decreases in amplitude from one peak to the next. Figure 1 shows the kinetic energy history of a typical $1 6 \times 1 6$ calculation. After some initial irregularity it settles into a weakly amplitude-dependent periodicity and decay rate. Each calculation was continued until the wave’s periodicity and decay rate became—to at least 4 figures—time independent. This typically took about 60 periods, during which the kinetic energy would decay by about 8 orders of magnitude. 

Both the eigenvalue and two-dimensional calculations were made using $\Psi ( C ) =$ $( C + 1 / 2 ) ^ { 2 } ( C - 1 / 2 ) ^ { 2 }$ . Tables VI and VII show exact diffuse-interface eigenfrequencies as a function of ² as $\epsilon  0$ for, respectively, $\kappa \propto \epsilon$ and $\kappa \propto \epsilon ^ { 2 }$ . The eigenfrequencies are shown in terms of their two components, the real frequency in radians/s, and the damping rate. The two sequences converge about equally fast down to $\epsilon \simeq . 0 1$ , at which point the relative error is less than 0.1%. In the asymptotic regime convergence is $O ( \epsilon )$ for $\kappa \propto \epsilon$ and $O ( \epsilon ^ { 2 } )$ for $\kappa \propto \epsilon ^ { 2 }$ . For $\kappa \propto \epsilon$ the asymptotic regime begins at about $\epsilon \simeq . 0 0 1$ . (The grid required to resolve this regime, if uniform, would have something like 5000 × 5000 points.) The asymptotic regime begins sooner for $\kappa \propto \epsilon ^ { 2 }$ because asymptotic surface tension errors, convection errors, and diffusive errors are then all the same order. In general, for $\kappa = O ( \epsilon ^ { \delta } ) , 1 \leq \delta < 2$ , the asymptotic rate of convergence is controlled by the rate of diffusion and is $O ( \kappa )$ . The comparatively high accuracy of these results is most likely due to the absence of solubility effects, because of the plane interface. It may also be due, because of the symmetry of v at the interface, to the absence of significant interface convective straining. This lack of straining makes it possible to obtain convergence for 

![](images/14bbcee32624a7ff04764d35c9bcf92cc7987ea1ccd22d52cea71dcedb736576.jpg)



FIG. 1. Time-dependent kinetic energy of a decaying finite-amplitude capillary wave. The wave was calculated on a $1 6 \times 1 6$ grid.


$2 \leq \delta \leq 3$ . This regime is unuseable for simulations because of effects from numerical convection. 

Tables VIII through XII show numerical results using the 4th order mehrstellungen method for different $\epsilon ( h )$ and $\kappa ( h )$ . The tables show the numerical complex frequencies compared to exact diffuse-interface frequencies. ² is proportional to $h ^ { 2 / 3 }$ in Table VII, to $h ^ { 4 / 5 }$ in Table VIII then, in order, to $h ^ { 6 / 7 } , h ^ { 8 / 9 }$ , and h. κ is proportional to $\epsilon ^ { 2 } \propto h ^ { 4 / 3 }$ in 


TABLE VI Diffuse-Interface Frequency and Damping Rate as a Function of ² for $\kappa = O ( \epsilon )$


<table><tr><td><eq>\epsilon</eq></td><td><eq>\kappa \times 10^{7}</eq></td><td>Freq</td><td>Freq error</td><td>Damp</td><td>Damp error</td></tr><tr><td>3/16</td><td>405.51</td><td>18.739</td><td></td><td>-.7666</td><td></td></tr><tr><td>3/32</td><td>202.75</td><td>19.498</td><td></td><td>-.5227</td><td></td></tr><tr><td>3/64</td><td>101.38</td><td>19.891</td><td></td><td>-.5154</td><td></td></tr><tr><td>3/128</td><td>50.688</td><td>20.03958</td><td></td><td>-.553939</td><td></td></tr><tr><td>3/256</td><td>25.344</td><td>20.08441</td><td><eq>-1.872 \times 10^{-2}</eq></td><td>-.573489</td><td><eq>6.370 \times 10^{-3}</eq></td></tr><tr><td>3/512</td><td>12.672</td><td>20.09723</td><td><eq>-5.897 \times 10^{-3}</eq></td><td>-.579231</td><td><eq>6.288 \times 10^{-4}</eq></td></tr><tr><td>3/1024</td><td>6.3360</td><td>20.10106</td><td><eq>-2.069 \times 10^{-3}</eq></td><td>-.580318</td><td><eq>-4.580 \times 10^{-4}</eq></td></tr><tr><td>3/2048</td><td>3.1680</td><td>20.10231</td><td><eq>-8.137 \times 10^{-4}</eq></td><td>-.580300</td><td><eq>-4.402 \times 10^{-4}</eq></td></tr><tr><td>3/4096</td><td>1.5840</td><td>20.10278</td><td><eq>-3.516 \times 10^{-4}</eq></td><td>-.580135</td><td><eq>-2.753 \times 10^{-4}</eq></td></tr><tr><td>3/8192</td><td>0.7920</td><td>20.10297</td><td><eq>-1.620 \times 10^{-4}</eq></td><td>-.580012</td><td><eq>-1.517 \times 10^{-4}</eq></td></tr><tr><td>0</td><td>0</td><td>20.10313</td><td>0</td><td>-.579860</td><td>0</td></tr></table>


TABLE VII Diffuse-Interface Frequency and Damping Rate as a Function of ² for $\kappa = O ( \epsilon ^ { 2 } )$


<table><tr><td><eq>\epsilon</eq></td><td><eq>\kappa \times 10^{7}</eq></td><td>Freq</td><td>Freq error</td><td>Damp</td><td>Damp error</td></tr><tr><td>3/16</td><td>405.51</td><td>18.739</td><td></td><td>-.7666</td><td></td></tr><tr><td>3/32</td><td>101.38</td><td>19.506</td><td></td><td>-.5870</td><td></td></tr><tr><td>3/64</td><td>25.344</td><td>19.887</td><td></td><td>-.5434</td><td></td></tr><tr><td>3/128</td><td>6.3360</td><td>20.03927</td><td></td><td>-.556827</td><td></td></tr><tr><td>3/256</td><td>1.5840</td><td>20.08641</td><td><eq>-1.672 \times 10^{-2}</eq></td><td>-.571300</td><td><eq>8.560 \times 10^{-3}</eq></td></tr><tr><td>3/512</td><td>0.3960</td><td>20.09890</td><td><eq>-4.227 \times 10^{-3}</eq></td><td>-.577348</td><td><eq>2.512 \times 10^{-3}</eq></td></tr><tr><td>3/1024</td><td>0.0990</td><td>20.10205</td><td><eq>-1.080 \times 10^{-3}</eq></td><td>-.579206</td><td><eq>6.540 \times 10^{-4}</eq></td></tr><tr><td>3/2048</td><td>0.0248</td><td>20.10285</td><td><eq>-2.822 \times 10^{-4}</eq></td><td>-.579702</td><td><eq>1.576 \times 10^{-4}</eq></td></tr><tr><td>3/4096</td><td>0.0062</td><td>20.10305</td><td><eq>-7.575 \times 10^{-5}</eq></td><td>-.579824</td><td><eq>3.517 \times 10^{-5}</eq></td></tr><tr><td>3/8192</td><td>0.0015</td><td>20.10311</td><td><eq>-2.091 \times 10^{-5}</eq></td><td>-.579852</td><td><eq>6.990 \times 10^{-6}</eq></td></tr><tr><td>0</td><td>0</td><td>20.10313</td><td>0</td><td>-.579860</td><td>0</td></tr></table>


TABLE VIII Plane Interface Capillary Wave Oscillation Frequency and Damping



Rate: ² h2/3, κ h2 h4/3


<table><tr><td>No. pts</td><td><eq>\epsilon</eq></td><td>Cell width</td><td><eq>\kappa \times 10^{7}</eq></td><td>Freq exact</td><td>Freq numerical</td><td>Damp exact</td><td>Damp numerical</td></tr><tr><td><eq>16 \times 16</eq></td><td>.1875</td><td>3.000</td><td>405.5</td><td>18.739</td><td>18.607</td><td>-.7666</td><td>-.6926</td></tr><tr><td><eq>32 \times 32</eq></td><td>.1181</td><td>3.780</td><td>160.9</td><td>19.302</td><td>19.286</td><td>-.6294</td><td>-.6178</td></tr><tr><td><eq>64 \times 64</eq></td><td>.07441</td><td>4.762</td><td>63.86</td><td>19.668</td><td>19.663</td><td>-.5607</td><td>-.5585</td></tr><tr><td><eq>128 \times 128</eq></td><td>.04687</td><td>6.000</td><td>25.34</td><td>19.887</td><td>19.886</td><td>-.5434</td><td>-.5429</td></tr><tr><td><eq>256 \times 256</eq></td><td>.02953</td><td>7.560</td><td>10.06</td><td>20.006</td><td></td><td>-.5507</td><td></td></tr><tr><td><eq>512 \times 512</eq></td><td>.01860</td><td>9.254</td><td>3.991</td><td>20.062</td><td></td><td>-.5626</td><td></td></tr><tr><td><eq>\infty \times \infty</eq></td><td>0</td><td><eq>\infty</eq></td><td>0</td><td>20.103</td><td></td><td>-.5799</td><td></td></tr></table>


TABLE IX Plane Interface Capillary Wave Oscillation Frequency and Damping



Rate: ² ∝ h4/5, κ ∝ ²3/2 ∝ h6/5


<table><tr><td>No. pts</td><td><eq>\epsilon</eq></td><td>Cell width</td><td><eq>\kappa \times 10^{7}</eq></td><td>Freq exact</td><td>Freq numerical</td><td>Damp exact</td><td>Damp numerical</td></tr><tr><td><eq>16 \times 16</eq></td><td>.1875</td><td>3.000</td><td>405.5</td><td>18.739</td><td>18.607</td><td>-.7666</td><td>-.6926</td></tr><tr><td><eq>32 \times 32</eq></td><td>.1077</td><td>3.446</td><td>176.5</td><td>19.384</td><td>19.342</td><td>-.5743</td><td>-.5615</td></tr><tr><td><eq>64 \times 64</eq></td><td>.06185</td><td>3.959</td><td>76.83</td><td>19.773</td><td>19.765</td><td>-.5209</td><td>-.5186</td></tr><tr><td><eq>128 \times 128</eq></td><td>.03552</td><td>4.547</td><td>33.44</td><td>19.972</td><td>19.972</td><td>-.5346</td><td>-.5342</td></tr><tr><td><eq>256 \times 256</eq></td><td>.02040</td><td>5.223</td><td>14.56</td><td>20.055</td><td></td><td>-.5577</td><td></td></tr><tr><td><eq>512 \times 512</eq></td><td>.01172</td><td>6.000</td><td>6.336</td><td>20.086</td><td></td><td>-.5713</td><td></td></tr><tr><td><eq>\infty \times \infty</eq></td><td>0</td><td><eq>\infty</eq></td><td>0</td><td>20.103</td><td></td><td>-.5799</td><td></td></tr></table>


TABLE X Plane Interface Capillary Wave Oscillation Frequency and Damping Rate: $\epsilon \propto h ^ { 6 / 7 } , \kappa \propto \epsilon ^ { 4 / 3 } \propto h ^ { 8 / 7 }$


<table><tr><td>No. pts</td><td><eq>\epsilon</eq></td><td>Cell width</td><td><eq>\kappa \times 10^{7}</eq></td><td>Freq exact</td><td>Freq numerical</td><td>Damp exact</td><td>Damp numerical</td></tr><tr><td><eq>16 \times 16</eq></td><td>.1875</td><td>3.000</td><td>405.5</td><td>18.739</td><td>18.607</td><td>-.7666</td><td>-.6926</td></tr><tr><td><eq>32 \times 32</eq></td><td>.1035</td><td>3.312</td><td>183.6</td><td>19.417</td><td>19.342</td><td>-.5561</td><td>-.5420</td></tr><tr><td><eq>64 \times 64</eq></td><td>.05714</td><td>3.657</td><td>83.16</td><td>19.812</td><td>19.777</td><td>-.5146</td><td>-.5116</td></tr><tr><td><eq>128 \times 128</eq></td><td>.03154</td><td>4.038</td><td>37.66</td><td>19.997</td><td>19.986</td><td>-.5388</td><td>-.5382</td></tr><tr><td><eq>256 \times 256</eq></td><td>.01741</td><td>4.458</td><td>17.06</td><td>20.067</td><td></td><td>-.5630</td><td></td></tr><tr><td><eq>512 \times 512</eq></td><td>.00961</td><td>4.922</td><td>7.724</td><td>20.091</td><td></td><td>-.5744</td><td></td></tr><tr><td><eq>\infty \times \infty</eq></td><td>0</td><td><eq>\infty</eq></td><td>0</td><td>20.103</td><td></td><td>-.5799</td><td></td></tr></table>

Table VII then to $\epsilon ^ { 3 / 2 } \propto h ^ { 6 / 5 } , \epsilon ^ { 4 / 3 } \propto h ^ { 8 / 7 } , \epsilon ^ { 5 / 4 } \propto h ^ { 9 / 8 }$ , and $\epsilon \propto h$ . The theoretical asymptotic errors relative to diffuse-interface results are, respectively $O ( h ^ { 4 / 3 } )$ (because of the absence of interace straining), $O ( h ^ { 4 / 5 } ) , O ( h ^ { 4 / 7 } ) , O ( h ^ { 4 / 9 } )$ , and $O ( 1 )$ . Except for $\epsilon \propto h$ , the results show that the practical convergence of numerical frequencies to diffuse-interface frequencies is much faster than asymptotic. Convergence of numerical frequencies to diffuse-interface frequencies becomes slower as $\epsilon$ becomes more nearly proportional to h but it remains very rapid, much faster than linear, even for $\epsilon \propto h ^ { 8 / 9 }$ . The overall convergence to the exact sharp-interface frequency is also much faster than linear. The case $\epsilon \propto h ^ { 6 / 7 }$ shows the fastest frequency convergence; $\epsilon \propto h ^ { 8 / 9 }$ shows the fastest damping convergence. This would no longer hold true it the computations were made on finer and finer grids and the true asymptotic regime were reached. The calculations $\epsilon \propto h$ may be in the asymptotic regime. The results for this show an $O ( 1 )$ error and indicate that there may be divergence. This error is dominated by numerical grid and convection effects. 

These same calculations have been made with a number of other second- and fourthorder-accurate discretizations of the surface tension forcing and of the $\nabla ^ { 2 } \phi$ term. So far, the discretizations given in Section 7 have been found to be both the simplest and the best. The computed damping rate has been found to be much more sensitive than the real frequency to the discretization and also to such factors as number of pressure iterations or the number of iterations used to solve the discrete chemical potential equation (Eq. (7.10)). The computed damping rate seems to be most dependent on the discretization used for the surface tension forcing. It is significantly but secondarily affected by the discretization of $\nabla ^ { 2 } \phi$ . 


TABLE XI Plane Interface Capillary Wave Oscillation Frequency and Damping $\mathbf { R a t e } \colon \epsilon \propto h ^ { 8 / 9 } , \kappa \propto \epsilon ^ { 5 / 4 } \propto h ^ { 1 0 / 9 }$


<table><tr><td>No. pts</td><td><eq>\epsilon</eq></td><td>Cell width</td><td><eq>\kappa \times 10^{7}</eq></td><td>Freq exact</td><td>Freq numerical</td><td>Damp exact</td><td>Damp numerical</td></tr><tr><td><eq>16 \times 16</eq></td><td>.1875</td><td>3.000</td><td>405.5</td><td>18.739</td><td>18.607</td><td>-.7666</td><td>-.6926</td></tr><tr><td><eq>32 \times 32</eq></td><td>.1013</td><td>3.240</td><td>187.7</td><td>19.436</td><td>19.327</td><td>-.5473</td><td>-.5322</td></tr><tr><td><eq>64 \times 64</eq></td><td>.05468</td><td>3.500</td><td>86.90</td><td>19.831</td><td>19.745</td><td>-.5131</td><td>-.5091</td></tr><tr><td><eq>128 \times 128</eq></td><td>.02953</td><td>3.780</td><td>40.23</td><td>20.009</td><td>19.946</td><td>-.5419</td><td>-.5406</td></tr><tr><td><eq>256 \times 256</eq></td><td>.01595</td><td>4.082</td><td>18.62</td><td>20.072</td><td></td><td>-.5657</td><td></td></tr><tr><td><eq>512 \times 512</eq></td><td>.00861</td><td>4.409</td><td>8.622</td><td>20.093</td><td></td><td>-.5758</td><td></td></tr><tr><td><eq>\infty \times \infty</eq></td><td>0</td><td><eq>\infty</eq></td><td>0</td><td>20.103</td><td></td><td>-.5799</td><td></td></tr></table>


TABLE XII Plane Interface Capillary Wave Oscillation Frequency and Damping Rate: ² ∝ h, κ ∝ h


<table><tr><td>No. pts</td><td><eq>\epsilon</eq></td><td>Cell width</td><td><eq>\kappa \times 10^{7}</eq></td><td>Freq exact</td><td>Freq numerical</td><td>Damp exact</td><td>Damp numerical</td></tr><tr><td><eq>16 \times 16</eq></td><td>.1875</td><td>3.000</td><td>405.5</td><td>18.739</td><td>18.607</td><td>-.7666</td><td>-.6926</td></tr><tr><td><eq>32 \times 32</eq></td><td>.09375</td><td>3.000</td><td>202.8</td><td>19.498</td><td>19.101</td><td>-.5227</td><td>-.5009</td></tr><tr><td><eq>64 \times 64</eq></td><td>.04688</td><td>3.000</td><td>101.4</td><td>19.891</td><td>18.256</td><td>-.5154</td><td>-.4913</td></tr><tr><td><eq>128 \times 128</eq></td><td>.02344</td><td>3.000</td><td>50.69</td><td>20.040</td><td>12.121</td><td>-.5539</td><td>-.4424</td></tr><tr><td><eq>\infty \times \infty</eq></td><td>0</td><td>3.000</td><td>0</td><td>20.103</td><td></td><td>-.5799</td><td></td></tr></table>

# 11. RAYLEIGH–TAYLOR INSTABILITIES

The Rayleigh–Taylor instability can occur when a denser fluid lies over a lighter. Waves form on the interface, increase in amplitude, and transform into plumes. Computations in this section are of high-capillary-number, large-deformation flows. Analytic solutions are not available; convergence is examined through visual inspection of grid-refined results. The computed flows include near-singularities such as plume breaking, droplet formation, droplet coalescence, contact line flow, and the formation of wall films. The calculations in this section are fourth-order and use $\epsilon \propto h ^ { 4 / 5 }$ and $\kappa \propto \epsilon ^ { 3 / 2 } \propto h ^ { 6 / 5 }$ . This choice gives optimal asymptotic convergence $( O ( h ^ { 4 / 5 } )$ , see Section 5) while minimizing κ. As discussed in Section 7, the asymptotic error of the fourth-order compact method is the maximum of $\epsilon , h ^ { 4 } / \epsilon ^ { 4 }$ , and $h ^ { 2 } / \kappa$ . 

Figures 2a and 2b show the evolution of a Rayleigh–Taylor instability contained in a square box with no-slip walls. The box is ten by ten centimeters. The dense fluid has a density of 1.0, the light 0.9. Ths viscosity of both fluids is 1 poise. The less dense fluid occupies one-eighth of the box. The initial interface is flat except for a small perturbation at the box’s left wall. The figure shows results computed using four grids, $4 8 \times 4 8$ (top row), $9 6 \times 9 6$ (top middle), $1 9 2 \times 1 9 2$ (bottom middle), and 384 × 384 (bottom). Their 90% interface thicknesses are, respectively, 0.328 cm = 1.57 cellwidths, 0.188 cm = 1.81 cellwidths, $0 . 1 0 8 \mathrm { c m } = 2 . 0 8$ cellwidths, and 0.0621 cm = 2.38 cellwidths. The mobilities are $1 . 6 3 \times 1 0 ^ { - 4 } , 7 . 1 0 \times 1 0 ^ { - 5 } , 3 . 0 9 \times 1 0 ^ { - 5 }$ , and $1 . 3 4 \times 1 0 ^ { - 5 }$ . Contours are at intervals of 0.01, at $C = \pm 0 . 0 0 5 , C = \pm 0 . 0 1 5 , \mathrm { e t c }$ . This is done to make any fluid intermixing or interface profile deformation clearly visible. Except during interface breakups and coalescences interface profiles are deformed very little. Because of the high contour density interface regions are marked by solid back. Usually the black extends from C = −0.495 to $C = + 0 . 4 9 5$ , so 99% thickness is shown. 

Results are shown for 6 different times, at 1.1, 1.65, 2.2, 2.75, 3.3, and 3.85 seconds. The first shows the instability as it initially amplifies and propagates toward the right wall. The instability grows at the left wall (second column) and, because of no-slip effects there, begins to plume away from it. A plume begins to form midway between the two sidewalls. Between the two plumes the lower fluid has almost completely drained away from the lower wall, leaving a viscous film. The third column shows the right plume fully formed. In the fourth column this plume has transformed into a rising, highly asymmetric drop. The left-wall plume has also separated and become the roughly circular droplet on the left. At the fifth time the two droplets have risen to the top. The more resolved calculations show them nestled against the top wall (not attached). Fluid is rising slowly along the two side walls. The right plume filament is snapping back to the lower wall. In the two coarser calculations the lower-wall film, due to disjoining-pressure instabilities, has split into droplets. 

![](images/586922f901a384c9917aba1329d37119485e9134f41990978b69b1b97e60d305.jpg)



FIG. 2. (a) Rayleigh–Taylor instability of a two-phase fluid. The computation was made on four grids, $4 8 \times 4 8 ( \mathrm { t o p } ) , 9 6 \times 9 6$ (top middle), $1 9 2 \times 1 9 2$ (bottom middle), and $3 8 4 \times 3 8 4$ (bottom). Plots are in chronological order reading from left to right. Contours are at 0.01 intervals. Interface smearing and distortion errors are thus shown in detail. The lighter fluid is initially confined to the bottom. The initial interface is flat except for a small perturbation at the left wall. The interface moves upward at that location and a wave propagates out toward the right wall (first column). The disturbance amplifies (second column), plumes (third column), and then, in Fig. 2b, separates into droplets (first column, Fig. 2b). (b) Continuation of Fig. 2a. The plumes separate into drops, which rise to the top wall. A thin film of light fluid is left on the bottom. At the coarser resolutions this breaks into droplets.


At the first three times the convergence of the calculation as grid size is increased is very clear. The $9 6 \times 9 6 , 1 9 2 \times 1 9 2$ , and $3 8 4 \times 3 8 4$ results are all essentially the same. The $4 8 \times 4 8$ calculation takes longer to form the right plume and its lower-wall film breaks early. 

![](images/48b87b888fcb1c775da1091a5c7f894bd5a46d3492e51b89a0053a80d28b4256.jpg)



FIG. 2—Continued


The main differences in the finer calculations are that intermixing and interface deformation errors become much smaller. 

The later times are affected by the various flow near-singularities that occur. The more resolved the calculation the longer, in general, that it takes for wall film breakdown, for plume separation, and for droplet breakups and coalescences. Whether this matters or not depends on what is wanted from the calculation. The three finer calculations are in very good agreement on the speed of the instability and its rate of vertical mass and energy transfer. At the last time shown the configuration of the two upper drops is close to the same for all three. There are differences; the left droplet in the $9 6 \times 9 6$ calculation has begun to attach to the wall, it remains separate on the two finer grids; the 384 × 384 calculation shows a small third droplet in contact and about to coalesce with the right drop, this has already occurred on the $9 6 \times 9 6$ and $1 9 2 \times 1 9 2$ grids. These differences most likely have only local, in both space and time, effects. 

![](images/f48aeb03d34fdad8a98b24e1b8871a4c5a97e0a93f24a81f24f0bdc1573dc98c.jpg)



FIG. 3. The propagation of a Rayleigh–Taylor instability wave. Contours are at 0.1 intervals, from −0.55 to +0.55. There are 10 contours between the two $C _ { \mathrm { b u l k p h a s e } }$ values. The wave is initiated by a small disturbance in interface position at the box’s middle and then propagates symmetrically towards the two sidewalls. Times shown are 1.125, 1.6875, 3.375, 3.9375, 4.5, and 5.625 s.


The rightward propagating instability wave seen at times 1 and 2 is strongly affected by the presence of sidewalls. To see that wave and its manner of propagation more clearly the instability has been calculated in a much longer, 60 by 10 centimeters, box. Results are shown in Figs. 3–5 for grids $3 8 4 \times 6 4 , 7 6 8 \times 1 2 8$ , and $1 5 7 6 \times 2 5 6$ . Their 90% interface widths are $0 . 2 6 0 \mathrm { c m } = 1 . 6 7$ cellwidths, 0.150 cm = 1.91 cellwidths, and $0 . 0 8 5 9 \mathrm { c m } = 2 . 2 0$ cellwidths. The mobilities are $1 . 1 5 \times 1 0 ^ { - 4 } , 5 . 0 3 \times 1 0 ^ { - 5 }$ , and $2 . 1 9 \times 1 0 ^ { - 5 }$ . Results are shown at times 1.125, 1.6875, 3.375, 3.9375, 4.5, and 5.625 seconds. Contours in the figures are at 0.1 intervals, from $C = - 0 . 4 5 \mathrm { t o } + 0 . 4 5$ . The disturbance is initiated at the box’s midpoint and propagates toward each side wall. Half the container is calculated; the results for the other half found by symmetry. Material properties are the same as for the previous calculation and the lighter fluid, as before, fills one-eighth of the container. 

![](images/5a5ac204017cd0d8fe7c1da4946d2aede89ea7524a9ba05f095f9fbc31679767.jpg)



FIG. 4. The same as Fig. 3 but at a resolution of 768 × 128 and with smaller ² and κ.


The three calculations are in very good agreement on important macroscopic quantities such as propagation speed, wavelength, rates of mass transfer, and the form of the propagation. The propagation speed of the instability is about 5% faster for the finest resolution than for the coarsest. Frame 3 shows the nature of the propagation fairly clearly. Plumes 0 and ±1 have already detached, plumes ±2 are fully formed, and plumes ±3 are beginning to form. In the 1576 × 256 calculation the filaments of plumes ±1 are snapping back to the lower wall, while the snapped back filament of plume 0 has become a droplet ready to detach. Between plumes ±2 and ±3 the lighter fluid has been nearly drained, leaving, as with Fig. 2, a thin film. Ahead of plume ±3 the interface looks almost undisturbed. A precursor, very low amplitude, capillary wave has already propagated to the wall but it is invisible at the scale of the figure. 

![](images/77dda9d47b9b3e81b2ffc24b6c7c4b743ece04cfc495222a898d3c6e0a0a39ce.jpg)



FIG. 5. The same as Fig. 4 but at a resolution of 1576 × 256 and with smaller ² and κ.


Areas of disagreement in the three calculations include small scale quantities such as plume-filament and film breakup times, which clearly have almost no impact on the instability, and coalescence times of drops at the container top. The differences in coalescence times have an effect on the appearance and pattern of drops against the top wall but very little on their general distribution—in all three cases the flow tends to cluster the upper wall light fluid slightly toward the center of the box. The finer two calculations are in very good agreement in the upper part of the third frame. Between the third and fourth frames the upper center droplet, which is already very elongated at time 3, splits into 3 droplets in the $7 6 8 \times 1 2 8$ calculation while remaining whole in the $1 5 7 6 \times 2 5 6$ calculation. This then leads to further differences in coalescences and splittings in the next two frames. 

# 12. CONCLUSION

Fluid-dynamical phase-field modeling is a new numerical/modeling approach to the computation of two-phase flows and one with great promise. It allows the use of common, easily analyzable and easily useable centered finite-volume, finite-difference, or finiteelement convection schemes. One of the major disadvantages of phase-field models has been the relatively large width of their interfaces. This paper has introduced a new compact method that allows accurate computations, as has been shown in the section on high-capillary-number Rayleigh–Taylor instabilities, with interfacial thicknesses that can be less than two cellwidths. This is less than the mollified interface width used with the CSF model and also less than the width of interfacial force distribution used in tracking/distributed force methods. 

The overall accuracy of the method is a very complicated issue. It is a function of three parameters, the interface width ², the mobility κ and the mesh spacing h. It is also a function of the rate of convergence of the phase-field model to sharp-interface results. The analysis has for the most part assumed that this convergence is $O ( \epsilon )$ . If this is true then the compact method discussed herein can be optimized to an asymptotic overall accuracy of only $h ^ { 4 / 5 }$ Fortunately, it appears that practical convergence is generally much faster than asymptotic convergence. Also, there is reason to believe that $O ( \epsilon ^ { 2 } )$ models can be developed. If so, then better than $O ( h )$ asymptotic accuracy can be obtained. 

It might seem that the various difficult issues raised in the convergence analysis can be avoided with the CSF and distributed force models. This is probably not so. Both models have an implied analytic model of continuum forcing. The accuracy of these models is dependent on the model interface thickness ². The rate of convergence for these analytic models as $\epsilon  0$ has not been established for either case. The rate of convergence of the VOF-CSF method, particularly its curvature calculations, may be a function of $h / \epsilon$ rather than just h. Finally, there are aspects of the CSF method that are analogous to the phase-field method’s use of diffusivity. The effects of VOF and level-set convection on surface energy, how they control or don’t control it, have not yet been analyzed. 

The capillary wave test of Section 10 is an excellent and very difficult test for discrete diffuse interface surface tension models. The propagation rate and frequency of capillary waves are determined by the physics of energy transformations from kinetic to surface energy and back again. In real two-phase flows the creation of surface energy by convection is always equal and opposite to the creation of kinetic energy by surface tension. In diffuseinterface models interfacial energy creation occurs through the convection equation for the color function. Unfortunately, it is easy to formulate discrete systems in which the discrete convection of C and the discrete surface tension forcing by C have incompatibilities. The author’s experience is that the linear zero-amplitude limit brings these incompatibilities out. As discussed in Section 8, rapid convection is less susceptible to grid effects than is slow. In the same way, large-amplitude capillary waves are less susceptible to these effects than are small-amplitude waves. Also, any tendency to form parasitic flows becomes more noticeable as wave amplitude is reduced. 

Section 11 touched on the difficulties involved in calculating flow near-singularities such as interface breakups and coalescences. Times at which these events occur can be very sensitive to ². The actual events happen quickly but the time leading up them, for example, viscous drainage leading up to coalescence, can scale like $1 / \epsilon$ . This is not just a numerical or modeling problem but is a difficulty that is observed in real flows. The time scales for coalescences and interface break-ups are in reality highly variable from fluid system to fluid system. For example, air bubbles often coalesce fairly quickly, lava lamp drops never do. Coalescence times depend not just on viscous drainage times but on both longrange (micro-scale) and short range (nano-scale) electrostatic and molecular interactions. An advantage of the phase-field approach is that, if need be, material-dependent models of these small scale electrostatic energies and potentials can be included in its energy formulations. 

The ability to do such micro-scale modeling is one of the great strengths of the phasefield approach. For example, phase-field models are applicable to the simulation of complex fluids such as micro-emulsions (Lamura et al. [15]) and they may be useful for studying two-phase micro-fluidic flows in which electrophoretic or other effects play a part. Also, the phase-field method’s relatively analytical grounding makes it useful for the study of twophase flow singularities. Initial work in this area has been done by Lowengrub and Trusinovsky [16] for interface coalescences and break-ups and by the author [10] for moving contact lines. Finally, the method can be of use for “DNS” studies of 10-100 nanometer flows (Jacqmin [10]). At this scale, actual interface thicknesses (about one to two nanometers) can be included. 

# REFERENCES



1. D. M. Anderson and G. B. McFadden, A diffuse-interface description of internal waves in a near-critical fluid, Phys. Fluids 9, 1870 (1997). 





2. L. K. Antanovskii, A phase field model of capillarity, Phys. Fluids 7, 747 (1995). 





3. P. W. Bates and P. C. Fife, The dynamics of nucleation for the Cahn–Hilliard equation, SIAM J. Appl. Math. 53, 990 (1993). 





4. J. Brackbill, D. B. Kothe, and C. Zemach, A continuum method for modeling surface tension, J. Comput. Phys. 100, 335 (1992). 





5. J. W. Cahn and J. E. Hilliard, Free energy of a nonuniform system. III. Nucleation in a two-component incompressible fluid, J. Chem. Phys. 31, 688 (1959). 





6. R. Chella and J. Vi˜nals, Mixing of a two-phase fluid by cavity flow, Phys. Rev. E 53, 3832 (1996). 





7. D. Jacqmin, Three-Dimensional Computations of Droplet Collisions, Coalescence, and Droplet/Wall Interactions Using a Continuum Surface Tension Method, AIAA 95-0883, presented at the 33rd Aerospace Sciences Meeting, Reno, NV, 1995. 





8. D. Jacqmin, An Energy Approach to the Continuum Surface Tension Method: Application to Droplet Coalescences and Droplet/Wall interactions, presented at the ASME IMECE, San Francisco, CA, 1995. 





9. D. Jacqmin, An Energy Approach to the Continuum Surface Method, AIAA 96-0858, presented at the 34th Aerospace Sciences Meeting, Reno, NV, 1996. 





10. D. Jacqmin, Contact line dynamics of a diffuse fluid interface, J. Fluid Mech., in press. 





11. D. Jasnow and J. Vi˜nals, Coarse-grained description of thermo-capillary flow, Phys. Fluids 8, 660 (1996). 





12. D. D. Joseph, Fluid dynamics of two miscible liquids with diffusion and gradient stresses, Eur. J. Mech. B/Fluids 9, 565 (1990). 





13. D. B. Kothe, W. J. Rider, S. J. Mosso, and J. S. Brock, Volume Tracking of Interfaces Having Surface Tension in Two and Three Dimensions, AIAA 96-0859, presented at the 34th Aerospace Sciences Meeting, Reno, NV, 1996. 





14. B. Lafaurie, C. Nardone, R. Scardoveli, S. Zaleski, and G. Zanetti, Modeling merging and fragmentation in multiphase flows with SURFER, J. Comput. Phys. 113, 134 (1994). 





15. A. Lamura, G. Gonnela, and J. Yeomans, Modeling the dynamics of amphiphilic fluids, Int. J. Mod. Phys. C 9, 1469 (1998). 





16. J. Lowengrub and I. Truskinovsky, Cahn–Hilliard fluids and topological transitions, Proc. R. Soc. London A 454, 2617 (1998). 





17. B. T. Nadiga and S. Zaleski, Investigations of a two-phase fluid model, Eur. J. Mech. B/Fluids 15, 885 (1996). 





18. R. H. Nochetto, M. Paolini, and C. Verdi, A dynamic mesh algorithm for curvature dependent evolving interfaces, J. Comput. Phys. 123, 296 (1996). 





19. Y. Oono and S. Puri, Study of phase-separation dynamics by use of cell dynamical systems. I. Modeling, Phys. Rev. A 38, 434 (1988). 





20. W. J. Rider, D. B. Kothe, S. J. Moddo, and J. H. Cerutti, Accurate Solution Algorithms for Incompressible Multiphase Flows, AIAA 95-0699, presented at the 33rd Aerospace Sciences Meeting, Reno, NV, 1995. 





21. M. Sussman, P. Smereka, and S. J. Osher, A level set approach for computing solutions to incompressible two-phase flow, J. Comput. Phys. 114, 146 (1994). 





22. M. Sussman, E. Fatemi, P. Smereka, and S. J. Osher, An improved level set method for incompressible two-phase flows, J. Comput. Fluids 27, 663 (1998). 





23. S. O. Unverdi and G. Tryggvason, A front-tracking method for viscous, incompressible, multi-fluid flows, J. Comput. Phys. 100, 25 (1992). 





24. J. D. van der Waals, The Thermodynamic Theory of Capillarity Flow under the Hypothesis of a Continuous Variation of Density (Verhandel/Konink. Akad. Weten., 1893), Vol. 1; English translation, J. Statist. Phys. 20, 197. 

