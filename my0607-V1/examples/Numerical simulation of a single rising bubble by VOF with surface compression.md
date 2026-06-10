# Numerical simulation of a single rising bubble by VOF with surface compression

J. Klostermann*,†, K. Schaake and R. Schwarze 

Institute of Mechanics and Fluid Dynamics, TU Bergakademie Freiberg, Freiberg, Germany 

## SUMMARY

The capability of the direct volume of fluid method for describing the surface dynamics of a free twodimensional rising bubble is evaluated using quantities of a recently published benchmark. The model equations are implemented in the open source computational fluid dynamics library OpenFOAM®. Here, a main ingredient of the numerical method is the so-called surface compression that corrects the fluxes near the interface between two phases. The application of this method with respect to two test cases of a benchmark is considered in the main part. The test cases differ in physical properties, thus in different surface tension effects. The quantities centre of mass position, circularity and rise velocity are tracked over time and compared with the ones given in the benchmark. For test case one, where surface tension effects are more pre-eminent, deviations from the benchmark results become more obvious. However, the flow features are still within reasonable range. Nevertheless, for test case two, which has higher density and viscosity ratios and above all a lower influence of the surface tension force, good agreement compared with the benchmark reference results is achieved. This paper demonstrates the good capabilities of the direct volume of fluid method with surface compression with regard to the preservation of sharp interfaces, boundedness, mass conservation and low computational time. Some limitation regarding the occurrence of parasitic currents, bad pressure jump prediction and bad grid convergence have been observed. With these restrictions in mind, the method is suitable for the simulation of similar two-phase flow configurations. Copyright © 2012 John Wiley & Sons, Ltd. 

Received 1 September 2010; Revised 16 April 2012; Accepted 30 April 2012 

KEY WORDS: VOF; surface tension force; parasitic/spurious currents; direct method; free surface; rising bubble benchmark 

## 1. INTRODUCTION

Interfacial flows of two or more immiscible fluids are widely used in many industrial processes, for example, in casting operations, food processing, fibre coating, emulsification, spraying and so on. The kinematics and dynamics of the fluid interfaces play an important role in many of these processes. For example, the interface between the covering slag and the liquid steel in a continuous casting mould can strongly affect the quality of the final product, when interfacial instabilities are induced by the flow [1]. 

Fast, robust and accurate numerical models are necessary for the investigation of the interfacial fluid dynamics by means of computational fluid dynamics (CFD) on fixed computing grids. Prominent approaches are marker and cell [2], volume of fluid method (VOF) [3] and level set method (LS) [4, 5]. A detailed discussion of these approaches is beyond the scope of this paper; interested readers are referred to review papers, for example, by Rider and Kothe [6], Scardovelli and Zaleski [7] or Sethian and Smereka [8]. 

There are numerous papers available which demonstrate that marker and cell , VOF and LS are indeed fast and robust tools. However, as recently noticed by Hysing et al. [9], the accuracy of the numerical models is only seldom checked. Validation of two-phase flow simulations is often performed by qualitative comparison of interface locations in CFD simulations with pictures from corresponding experiments, see, for example, Ferziger and Peri [10] or Yeoh and Barber [11]. 

Therefore, Hysing et al. [9] have proposed a quantitative benchmark case, the two-dimensional rising bubble problem for two different set-ups. The benchmark is deduced from the numerical results of independent CFD simulations with three different codes. Two of them are based on the LS approach on fixed grids. The third one resolves the two-phase flow with an arbitrary Lagrangian– Eulerian approach where the computing mesh follows the phase interface. This is often referred to as interface tracking approach. The results of the three different codes are very close to each other, which suggests that the modelling errors should be very small. Hysing et al. [9] deduced measures from their results, which can be used to allow a quantitative comparison with results from other two-phase flow models. 

The benchmark of Hysing et al. [9] has been already employed by [12] to check the accuracy of a two-fluid approach. However, a quantitative comparison of the VOF model data with the LS and arbitrary Lagrangian–Eulerian results in Hysing et al. [9] is missing up to now. 

In the VOF approach, a phase transport equation is used to capture the interphase. There are basically two methods for evaluating phase transport equation: 

- direct methods: which discretise the phase transport equations directly by using high resolution discretisation schemes for the convection parts. Schemes specifically developed for interface phenomena such as, for example, compressive interface capturing scheme for arbitrary meshes [13] or high-resolution interface-capturing scheme [14]. There are no limitations for the underlying computing mesh for the direct methods. 

- reconstruction methods: where the phase transport equation is approximated typically in two steps, a geometric interface reconstruction step and an interface propagation step. Methods following this approach include donor–acceptor method [3], simple line interface calculation, piecewise linear interface construction [15] or efficient least-squares volume of fluid interface reconstruction algorithm [16]. All these methods are restricted to hexahedral (or quadratic in two dimensions) shaped mesh cells. The recently presented moments of fluids method by Dyadechko and Shashkov [17,18] overcomes this limit and extends the reconstruction methods to arbitrary cell shapes. However, these methods are more computational demanding. 

Özkan et al. [19] evaluated the two different VOF methods for a three-dimensional bubble-train flow. Their study reveals that there are some deficiencies for the evaluated direct methods (implemented in the commercial CFD-solvers STAR-CD, CFX and FLUENT). 

In this paper, we evaluate the VOF approach implemented in the CFD library OpenFOAM®version 1.5.1 through version 2.1.0 [20,21]. In OpenFOAM®, the phase transport equation is solved directly, which is described in detail in Section 2 

The rest of the paper is outlined as follows: A short summary of the test cases and the numerical set-up of the present investigation is given in Section 3. The results of the VOF model applied to the two test cases are presented, compared with the findings of Hysing and co-workers [9] and discussed in detail in Section 4. 

## 2. MATHEMATICAL MODEL AND COMPUTATIONAL METHODS

## 2.1. Governing equations

The rising bubble problem investigates the ascending motion of a bubble of fluid 2 that is surrounded by fluid 1. To track the two-phase flow of the rising bubble, the VOF method of Hirt and Nichols [3] is adopted. Here, a volume fraction $\gamma _ { 1 }$ is used to mark fluid 1 with the definition 

$$
\gamma_ {1} = \left\{ \begin{array}{l l} 1 & \text { in   fluid } 1 \\ 0 & \text { in   fluid } 2 \\ 0 <   \gamma_ {1} <   1 & \text { at   the   interface } \end{array} \right. \tag {1}
$$

The governing equations of the unsteady, laminar and incompressible flow are the equation of continuity (Equation (2)) and the Navier–Stokes equation (Equation (3)) 

$$
\nabla \cdot \vec {u} = 0 \tag {2}
$$

$$
\frac {\partial}{\partial t} (\rho \vec {u}) + \nabla \cdot (\rho \vec {u} \vec {u}) = - \nabla p + 2 \nabla \cdot \left[ \mu \stackrel {\leftrightarrow} {D} \right] + \rho \vec {g} + \sigma \kappa \nabla \gamma_ {1} \tag {3}
$$

Here, uE is the velocity, $p$ is the pressure and $\vec { g }$ is the gravitational acceleration. The mixed fluid properties density $\rho$ and viscosity $\mu$ are weighted by the volume fractions $\gamma _ { 1 }$ and $\gamma _ { 2 } = 1 - \gamma _ { 1 }$ of the two fluids 

$$
\rho = \gamma_ {1} \rho_ {1} + (1 - \gamma_ {1}) \rho_ {2}, \quad \mu = \gamma_ {1} \mu_ {1} + (1 - \gamma_ {1}) \mu_ {2} \tag {4}
$$

The rate of strain tensor is defined by 

$$
\overleftrightarrow {D} = \frac {1}{2} (\nabla \vec {u} + (\nabla \vec {u}) ^ {T}) \tag {5}
$$

Because of the continuous surface force (CSF) approach of Brackbill et al. [22], the surface tension force is included as an additional source term  $\cdot \nabla \gamma _ { 1 }$ in the Navier–Stokes equation. The curvature $\kappa$ is defined as 

$$
\kappa = - \nabla \cdot (\vec {n}) \tag {6}
$$

the divergence of the normal vector $\vec { n }$ of the interface. Its definition is given in Section 2.2.2. 

Solving the phase transport equation by the surface compression approach. The transport equation of each volume fraction $\gamma _ { 1 }$ and $\gamma _ { 2 }$ in a incompressible two-fluid system is given by 

$$
\frac {\partial \gamma_ {i}}{\partial t} + \nabla \cdot (\vec {u} _ {i} \gamma_ {i}) = 0, i = 1, 2 \tag {7}
$$

For the derivation of the surface compression approach, it is sufficient to consider the transport equation of the volume fraction $\gamma _ { 1 }$ only 

$$
\frac {\partial \gamma_ {1}}{\partial t} + \nabla \cdot (\vec {u} _ {1} \gamma_ {1}) = 0 \tag {8}
$$

To solve this transport equation, the velocity $\vec { u } _ { 1 }$ of fluid 1 is needed. In the widely used original VOF method by Hirt and Nichols [3], the velocity $\vec { u } _ { 1 }$ is assumed to be equal to the mixed velocity $\vec { u } = \vec { u } _ { 1 }$ 

$$
\frac {\partial \gamma_ {1}}{\partial t} + \nabla \cdot (\vec {u} \gamma_ {1}) = 0 \tag {9}
$$

This is only valid if $\gamma _ { 1 }$ is maintained as a step function throughout the domain, for example, numerical diffusion at the interface is not allowed. Furthermore, the upper boundedness $\gamma _ { 1 } \leqslant 1$ for Equation (9) is not guaranteed because this formulation is not conservative [23]. 

According to Rusche [23], it was Weller who first developed and implemented a conservative form in the CFD library OpenFOAM® [20, 21]. He defined the mixed uE and the relative velocities $\vec { u } _ { r }$ between phases $\gamma _ { 1 }$ and $\gamma _ { 2 }$ , 

$$
\vec {u} = \gamma_ {1} \vec {u} _ {1} + \gamma_ {2} \vec {u} _ {2} = \gamma_ {1} \vec {u} _ {1} + (1 - \gamma_ {1}) \vec {u} _ {2} \tag {10}
$$

$$
\vec {u} _ {r} = \vec {u} _ {1} - \vec {u} _ {2} \tag {11}
$$

Then, the addition of Equations (10) (multiplied by $\gamma _ { 1 } )$ and (11) (multiplied by $1 - \gamma _ { 1 } )$ yields 

$$
\gamma_ {1} \vec {u} _ {1} = \gamma_ {1} \vec {u} + \gamma_ {1} (1 - \gamma_ {1}) \vec {u} _ {r} \tag {12}
$$

Finally, the relative velocity formulation of the transport of $\gamma _ { 1 }$ in the surface compression approach is obtained by inserting Equation (12) in Equation (8) 

$$
\frac {\partial \gamma_ {1}}{\partial t} + \nabla \cdot (\vec {u} \gamma_ {1}) + \nabla \cdot [ \vec {u} _ {r} \gamma_ {1} (1 - \gamma_ {1}) ] = 0 \tag {13}
$$

with the explicitly fixed relative velocity $\vec { u } _ { r }$ , which is defined in Section 2.2.2. The surface compression term $\nabla \cdot [ \vec { u } _ { r } \gamma _ { 1 } ( 1 - \gamma _ { 1 } ) ]$ contributes only in the region of the interface $( 0 < \gamma _ { 1 } < 1 )$ and limits the smearing of the interface because of the compensation of the diffusive fluxes. The convective term in the surface compression approach (Equation (13)) differ from the ones used in the other direct VOF methods (e.g. compressive interface capturing scheme for arbitrary meshes [13] and [14]) where the third term on the left-hand side (LHS) of Equation (13) is neglected. 

## 2.2. Discretised model equations

Equations (2), (3) and (13) are integrated with the CFD library $\mathrm { O p e n F O A M ^ { \otimes } \ [ 2 0 , 2 1 ] . O p e n F O A M ^ { \otimes } }$ uses the finite volume method in cell-centred formulation to solve systems of partial differential equations on three-dimensional block-structured or unstructured meshes consisting of arbitrarily shaped convex cells. With Gauss’s theorem applied to the convective and diffusive ‡ terms and using the Euler implicit time scheme § of Equation (3), the following semi-discretised system of equations is derived 

$$
\begin{array}{l} \int_ {V _ {P}} \left(\frac {\rho}{\Delta t} (\vec {u} ^ {n} - \vec {u} ^ {o})\right) d V + \sum_ {f} \rho_ {f} \phi \vec {u} _ {f} ^ {n} = \sum_ {f} \mu_ {f} \vec {S} \cdot \nabla_ {f} ^ {\perp} \vec {u} ^ {n} + \left(\frac {1}{V _ {p}} \sum_ {f} \vec {u} ^ {n - 1} \vec {S}\right) \cdot \sum_ {f} \mu_ {f} ^ {n - 1} \vec {S} \\ - \int_ {V _ {P}} (\nabla p _ {m} - \vec {g} \cdot \vec {x} \nabla \rho + \sigma \kappa \nabla \gamma) d V \tag {14} \\ \end{array}
$$

where $\phi = \vec { u } _ { f } \cdot \vec { S }$ represents the volumetric flux through the cell face $S .$ The subscript $( ) _ { f }$ denotes face values, which are interpolated from the cell values indicated by $( ) _ { p }$ . The superscript $( ) ^ { o }$ indicates previous (old) time step. The superscript $( ) ^ { n }$ denotes the value of the actual iteration within the current time step, thus indicating implicit treatment. Respectively, the superscript $0 ^ { n - 1 }$ gives the value of the previous iteration, thus indicating explicit treatment. Note that $\phi$ is always a face value, thus the subscript is omitted here. For simplification, the subscript of the volume fraction is also omitted, thus $\gamma = \gamma _ { 1 }$ . The operator $\nabla _ { f } ^ { \perp }$ denotes a surface normal gradient, which is the component normal to the cell face of the gradient. This term is part of the Laplacian term and is discretised as follows $\begin{array} { r } { \nabla _ { f } ^ { \perp } \vec { u } ^ { n } = \frac { \vec { u } _ { N } ^ { n } - \vec { u } _ { P } ^ { n } } { | d | } } \end{array}$ . The index P denotes the cell of interest, the index N the neighbouring cell and $| \check { d } |$ the distance between both neighbouring cells. The pressure $p$ has been replaced by the modified pressure $p _ { m } = p - \rho \vec { g } \cdot \vec { x }$ with $\bar { \vec { x } }$ as the position vector. Its derivation is given by Rusche [23]. He states that this treatment enables an efficient numerical treatment of the steep density jump at the interface by including the hydrostatic term $\vec { g } \cdot \vec { x } \nabla \rho$ into the Rhie–Chow correction. 

## 2.2.1. Segregated pressure-based solver (velocity–pressure coupling). For the description of the velocity–pressure coupling, we follow the procedure of Rusche [23] and Medina [24].

Equations (2) and (3) are solved by the segregated pressure correction method, which will be outlined in this section. If we take only the terms containing the velocities uE from Equation (14), the following system of linear algebraic equations can be defined. 

$$
\mathcal {A} := \int_ {V _ {P}} \frac {\rho}{\Delta t} \left(\mathbf {U} ^ {n}\right) d V + \underbrace {\sum_ {f} \rho_ {f} \phi \mathbf {U} _ {f} ^ {n}} _ {\nabla \cdot (\rho \vec {u} \vec {u})} + \underbrace {\sum_ {f} \mu_ {f} \cdot \nabla_ {f} ^ {\perp} \mathbf {U} _ {f} ^ {n} + \sum_ {f} \mathbf {U} _ {f} ^ {n - 1} \cdot \sum_ {f} \mu_ {f} ^ {n - 1} | \vec {S} |} _ {\nabla \cdot (\mu \nabla \vec {u} + (\nabla \vec {u}) ^ {T})} \tag {15}
$$

This (A) represents the full system of equations $\mathbf { A U } = \mathbf { b } $ , with A being the momentum coefficient matrix, U the vector of the unknown velocity vectors uE and b the vector containing the source terms. The following operators acting on $\mathcal { A }$ can be defined and are implemented in $\mathrm { O p e n F O A M ^ { \circledast } } ,$ : the $ { \mathrm { ^ { 6 } D } } ^ {  { \prime }  { \tau } }$ -operator $A _ { \mathcal { D } }$ gives the diagonal part of $\mathbf { A } ,$ the $\mathbf { \bar { \Psi } } ^ { \bullet } \mathbf { N } ^ { \mathbf { \vec { \curlyeq } } }$ -operator $\mathcal { A } _ { \mathcal { N } }$ representing the lower and upper triangle (off-diagonal) part of $\mathbf { A } ,$ the $\mathbf { \bar { \Sigma } } ^ { 6 6 } \mathbf { S } ^ { \prime }$ -operator $\mathcal { A } _ { \mathcal { S } }$ extracting the source vector b and the $^ { \mathrm { \mathfrak { s } } } \mathrm { H } ^ { \mathrm { \mathfrak { s } } }$ -operator $\boldsymbol { A } _ { \mathcal { H } }$ defined by $\mathcal { A } _ { \mathcal { H } } : = \mathcal { A } _ { S } - \mathcal { A } _ { N } \mathbf { U } = \mathcal { A } _ { D } \mathbf { U }$ . With the operators given earlier, the equation system $\mathcal { A }$ can be rewritten as $( \mathcal { A } _ { \mathcal { D } } + \mathcal { A } _ { \mathcal { N } } ) \mathbf { U } \ = \ \mathcal { A } _ { \mathcal { S } }$ , which allows us to rewrite Equation (14) as 

$$
\mathcal {A} _ {\mathcal {D}} \mathbf {U} ^ {n} = \mathcal {A} _ {\mathcal {H}} + \int_ {V _ {P}} \left(- \nabla p _ {m} - \vec {g} \cdot \vec {x} \nabla \rho + \sigma \kappa \nabla \gamma + \frac {\rho}{\Delta t} \mathbf {U} ^ {o}\right) d V \tag {16}
$$

Assembling a block matrix system of Equations (16) and (2) 

$$
\left[ \begin{array}{c c} \mathcal {A} _ {\mathcal {D}} & \nabla \\ \nabla \cdot & 0 \end{array} \right] \left[ \begin{array}{l} \mathbf {U} ^ {n} \\ p _ {m} \end{array} \right] = \left[ \begin{array}{c} \mathcal {A} _ {\mathcal {H}} - \int_ {V _ {P}} \left(\nabla p _ {m} - \vec {g} \cdot \overrightarrow {z} \nabla \rho + \sigma \kappa \nabla \gamma + \frac {\rho}{\Delta t} \mathbf {U} ^ {o}\right) d V \\ 0 \end{array} \right] \tag {17}
$$

and taking the Schur complement of $\mathrm { i t } ^ { \ P }$ gives the pressure equation: 

$$
\nabla \cdot \left[ \left(\mathcal {A} _ {\mathcal {D}} ^ {- 1}\right) _ {f} \nabla [ p _ {m} ] \right] = \nabla \cdot \phi^ {*} \tag {18}
$$

with flux predictor $\phi ^ { * }$ 

$$
\begin{array}{l} \phi^ {*} = (\mathcal {A} _ {\mathcal {D}} ^ {- 1} \mathcal {A} _ {\mathcal {H}}) _ {f} \cdot \vec {S} + \underbrace {(\mathcal {A} _ {\mathcal {D}} ^ {- 1}) _ {f} \left[ - (\vec {g} \cdot \vec {x}) _ {f} | \vec {S} | \nabla_ {f} ^ {\perp} \rho + (\sigma \kappa) _ {f} | \vec {S} | \nabla_ {f} ^ {\perp} \gamma \right]} _ {\text { source   term }} \tag {19} \\ + \underbrace {(\mathcal {A} _ {\mathcal {D}} ^ {- 1}) _ {f} \left(\frac {\rho}{\Delta t} \mathbf {U} ^ {o}\right) _ {f} | \vec {S} |} _ {\text { extra   source   from   time   derivative }} \\ \end{array}
$$

The flux $\phi$ is corrected by 

$$
\phi = \phi^ {*} - (\mathcal {A} _ {\mathcal {D}} ^ {- 1}) _ {f} \vec {S} \cdot \nabla_ {f} ^ {\perp} p _ {m} \tag {20}
$$

after the calculation of the pressure with Equation (18). 

This method gives an oscillation-free velocity field in line with the Rhie–Chow correction [25], even though there is no explicit Rhie–Chow correction. 

To obey the conservation of the fluxes on the right-hand side of equations ((18), (19) and (20)), one has to assure consistent treatment of all terms because they are assembled explicit. 

2.2.2. Surface compression. In $\mathrm { O p e n F O A M ^ { \textregistered } }$ [20, 21], the transport equation for the volume fraction - is given by Equation (13). By default, a relative velocity $\vec { u } _ { r }$ cannot be determined in the VOF method because only a single velocity uE for both fluids is considered in the whole domain. Thus, $u _ { r }$ or rather the cell face flux $\phi _ { r }$ velocity has to be approximated as given in Equation (25). The semi-discretised form of Equation (13) is 

$$
[ \gamma ] _ {i} ^ {n + 1} = [ \gamma ] _ {i} ^ {n} - \Delta t \{\nabla \cdot [ (\gamma \vec {u}) + \vec {u} _ {r} \gamma (1 - \gamma) ] \} \tag {21}
$$

Gauss’s theorem and integration over a control volume $V _ { p }$ gives 

$$
\int_ {V _ {P}} \nabla \cdot [ (\gamma \vec {u}) + \vec {u} _ {r} \gamma (1 - \gamma) ] d V = \sum_ {f} \{(\gamma \phi) _ {f} + (\gamma \phi_ {r} (1 - \gamma)) _ {f} \} \tag {22}
$$

Introducing $\phi _ { \gamma }$ as 

$$
\phi_ {\gamma} := \{(\gamma \vec {u}) + [ \vec {u _ {r}} \gamma (1 - \gamma) ] \} _ {f} \cdot \vec {S} = (\gamma \phi) _ {f} + (\gamma \phi_ {r} (1 - \gamma)) _ {f} \tag {23}
$$

which results in the discretised transport equation for the phase indicator 

$$
\int_ {V _ {P}} [ \gamma ] _ {i} ^ {n + 1} d V = \int_ {V _ {P}} [ \gamma ] _ {i} ^ {n} d V - \Delta t \sum_ {f} \phi_ {\gamma} \tag {24}
$$

The compression flux, $\phi _ { r }$ in Equation (23), is defined by the magnitude of the face flux velocity in the transition region $\left( u _ { c } \right) _ { f }$ multiplied with the normal face vector of the interface $( \widehat { n } _ { i } ) _ { f }$ , which yields only the contribution normal to the interface. 

$$
\phi_ {r} := \left(u _ {c}\right) _ {f} \cdot \left(\widehat {n} _ {i}\right) _ {f} \tag {25}
$$

The face flux velocity in the transition region $( u _ { c } ) _ { f }$ is then evaluated by 

$$
\left(u _ {c}\right) _ {f} = \min \left(c _ {\gamma} \cdot \left| \frac {\phi}{| \vec {S} |} \right|, \left| \frac {\phi}{| \vec {S} |} \right| _ {\max}\right) \tag {26}
$$

which is bounded to the maximum face flux velocity $\left| \left. \frac { \phi } { | \vec { S } | } \right| _ { \mathrm { m a x } } \right|$ in the flow field. The coefficient $c _ { \gamma }$ controls the weight of the compression flux $\phi _ { r }$ and should be in the range of unity [24]. By choosing $c _ { \gamma } = 0$ , the compression flux can be forced to $\phi _ { r } = 0$ . The contribution to the interface, which is normal at the cell face, is derived by the scalar projection $( \widehat { n } _ { i } ) _ { f }$ of the normal vector of the interface pointing from fluid $2 \left( \gamma = 0 \right)$ to fluid 1 $( \gamma = 1 )$ onto the cell face vector ${ \vec { S } } .$ This term is determined by the following inner product: 

$$
(\widehat {n} _ {i}) _ {f} = \vec {n} \cdot \vec {S} = \frac {\nabla_ {f} ^ {\perp} \gamma}{| \nabla_ {f} ^ {\perp} \gamma + b |} \cdot \vec {S} \quad \text { with } \quad b = \frac {1 0 ^ {- 8}}{\sqrt [ 3 ]{\frac {\sum V _ {P}}{N}} N} \tag {27}
$$

1.set initial and boundary conditions for the fields of $\vec { u } , p , \gamma$ 

2.set time step, if time step is variable in accordance to the CFL criteria 

3.calculate viscosity μ and ρ (Eq. 4) with an explicit limiting of $0 \leq \gamma \leq 1$ 

4. perform γ -sub cycle $N _ { \gamma - s u b }$ -times (explicit sub time stepping ) 

(a）update/define $\phi _ { r }$ (Eq. 25) 

(b）perform $\gamma$ correction $N _ { \gamma - c o r r }$ -times 

i. update/define  (Eq. 23) 

ii. solve Eq.24 explicitly for γ with MULES 

ii. update/define $\bar { \rho \phi ^ { i } }$ (Eq.29) 

(c) calculate pΦ (Eq. 30) 

(d) calculate $( \widehat { n } _ { i } ) _ { f }$ (Eq.27) and the curvature k (Eq. 28) 

(e) determine density p (Eq. 4) 

5. update mass flux (Eq. 29) again 

6. construct ${ \mathcal { A } } ,$ equation (Eq. 15) 

7. if (mom $P r e d = t r u e )$ do a momentum prediction step with the old pressure, (Eq. 16) 

8．perform PISO loop max. $\boldsymbol { N _ { c o r r } + 1 - } \mathrm { t i m e s }$ 

(a) (re-) assemble $A _ { \mathcal { D } } , A _ { \mathcal { H } }$ coefficient and calculate $\mathbf { U } = { \mathcal { A } _ { D } } ^ { - 1 } \mathbf { \mathcal { A } } _ { \mathcal { H } }$ 

(b) optional momentum predictor (Eq. 16) with the old pressure 

(c） solve fux predictor(Eq.19) with U from previous step 

(d) solve Eq. 18 Nnonortho-corr + 1 -times (nonorthogonal corrector loop) 

i. correct velocities with the new pressure field $( \mathrm { E q } , 2 0 )$ 

ii. based on the pressure solution,assemble conservative face flux 

(e) update cell-centred velocity field U with assembled momentum coefficient A 

(f） update boundary conditions 

(g）check for convergence 

9.begin from (2) or finish calculation 

Algorithm 1. Explicit two-phase numerical solution procedure. 

to prevent division by zero with $V _ { P }$ denoting the cell volume and N the number of cells. The scalar $( \widehat { n } _ { i } ) _ { f }$ is also used for the determination of the curvature in Equation (6); here, the divergence term is rewritten with Gauss’s theorem: 

$$
\kappa = - \nabla \cdot (\vec {n}) = - \sum_ {f} (\widehat {n} _ {i}) _ {f} \tag {28}
$$

2.2.3. Solution algorithm. The overall solution cycle is given in Algorithm 1. Here, the distribution of the volume fraction field ( --subcyle at step 4) is calculated ahead of the well-known PISO algorithm [26] (step 8) that updates the velocity and pressure fields. 

The key element of the --subcyle is the solution of Equation (24). To guarantee a sharp and bounded solution, a numerical scheme designed for the multi-dimensional advection equation should be employed (e.g. [27, 28]). For this study, the explicit multidimensional universal limiter with explicit solution (MULES) implemented in $\mathrm { O p e n F O A M ^ { \textregistered } }$ [20, 21] is used to integrate the --equation Equation (24) explicitly. Therefore, Equations (23) and (24) are solved together with the mass flux $\rho _ { f } \phi ^ { i }$ 

$$
\rho_ {f} \phi^ {i} := \rho_ {f} \vec {u} _ {f} \cdot \vec {S} = (\gamma \rho_ {1}) _ {f} \phi + [ (1 - \gamma) \rho_ {2} ] _ {f} \phi = \phi_ {\gamma} (\rho_ {1} - \rho_ {2}) _ {f} + (\rho_ {2}) _ {f} \phi \tag {29}
$$

for a number of $N _ { \gamma - \mathrm { c o r r } }$ corrector steps. During this correction, $\phi _ { r }$ stays constant, see Algorithm 1. 

Finally, the mass flux is calculated according to Equation (24) as 

$$
\rho_ {f} \phi = \sum_ {i = 1} ^ {N _ {\gamma - \mathrm{sub}}} \frac {\delta t ^ {i}}{\Delta t} \left\{\rho_ {f} \phi^ {i} \right\} \tag {30}
$$

with $\delta t ^ { i } = \Delta t / N _ { \gamma - \mathrm { s u b } }$ . Equations (4), (6), (25) and (30) are iterated inside the so-called - -sub cycle. This sub cycle consists of solving the equations $N _ { \gamma - \mathrm { s u b } ^ { - } } \mathrm { t i m e s }$ (sub time-stepping; here, also $\phi _ { r }$ is updated, see in Algorithm 1) and then averaging the mass flux. 

During the PISO loop, the mass flux $\rho _ { f } \phi$ in the advection term is maintained constant. Once the pressure–velocity equations are solved, the volumetric flux $\phi$ is updated. The overall numerical procedure is sketched in Algorithm 1. 

2.2.4. Discretisation schemes and used parameters. The applied discretisation schemes and the parameters of the numerical model are summarised in Table I. For convenience, the corresponding terminology of $\mathrm { O p e n F O A M ^ { \textregistered } }$ is also given. The transient term in the momentum equation, Equation (14), is discretised by the Crank–Nicolson scheme [29]. Different interpolation schemes are applied to determine cell face values in the convective terms. For the convective part in the momentum equation, a total variation diminishing scheme with the flux limiter function $\psi ( r ) =$ maxŒ0, min.2  r, 1/ is used. The smoothness parameter r is defined by the ratio of consecutive gradients. The van Leer [30] scheme with the flux limiter function $\dot { \mathbf { \rho } } ( r ) = ( r + | r | ) / ( 1 + | r | )$ / is used for the discretisation of the $\nabla \cdot ( \vec { u } \gamma )$ term in Equation (13). Here, the flux is bounded between 0 and 1. For the second convective part $\nabla \cdot ( \vec { u } _ { r } \ : ( 1 - \gamma ) \gamma )$ in Equation (13), the so-called interfaceCompression scheme is used. Here, the limiter function reads 

$$
\psi (\phi_ {P}, \phi_ {N}) = \min
$$

$$
\left(\max \left\{1 - \max \left[ \sqrt {1 - 4 \cdot \phi_ {P} \cdot (1 - \phi_ {p})}, \sqrt {1 - 4 \cdot \phi_ {N} \cdot (1 - \phi_ {N})} \right], 0 \right\}, 1\right) \tag {31}
$$

which is bounded between 0 and 1, with the flux $\phi _ { P }$ at the evaluated cell and the flux $\phi _ { N }$ at the neighbour cell. In our study, we included the optional momentum predictor step that brought small but noticeable improvements in the results. This will be discussed later. 


Table I. Used discretisation and interpolation schemes, parameters of the numerical model.


<table><tr><td>Term</td><td>Discretisation scheme (OpenFOAM® terminology)</td><td>Method</td></tr><tr><td><eq>\frac{\partial}{\partial t}(\rho \vec{u})</eq></td><td>CrankNicolson 1</td><td>Crank–Nicolson scheme [29]</td></tr><tr><td><eq>\nabla \cdot (\rho \vec{u}\vec{u})</eq></td><td>Limited linearV 1</td><td>TVD scheme, limiter function see text</td></tr><tr><td><eq>\nabla \cdot (\vec{u}\gamma)</eq></td><td>vanLeer 01</td><td>TVD scheme, van Leer limiter see text</td></tr><tr><td><eq>\nabla \cdot (\vec{u}_{r}(1-\gamma)\gamma)</eq></td><td>InterfaceCompression</td><td>Bounded limited scheme, limiter see text</td></tr><tr><td><eq>\nabla \vec{u}, \nabla \gamma</eq></td><td>Linear</td><td>CDS</td></tr><tr><td><eq>\nabla \frac{1}{f}\chi^{**}</eq></td><td>Linear corrected</td><td>Surface normal gradient: CDS from the neighbouring cell values with explicit correction on nonorthogonal meshes [31]</td></tr><tr><td><eq>\nabla \cdot (\mu \nabla \phi_{r\gamma})</eq></td><td>Linear corrected</td><td>Face values (<eq>\mu</eq>) approximated by CDS, and the resulting surface normal gradient is evaluated as given earlier by CDS with nonorthogonal correction [31]</td></tr><tr><td>Term</td><td>Interpolation scheme</td><td>Method</td></tr><tr><td><eq>(\chi)_{f}</eq></td><td>Linear</td><td>Default interpolation scheme for getting face values from cell values, for example, <eq>\rho, \mu</eq> and CDS from the neighbouring cell values</td></tr><tr><td>Parameter</td><td>Value</td><td>Notes</td></tr><tr><td><eq>c_{\gamma}</eq></td><td>1</td><td>Interface compression valid values between 0 and 1</td></tr><tr><td><eq>N_{\gamma -corr}</eq></td><td>1</td><td>Number of <eq>\gamma</eq> corrector steps</td></tr><tr><td><eq>N_{\gamma -sub}</eq></td><td>2</td><td>Number of <eq>\gamma</eq>-sub cycles</td></tr><tr><td><eq>N_{\text{nonortho-corr}}</eq></td><td>0</td><td>Only used for nonorthogonal meshes</td></tr><tr><td><eq>N_{\text{corr}}</eq></td><td>15</td><td>Max. no. of PISO loops</td></tr><tr><td>momPred</td><td>true</td><td>Momentum prediction step</td></tr></table>


TVD, total variation diminishing; CDS, central differencing scheme. 



The symbol $\chi$ stands for an arbitrary scalar or vector. 


## 3. NUMERICAL SET-UP OF THE SIMULATIONS

## 3.1. Computational domain

On the basis of the benchmark definition of Hysing et al. [9], a two-dimensional computational domain with an aspect ratio $x : y = 1 : 2$ is employed, see Figure 1. The bubble is initially centred at $( x , y ) = ( 0 . 5 , 0 . 5 )$ with $r _ { b 0 } = 0 . 2 5$ as the initial radius. The viscosity and density of fluid 2 (bubble $\Omega _ { 2 } )$ are smaller than those of the surrounding fluid 1 $( \Omega _ { 1 } )$ . The domain is fully enclosed by no-slip walls $( \vec { u } = ( 0 , 0 ) )$ at the top and the bottom and free slip walls $( u _ { x } = 0$ and $\dot { \nabla } _ { f } ^ { \perp } \vec { u } = 0 )$ on the left and the right. The gravity vector $\vec { g }$ points towards the bottom of the domain. 

Hysing et al. [9] distinguished two different set-ups of the numerical experiment: ellipsoidal bubble (TC1) and skirted bubble (TC2). They differ by the density ratio $\rho _ { 1 } / \rho _ { 2 }$ and the ratio of the viscosity $\mu _ { 1 } / \mu _ { 2 }$ between both phases (index 1 for heavy and index 2 for light fluid). Physical parameters of TC1 and TC2 are given in Table II. 

All parameters in Table II are made dimensionless with characteristic scales for the length $L = 2 r _ { b 0 }$ , for the time $t = L / U _ { g }$ and for the rising velocity $U _ { g } = \sqrt { g 2 r _ { b 0 } }$ . The Reynolds number $R e$ , the Etvs number $E o$ and the capillary number Ca are defined as 

$$
R e = \frac {\rho_ {1} U _ {g} L}{\mu_ {1}}, \quad E o = \frac {\rho_ {1} U _ {g} ^ {2} L}{\sigma}, \quad C a = \frac {\mu_ {1} U _ {g}}{\sigma} = \frac {E o}{R e} \tag {32}
$$

Thus, the surface tension force is more prominent for TC1 in comparison with TC2, see Table II. 

The simulations of TC1 and TC2 are performed with four different grid spacings $\begin{array} { r l } { h } & { { } = } \end{array}$ $1 / [ 4 0 , 8 0 , 1 6 0 , 3 2 0 ]$ , respectively. The simulation time is $t _ { \mathrm { f i n a l } } = 3$ with a grid size-dependent time step of $\Delta t = h / 2$ . 

![](images/b7edbbf7dd2b29f5d9856c0061c8bd6f63d7e968301a49347d133b78cac6cebd.jpg)



Figure 1. Simulation domain, boundary conditions and initial configuration of the rising bubble problem due to Hysing et al. [9], length ratio $\mathrm { \dot { X } } { : } \mathrm { y } { = } 1 { : } 2$ , gravity $\vec { g } , \Omega _ { 1 }$ and $\bar { \Omega } _ { 2 }$ regions of heavy and light phase, respectively, in green free surface $( \gamma = 0 . 5 )$ , bubble centre at $( x , y ) = ( 0 . 5 , 0 . 5 )$ and initial bubble diameter $\bar { 2 } r _ { b 0 } = 0 . 5 .$ . The boundary conditions are slip (blue) and no-slip walls (red).



Table II. TC1 and TC2: Physical properties and similarity parameters.


<table><tr><td>Case</td><td><eq>\rho_1</eq></td><td><eq>\rho_2</eq></td><td><eq>\mu_1</eq></td><td><eq>\mu_2</eq></td><td>g</td><td>σ</td><td>Re</td><td>Eo</td><td>Ca</td><td><eq>\rho_1/\rho_2</eq></td><td><eq>\mu_1/\mu_2</eq></td></tr><tr><td>TC1</td><td>1000</td><td>100</td><td>10</td><td>1</td><td>0.98</td><td>24.5</td><td>35</td><td>10</td><td>0.286</td><td>10</td><td>10</td></tr><tr><td>TC2</td><td>1000</td><td>1</td><td>10</td><td>0.1</td><td>0.98</td><td>1.96</td><td>35</td><td>125</td><td>3.571</td><td>1000</td><td>100</td></tr></table>


Table III. TC1 and TC2: Simulation statistics, grid refinement $1 / h ,$ number of elements $n _ { \mathrm { e l e m e n t s } } ,$ total number of time steps $n _ { \Delta t }$ and computation time CPU on one processor (dual-core AMD Opteron 2.4 GHz).


<table><tr><td>1/h</td><td><eq>n_{\text{elements}}</eq></td><td><eq>n_{\Delta t}</eq></td><td>CPU TC1</td><td>CPU TC2</td></tr><tr><td>40</td><td>3 200</td><td>240</td><td>1 min</td><td>1 min</td></tr><tr><td>80</td><td>12 800</td><td>480</td><td>10 min</td><td>10 min</td></tr><tr><td>160</td><td>51 200</td><td>960</td><td>2.5 h</td><td>2.5 h</td></tr><tr><td>320</td><td>204 800</td><td>1 920</td><td>35 h</td><td>37 h</td></tr></table>

The resulting number of elements $n _ { \mathrm { e l e m e n t s } }$ and time steps $n _ { \Delta t }$ along with computing times CPU are given in Table III. The computing times are the measured wall clock times for serial simulations on a single 2.4 GHz dual-core AMD Opteron processor (AMD Corporation, Sunnyvale, California, USA). The computing times for both test cases differ only marginally and are comparable with the ones given by Hysing and co-workers [9] indicating that OpenFOAM® has a similar computing performance. 

## 3.2. Benchmark quantities

The following specific values of the benchmark quantities are used for the comparison of our simulation results with the reference data of Hysing et al. [9]: 

3.2.1. Final vertical position of the centre of mass $Y _ { c } ( t = 3 )$ . The final position of the centre of mass is obtained from the centre of mass position $\vec { X } _ { c } = ( X _ { c } , Y _ { c } )$ 

$$
\vec {X} _ {c} = \frac {\int_ {\Omega_ {1} \cap \Omega_ {2}} \gamma \vec {x} _ {c} d A}{\int_ {\Omega_ {1} \cap \Omega_ {2}} \gamma d A} \tag {33}
$$

Here, $\vec { x } _ { c } = ( x _ { c } , y _ { c } )$ is the centre of mass position of an individual cell in the computational mesh. 

3.2.2. Minimum circularity $C _ { m i n }$ with corresponding incidence time. The minimum circularity is deduced from the circularity 

$$
C = \frac {2 \pi r _ {b 0}}{P _ {b}} \tag {34}
$$

where $2 \pi r _ { b 0 }$ is the perimeter of the initial round bubble and $P _ { b }$ is the actual perimeter of the deformed bubble. The actual perimeter depends strongly on the position of the interface. Nevertheless, within the VOF method, the interface is defined somewhere in the $0 ~ < ~ \gamma ~ < ~ 1$ region. Theoretically, this region should tend to zero; in practice, the region is spread over a few cells so that the surface position can vary within this region. To pinpoint the surface, we chose a value of $\gamma = 0 . 5$ . 

3.2.3. Maximum rise velocity $V _ { m a x }$ with corresponding incidence time. The maximum rise velocity is deduced from the rise velocity $V ,$ which is evaluated for each time step as follows. 

$$
V = \frac {\int_ {\Omega_ {1} \cap \Omega_ {2}} \gamma v d A}{\int_ {\Omega_ {1} \cap \Omega_ {2}} \gamma d A} \tag {35}
$$

Here, v is the velocity in y-direction in an individual cell in the computational mesh. 

![](images/dfa0805f776c07efce9885c07838a663541bd94e441d3cf69f5f65d2c237dbbe.jpg)



(a) domain after patching


![](images/1f52b95ba33b563d2a909fdee7b5c7e027569bbec3e0a9767a1e4d4101e4099f.jpg)



(b) domain after relaxation



Figure 2. Distribution of both phases in the computational domain of the coarsest grid $( h = 1 / 4 0 )$ , yellow corresponds to $\gamma = 1$ , blue to $\gamma = 0$ and green to $0 ~ < ~ \gamma ~ < ~ 1$ , situation after patching and relaxation. (a) Domain after patching and (b) relaxation.


## 3.3. Initialisation

At the beginning of the simulations, the cell values for - are patched to either $\gamma = 1$ for the bubble or $\gamma = 0$ for the surrounding fluid. An example of the staircase profile resulting from the patching process can be seen in Figure $2 ( \mathrm { a } ) ;$ no cells with $0 < \gamma < 1$ have been defined. A transient simulation with surface tension but without gravity is employed to relax the unphysical staircase profile of the bubble surface. The results of $\gamma$ from the transient simulation with the zero gravity condition at $t = 3$ serve as initial values for the rising bubble simulation, see Figure 2(b). However, after the relaxation, the interface is smeared across four cells no matter which grid spacing was used. 

![](images/a662f1c98205feb6906277d2da4b679b10dabf7a9d0aec2e3cb1ceef537fbbae.jpg)



(a)


![](images/12b8b0706a5262470cb38b4f166c5e80a66a8e588d49a6bd3d9e4c5f51fe681f.jpg)



(b)


![](images/cd1a10d1713ab81a8705fd62d2680dfe3069460da19df0bb626939c12d406cfa.jpg)



(c)


![](images/58c8a1199ef7fe363741ebc8877b8904cd79e61da64bbda7060ceab50a5478d3.jpg)



(d)


![](images/590e26715c2b451519bb9fe14e787c6a5c4703be45b6f131562e9c1d89263da5.jpg)



Figure 3. The nondimensional velocity field during initialisation and zero gravity condition indicating the patterns of the parasitic currents for the conditions of TC1. The fields are nearly symmetric. For visualisation purposes, the pictures have been divided: the left sides show the flow patterns indicated by vector arrows, and the right sides show the coloured magnitude of the velocity.


## 4. RESULTS OF THE BENCHMARK TEST CASES

## 4.1. Zero gravity condition

After the initialisation with the zero gravity condition, some first errors and uncertainties of the numerical model can be estimated. Because of the action of the viscous forces, the momentum of the flow field should be damped out, and a final velocity of $\vec { u } = 0$ is expected in $\Omega _ { 1 }$ and $\Omega _ { 2 }$ . Thus, the bubble should keep its initial position $\vec { X } _ { c } ^ { i } = ( 0 . 5 , 0 . 5 )$ . However, small deviations from these assumed values are found in the zero gravity simulations at $t = 3$ . Spurious velocities are found on both sides of the interface, which we interpret as parasitic currents. These observations are in agreement with the findings in, for example, [12, 32–34] for a static viscous drop in equilibrium. The shape and the magnitude of the parasitic currents are shown in Figure 3 (for the parameters of TC1 at $t = 3 )$ . For all mesh resolutions, 16 counterrotating vortices can be identified in the flow field of these parasitic currents. These parasitic currents arise from small deviations in the curvature because of the implementation of the CSF model, see, for example, Lafaurie et al. [32]. 

In Figure 4, the influence of $c _ { \gamma }$ on the development of max $\left| u _ { p } \right|$ (local parasitic current maximum) over time is shown. Obviously, the parasitic currents are damped for both cases with and without surface compression until a small but finite value max $| u _ { p } | > 0$ . In detail, the development of $\operatorname* { m a x } \lvert u _ { p } \rvert$ shows small dependence on the grid spacing and on $c _ { \gamma }$ . In the case of $\phi _ { r } \neq 0 ( c _ { \gamma } = 1 )$ , a convergent behaviour for a resolution finer than $h \leqslant 1 / 4 0$ is noticeable. Note that the parasitic current is proportional to the kinetic energy in the system. To obey energy conservation, the kinetic energy and thus the parasitic currents should approach zero because of the inclusion of viscous forces. 

A snapshot of the pressure field for TC1 under zero gravity condition at $t = 3$ is depicted in Figure 5. The pressure is normalised by the pressure jump across the surface given by the Young– Laplace equation $\Delta p = 2 \sigma / r _ { b 0 }$ resulting in $p * = p / \Delta p = p r _ { b 0 } / ( 2 \sigma )$ . In the case of a static bubble, the normalised pressure should read in $\Omega _ { 1 } \colon p * = 0$ and $\Omega _ { 2 } \colon p * = 1$ with a sharp pressure jump at the surface. This is not the case because the values range from $p * = 0 . 8 3$ to $p * = 0 . 7 0$ (in $\Omega _ { 2 } )$ for the mesh resolution ranging from $h = 1 / 4 0 \mathrm { t o } h = 1 / 3 2 0$ , respectively. Furthermore, at the interface, some unphysical spikes (overshots and undershots) in the pressure field can be found. To reduce these spikes, we varied the interpolation scheme (least squares and cell-limited linear schemes) [20, 21] for the $\nabla \gamma$ term, which gave no improvements. 

The resulting errors of the zero gravity simulations for TC1 and TC2 on the coarsest and the finest grids are given in Table IV. The order of magnitude of $C a _ { p }$ remains nearly constant in both test cases. No grid convergence was found. Therefore, $C a _ { p } \sim 1 \bar { 0 } ^ { - 4 }$ is determined as a measure to estimate the influence of the parasitic currents in our simulations. 

However, the parasitic currents only have a marginal influence on the centre of mass position because this error is much smaller than the nominal displacement length of the parasitic currents, 

![](images/6e808e1f7e3c49144f100f3a5ad8da258e6ad14c47d45fb6d412d638014fa47b.jpg)



(a) ${ { \iota } _ { c _ { \gamma } } } = 0$


![](images/1cd53bb4a2a40b9f8c4ebccd877f2832f17d2a3ef84d0807a405728a4e26d977.jpg)



(b) $c _ { \gamma } = 1$



Figure 4. Influence of compression flux $\phi _ { r }$ on the parasitic currents. (a) $c _ { \gamma } = 0$ and (b) $c _ { \gamma } = 1$ .


![](images/87f20c9174f8c62ea8f29ac8155a094e1e2f3580f6e4a2231807ed54511be4a1.jpg)



(a) $h = 1 / 4 0$


![](images/f9a6280a0e2fb7ec948dea87c32c6b3a3411dc884be054b8bc55b00a4e5c30ca.jpg)



(b) $h = 1 / 8 0$


![](images/9fc2ecda02e37e0f092d0581050a8ab80814e00a1470fcf9eb546a9b8bb2f91b.jpg)



(c)


![](images/1fa33d55a0af234d05599ae69a6e18a3d51765ff1f2f3c315c900966874ab6da.jpg)



(d) $h = 1 / 3 2 0$


![](images/fc7bcf50b96f58961eff641537dfd2024aeff3f0927a22880f414357fdd299f0.jpg)



Figure 5. Pressure field after relaxation for initialisation of TC1. (a) $h = 1 / 4 0$ , (b) $h = 1 / 8 0 ,$ , (c) h D 1=160 and (d) $h = 1 / 3 2 0$ .



Table IV. The spatial errors in the zero gravity simulation at $t = 3$ for TC1 and TC2: centre of mass position $X _ { c } , Y _ { c }$ and difference $| \Delta \vec { X } _ { c } = \vec { X _ { c } } - \vec { X _ { c } ^ { i } } |$ with respect to the initial position $\vec { X } _ { c } ^ { i }$ , maximum of the parasitic current velocities max $| u _ { p } ( t = 3 ) |$ j and their capillary number $C a _ { p } = \eta | u _ { p } | _ { \operatorname* { m a x } } / \sigma$ .


<table><tr><td></td><td>1/h</td><td><eq>X_c</eq></td><td><eq>Y_c</eq></td><td><eq>|\Delta \vec{X}_c|</eq></td><td><eq>\max |u_p(t=3)|</eq></td><td><eq>Ca</eq></td></tr><tr><td rowspan="2">TC1</td><td>40</td><td>0.500011</td><td>0.499637</td><td><eq>3.6 \cdot 10^{-4}</eq></td><td>0.0112</td><td><eq>4.6 \cdot 10^{-4}</eq></td></tr><tr><td>320</td><td>0.500002</td><td>0.499948</td><td><eq>5.2 \cdot 10^{-4}</eq></td><td>0.00522</td><td><eq>2.1 \cdot 10^{-4}</eq></td></tr><tr><td rowspan="2">TC2</td><td>40</td><td>0.5</td><td>0.499776</td><td><eq>2.2 \cdot 10^{-4}</eq></td><td>0.00252</td><td><eq>1.3 \cdot 10^{-4}</eq></td></tr><tr><td>320</td><td>0.500007</td><td>0.499908</td><td><eq>9.0 \cdot 10^{-5}</eq></td><td>0.00575</td><td><eq>2.9 \cdot 10^{-4}</eq></td></tr></table>

![](images/0ec7b972f13038dc28c20760977d81accd667a93f6df05ee6feb308e12e77dc8.jpg)



(a) t = 0


![](images/15a89c89d8536f892badd63cc3e22d7a46a24bba2c267a73dc8af250e79f7c33.jpg)



(b) t = 1.5


![](images/df9d4921206a7b3be4ae887e29b3fbd88198b5cf41d788383fe2d0393781186c.jpg)



Figure 6. TC1: shapes of the rising bubble at different times. The bubble surface is indicated by a value of $\gamma = 0 . 5$ and fine grid h D 1=320. (a) t D 0, (b) t D 1.5 and (c) t D 3.


![](images/316dd87bdd2c53808831174b60c254d7b63579c70ba6eb08f4d48f9eae4dbf2c.jpg)



(a) overall situation


![](images/fe4b5f34e642209618d0fddd256ba3088671f3bd5876cc864a8ed33053536cd9.jpg)



(b) detail



Figure 7. TC1: final shape of the bubble at $t = 3 ,$ , volume of fluid simulations on the fine $h = 1 / 3 2 0$ (red) and the coarse grids $h \stackrel { \cdot } { = } 1 / 4 0$ (blue) and reference solutions by Hysing et al. [9] black (TP2D) and grey (FreeLIFE). (a) Overall situation and (b) detail.


$| \Delta \vec { X } _ { c } | \ll \operatorname* { m a x } | u _ { p } ( t = 3 ) | t$ . The effects of $| \Delta \vec { X } _ { c } |$ and max $| u _ { p } ( t = 3 ) |$ j on the results of the rising bubble simulations are discussed later. 

## 4.2. TC1

4.2.1. Simulation results. In Figure 6, typical results of the computed interface of the rising bubble for TC1 are shown. As illustrated for the stages $( t = 0 , t = 1 . 5$ and $t = 3 )$ , the bubble gets deformed, while staying compact as one connected region throughout the simulation. This is due to the low viscosity and density ratios and low Eo and Ca numbers indicating the high influence of the surface tension force. 

Figure 7 compares the bubble shape at the final time $\left( t \right. \ = \ 3 )$ in our simulations for the coarsest $( h = 1 / 4 0 )$ and the finest grids $( h = 1 / 3 2 0 )$ with the corresponding results for the codes 

![](images/76e45ef2b2e14fb04368cdfdfb9949e82a4873595b0aa479f03cf7bb55c89526.jpg)



(a) complete simulation period


![](images/b1ff84e2e6f0b4ac7cfc9c9375c047e43e151c6ae6351663ef29b990e8c63e91.jpg)



(b) detail with benchmark point $Y _ { c } ( t = 3 )$



Figure 8. TC1: temporal evolution of centre of mass position $Y _ { c } ,$ volume of fluid simulations on grids $h \stackrel { \sim } { = } 1 / [ 4 0 , 8 0 \AA$ , 160, 320 (coloured) and reference solutions from Hysing et al. [9] in black (TP2D) and grey (FreeLIFE). (a) Complete simulation period and (b) detail with benchmark point $Y _ { c } ( t = 3 )$ /.


![](images/5251284549ca33a99cc3d7463d008819f5e5d121443c4d51bd198a16ee8d5404.jpg)



(a) complete simulation period


![](images/b70ed29b095e731001a374fc7115a32d33f8c5b422faa35130ac7da8bc6a7305.jpg)



(b) detail with benchmark point $C _ { m i n } ( t _ { m i n } )$



Figure 9. TC1: temporal evolution of circularity $C ,$ volume of fluid simulations on grids $\begin{array} { r l } { h } & { { } = } \end{array}$ 1=Œ40, 80, 160, 320 (coloured) and reference solutions in black (TP2D) and grey (FreeLIFE). (a) Complete simulation period and (b) detail with benchmark point $\bar { C } _ { \mathrm { m i n } } ( t _ { \mathrm { m i n } } )$ .


TP2D/MooNMD and FreeLIFE|| given by Hysing et al. [9]. Note that the results of the TP2D and the MooNMD codes are indistinguishable for TC1; thus, we use only TP2D to refer both of them in this test case. Obviously, the VOF solutions on the coarsest and the finest grids are very close to each other, but there is a small but noticeable difference to the final bubble shape of the reference results, see Figure 7 (left). A detailed view at a typical section reveals that congruence between the VOF simulations and the reference results seems to be unreachable, not even with much finer grid resolutions. 

Figures 8–10 illustrate the temporal development of the quantities centre of mass position $Y _ { c } .$ , circularity C and rising velocity V from which the benchmark quantities are deduced. Again, the results of the VOF simulations for the four-grid resolutions h D 1=Œ40, 80, 160, 320 are compared with the corresponding results of the codes TP2D and FreeLIFE given by Hysing et al. [9]. 

The temporal evolution of $Y _ { c }$ in our simulations is nearly identical to the reference results, see Figure 8 (left), but the close-up showing the range $t = 2 . 5 - 3$ in Figure 8 (right) shows that complete convergence towards the benchmark results is again not achieved. On the contrary, the temporal evolution of C in Figure 9 (left) and V in Figure 10 (left) shows more apparent differences between our findings and the corresponding benchmark data. Here, the close-ups near the minimum circularity $C _ { \mathrm { m i n } }$ in Figure 9 (right) and the maximum rising velocity $V _ { \mathrm { m a x } }$ in Figure 10 (right) reveal that the results of our simulations do not show any convergence towards the reference. Instead of this result, our simulations scatter slightly around the results of the benchmark. 

![](images/7a588e5a3f0bf835699909d4426344ce8d1d8eed2a546a8750b291ccaf2aa3d6.jpg)



(a) complete simulation period


![](images/f8c3469822d84726e666255a881d2c790b7fbd7fbf31cece206779d05ca87b29.jpg)



(b) detail with benchmark point $V _ { m a x } \left( t _ { m a x } \right)$



Figure 10. TC1: temporal evolution of rising velocity $V ,$ volume of fluid simulations on grids $h \ =$ 1=Œ40, 80, 160, 320 (coloured) and reference solutions in black (TP2D) and grey (FreeLIFE). (a) Complete simulation period and (b) detail with benchmark point $\bar { V _ { \mathrm { m a x } } ( t _ { \mathrm { m a x } } ) }$ .



Table V. TC1: the benchmark quantities minimum circularity $C _ { \mathrm { m i n } }$ with corresponding incidence times $t ( C _ { \mathrm { m i n } } )$ , maximum rising velocity $V _ { \mathrm { m a x } }$ with $t ( V _ { \mathrm { m a x } } )$ and final position $Y _ { c }$ of centre of mass at $t \ = \ 3$ .


<table><tr><td>1/h</td><td>40</td><td>80</td><td>160</td><td>320</td><td>Benchmark Hysing et al. [9]</td><td>Štrubelj et al. [12]</td></tr><tr><td><eq>C_{\text{min}}</eq></td><td>0.8985</td><td>0.8966</td><td>0.8999</td><td>0.9044</td><td>0.9012 ± 0.0001</td><td>0.8876</td></tr><tr><td><eq>t(C_{\text{min}})</eq></td><td>1.8500</td><td>1.9625</td><td>1.9375</td><td>1.9625</td><td>1.9</td><td>1.8915</td></tr><tr><td><eq>V_{\text{max}}</eq></td><td>0.2364</td><td>0.2373</td><td>0.2365</td><td>0.2348</td><td>0.2419 ± 0.0002</td><td>0.2457</td></tr><tr><td><eq>t(V_{\text{max}})</eq></td><td>0.9250</td><td>0.9250</td><td>0.9219</td><td>0.9516</td><td>0.927 ± 0.011</td><td>0.9235</td></tr><tr><td><eq>Y_c(t=3)</eq></td><td>1.0616</td><td>1.0655</td><td>1.0668</td><td>1.0696</td><td>1.081 ± 0.001</td><td>1.0679</td></tr></table>

Table V summarises the benchmark quantities for TC1. Minimum circularity $C _ { \mathrm { m i n } }$ and maximum rising velocity $V _ { \mathrm { m a x } }$ , both with the corresponding incidence times $t ( C _ { \mathrm { m i n } } )$ and $t ( V _ { \mathrm { m a x } } )$ , respectively, and the final position $Y _ { c }$ of the centre of mass at $t = 3$ for TC1 are given. For convenience, the results of Štrubelj et al. [12] for TC1 are also included. The values of $h = 1 / 3 2 0$ serve as the basis. The maximum relative errors in our simulations for the grids $h = 1 / [ 4 0$ , 80, 160 are about 1% for $C _ { \mathrm { m i n } } , \ : V _ { \mathrm { m a x } }$ and $Y _ { c } ( t = 3 )$ , 3% for $t ( V _ { \mathrm { m a x } } )$ and 6% for $t ( C _ { \mathrm { m i n } } )$ . Comparing our results of grid $h = 1 / 3 2 0$ with the benchmark data of Hysing et al. [9], the relative errors are in the order of 1% for $C _ { \mathrm { m i n } }$ and $Y _ { c } ( t = 3 )$ and 3% for $V _ { \mathrm { m a x } } , t ( V _ { \mathrm { m a x } } )$ and $t ( C _ { \mathrm { m i n } } )$ . 

As explained in the previous section, the parasitic currents induce some uncertainties in the VOF simulations. From the findings of the zero gravity simulations for TC1, we estimate the magnitude of the parasitic currents to be in the order of u $\iota _ { P } \sim 0 . 0 1$ . Therefore, we expect a relative uncertainty in the flow field near the bubble interface in the order of $u _ { P } / | V | \sim 0 . 0 1 / 0 . 2 5 = 4 \%$ . Indeed, this is an upper limit of all the relative errors between our data and the benchmark results. So, we conclude that the performance of our VOF model for TC1 is mainly dominated by the existence of the parasitic currents. 

## 4.2.2. Influence of the model parameters.

Compression flux. For comparison, some simulations with zero compression flux $\phi _ { r } = 0$ have been carried out to reveal its influence on the sharpness of the surface and the temporal evolution of the circularity. 

In Figure 11, the influence of $\phi _ { r }$ on the surface is depicted. Two grid resolutions $( h \ =$ $1 / 3 2 0$ and $h = 1 / 4 0 )$ are investigated. For both grid resolutions, a clear smearing of the surface can be identified in the case of $\phi _ { r } = 0$ . In the case of $\phi _ { r } = 1$ , the smearing is clearly reduced for the coarse $( h = 1 / 4 0 )$ , and for the fine grid $( h = 1 / 3 2 0 )$ , it is invisible. 

![](images/5839271b57d73750d6c270640931b8634aefec7855f54754a9cb132e2e92db53.jpg)



(a) = 1/40


![](images/69b50a0b29d0da77e5fd84cf809aaa20cf2a2f4a2c6a476db5fed113136fe125.jpg)



(b) = 1/320



Figure 11. Influence of compression flux $\phi _ { r }$ on the surface, $\phi _ { r } ~ = ~ 0$ (left) and $\phi _ { r } ~ = ~ 1$ (right), for two different grid instances: $( \mathrm { a } ) h = 1 / 4 0$ and (b) $h = 1 / 3 2 0 .$ .


![](images/386f9ca65b190e9e7c2ead6407df11006aeb1f6f88c130ece772a9da47b236bb.jpg)



Figure 12. Influence of compression term on $C .$


The influence on the evolution of C shows exemplarily the influence of the compression flux, see Figure 12. In the case of the coarse mesh $( h = 1 / 4 0 )$ , the case of $\phi _ { r } = 0$ leads to a diverging behaviour at $t = 3$ . For the fine mesh case $( h = 1 / 3 2 0 )$ , a larger discrepancy to the reference solutions between $t = 1$ and $t = 2 . 5$ can be distinguished, whereas for $t \geqslant 2 . 5 .$ , the circularity converges towards the reference solutions. 

![](images/62a6e10219eb37343a8ab8b34838e1f6620c0dbb1f4aee80b2fe1b13fc07889c.jpg)



(a) Crank-Nicolson scheme


![](images/c6093777dfaf8654df795433848909672fc1aff335b8c628996aec57da080fd2.jpg)



(b) backward scheme



Figure 13. Influence of time step and time discretisation scheme on the evolution of C for the grid $h = 1 / 1 6 0 \AA$ . (a) Crank–Nicolson and (b) backward schemes.


![](images/d9ae159c4ff31f41a341110b20d22aff646ce628b267e72ff19cfb0979ff69d0.jpg)



Figure 14. Influence of the momentum predictor on the evolution of C .


Time schemes and time step. The time step $( d t = 1 / 1 6 0 , 1 / 3 2 0 , 1 / 6 4 0 )$ ) and the time discretisation scheme (Crank–Nicolson and backward) have been varied for TC1, with $h = 1 / 1 6 0$ to check for their influence. The temporal evolution of C was chosen exemplarily to show the negligible influence of these variations, see Figure 13. The choice of the time discretisation scheme had no influence on the evolution of C . 

Momentum predictor. In Figure 14, the influence of the momentum predictor step (see Section 2.2.3) on the temporal evolution of C for the case TC1, $h = 1 / 4 0$ and $1 / 3 2 0$ is depicted. In the low-resolution case $~ ( h ~ = ~ 1 / 4 0 )$ , influence of an additional momentum predictor step (momPred D 1) on the temporal evolution is negligible. For the finer resolution $( h \ = \ 1 / 3 2 0 )$ , we see that omitting the momentum predictor leads to larger deviations; hence, all simulations are carried out with the additional momentum predictor step. 

## 4.3. TC2

Figure 15 shows typical results for the rising bubble in TC2 at three different stages $( t = 0 , t =$ $1 . 5 \mathrm { a n d } t = 3 )$ . The distortion of the interface in TC2 is considerably larger than in TC1. The bubble is deformed to a skirted bubble and eventually forms thin filaments with smaller satellite bubbles that are detached. Compared with TC1, lower surface tension forces indicated by the higher Eo and Ca numbers as well as higher viscosity and density ratios are the reasons for the formation of a skirted bubble. The correct simulation of these physics is especially challenging as there is no unique numerical result for $t \geqslant 2$ available. For $Y _ { C }$ and V , the difference is less pronounced as for C . 

![](images/a568c461da5f3500b0c6a8f8249bbc947b5aa40a6ba2fd0aa03621aa178e4941.jpg)



(a) = 0


![](images/e5edece378fef703d9a1034c36b70c643f8ffe3ac84d918fa837b56378ce44d1.jpg)



(b) = 1.5


![](images/4183c5bf8925f76b9ba6ef5723f4d0245bc63a231b4d735a6de8afc685eabf4b.jpg)



(c) = 3



Figure 15. TC2: shapes of the rising bubble at different times, bubble surface indicated by a value of $\gamma = 0 . 5$ and simulation on grid $h = 1 / 3 2 0 .$ . (a) t D 0, (b) t D 1.5 and t D 3.


![](images/cf8397e8fdad9a04007da13f664a5c870bd7834f652d4a309f2b3e7662b909fc.jpg)



(a) overall situation


![](images/2e05739d44d6ab68af7b4f12ced97ff8dc361d9c4ba7518de45b209eea4edafd.jpg)



(b) in detail



Figure 16. TC2: final shape of the bubble at $t = 3 .$ , volume of fluid simulations on fine $h = 1 / 3 2 0$ (red) and coarse grids $h = 1 / 4 \hat { 0 }$ (blue) and reference solutions from Hysing et al. [9] in black (TP2D) and grey (FreeLIFE). (a) Overall situation and (b) in detail.


![](images/59c2b1c26394c3495042423831b4d6639107071fffa82b912708f957d052d3aa.jpg)



(a) complete simulation period


![](images/14cc1b3d858f152b13c22e5bff9a24307b96a792359ff09d858f597cd849a74b.jpg)



(b) detail with benchmark point $Y _ { c } ( t = 3 )$



Figure 17. TC2: temporal evolution of centre of mass position $\mathit { Y } _ { c } ,$ volume of fluid simulations on grids $h \doteq 1 / [ 4 0 , 8 0 , 1 6 0 , \dot { 3 } 2 0 ]$ (coloured) and reference solutions from Hysing et al. [9] in black (TP2D) and grey (FreeLIFE). (a) Complete simulation period and (b) detail with benchmark point $Y _ { c } ( t _ { 3 } )$ .


In Figure 16, the bubble shapes at the final time $\left( t \right. = \left. 3 \right)$ of our simulations for the coarsest $( h = 1 / 4 0 )$ and the finest grids $( h = 1 / 3 2 0 )$ with the corresponding results of the codes TP2D and FreeLIFE are compared. Note that we omit the results of MooNMD code because the applied arbitrary Lagrangian–Eulerian method does not allow for bubble breakup [9]. However, our solutions on the finest grid are very close to the final bubble shape of the reference calculations, see Figure 16 (left). The close-up at a range between $x \ = \ 0 . 6 5$ and $x = 0 . 7 5$ reveals that congruence between our VOF simulations and the reference result of FreeLIFE is achieved at this grid resolution, see Figure 16 (right). 

![](images/a76228a70288711a6f9fd0e336518b66d6856faeffbf298cbf1f2d8e5714c921.jpg)



(a) complete simulation period


![](images/8b0250bb3b576903ed184c8a8a7bbc588992e05981d10a4dcb5c300b87016bae.jpg)



(b) detail with benchmark point $C _ { m i n } ( t _ { m i n } )$



Figure 18. TC2: temporal evolution of circularity C , volume of fluid simulations on grids $\begin{array} { r l } { h } & { { } = } \end{array}$ 1=Œ40, 80, 160, 320 (coloured) and reference solutions in black (TP2D) and grey (FreeLIFE). (a) Complete simulation period and (b) detail with benchmark point $\breve { C } _ { \mathrm { m i n } } \not ( t _ { \mathrm { m i n } } )$ .


![](images/ee04e2764bef2e369723c7125056f56cfc2c4d50b1703e507987ee301e2a1b3a.jpg)



(a) complete simulation period


![](images/688de74ecc46ab87d94d4b333ff25d70be288203c2f75752705c3c3206112d24.jpg)



(b) with benchmark points $V _ { m a x 1 }$ $\left( t _ { m a x 1 } \right)$ $V _ { m a x 2 } \left( t _ { m a x 2 } \right)$



Figure 19. TC2: temporal evolution of rising velocity V , volume of fluid simulations on grids $h \ =$ 1=Œ40, 80, 160, 320 (coloured) and reference solutions in black (TP2D) and grey (FreeLIFE). (a) Complete simulation period and (b) with benchmark points $V _ { \operatorname* { m a x } 1 } ( t _ { \operatorname* { m a x } 1 } ) , \bar { V } _ { \operatorname* { m a x } 2 } ( t _ { \operatorname* { m a x } 2 } )$ .


The temporal development of the quantities $Y _ { c } , C$ and V is shown in Figures 17–19. 

Similar to TC1, the evolution of $Y _ { c }$ in our simulations is nearly identical to the reference results, see Figure 17 (left). The close-up of Figure 17 (right) confirms the good agreement between our results on grid $h = 1 / 3 2 0$ and the reference results; even so, convergence, that is a grid-independent solution, is again not achieved. 

Contrary to the findings for TC1, the evolution of C in Figure 18 (left) and V in Figure 19 (left) of our simulations are also closer to the corresponding results of FreeLIFE. Please note the differences in the benchmark results: in TP2D, there is a separation of satellite bubbles, whereas in FreeLIFE, the bubble remains skirted with elongated filaments. This difference between TP2D and FreeLIFE is documented by the very different temporal evolution of C . 

In our simulations, the rising bubble also remains skirted, although the separation of the satellite bubbles is found close to the final situation at $t = 3$ . The close-ups near the minimum circularity $C _ { \mathrm { m i n } }$ in Figure 18 (right) and the maximum rising velocities $V _ { \mathrm { m a x 1 } }$ and $V _ { \mathrm { m a x } 2 }$ in Figure 19 (right) confirm that the results of our simulations approximate the results of FreeLIFE, although complete convergence is again not achieved. 


Table VI. TC2: the benchmark quantities minimum circularity $C _ { \mathrm { m i n } }$ with corresponding incidence times $t ( C _ { \mathrm { m i n } } )$ , maximum rising velocities $V _ { \mathrm { m a x 1 / 2 } }$ with $t ( V _ { \mathrm { m a x } 1 / 2 } )$ and final position $Y _ { c }$ of centre of mass at $t = 3$ .


<table><tr><td>1/h</td><td>40</td><td>80</td><td>160</td><td>320</td><td>Benchmark(FreeLIFE) Hysing et al. [9]</td></tr><tr><td><eq>C_{\text{min}}</eq></td><td>0.5000</td><td>0.4999</td><td>0.5051</td><td>0.4945</td><td>0.4647</td></tr><tr><td><eq>t(C_{\text{min}})</eq></td><td>3.0000</td><td>2.8125</td><td>2.9875</td><td>2.9500</td><td>3.0000</td></tr><tr><td><eq>V_{\text{max1}}</eq></td><td>0.2388</td><td>0.2441</td><td>0.2431</td><td>0.2474</td><td>0.2514</td></tr><tr><td><eq>t(V_{\text{max1}})</eq></td><td>0.6875</td><td>0.7000</td><td>0.7250</td><td>0.7156</td><td>0.7281</td></tr><tr><td><eq>V_{\text{max2}}</eq></td><td>0.2134</td><td>0.2240</td><td>0.2302</td><td>0.2353</td><td>0.2440</td></tr><tr><td><eq>t(V_{\text{max2}})</eq></td><td>1.7000</td><td>1.8688</td><td>1.9594</td><td>2.0047</td><td>1.9844</td></tr><tr><td><eq>y_c(t=3)</eq></td><td>1.0862</td><td>1.1014</td><td>1.1089</td><td>1.1218</td><td>1.1249</td></tr></table>

Table VI summarises the benchmark quantities minimum circularity $C _ { \mathrm { m i n } }$ and maximum rising velocities $V _ { \operatorname* { m a x } 1 } , ~ V _ { \operatorname* { m a x } 2 }$ , with the corresponding incidence times $t ( C _ { \mathrm { m i n } } ) , t ( V _ { \mathrm { m a x 1 } } )$ and $t ( V _ { \operatorname* { m a x } 2 } )$ , respectively, and the final position $Y _ { c }$ of the centre of mass at $t = 3$ for TC2. The relative errors of our simulations decrease for an increasing grid resolution** up to the grid resolution of $h = 1 / 1 6 0$ , if $h = 1 / 3 2 0$ is taken as a basis. Nevertheless, complete convergence has not been achieved with the finest grid. Comparing our results of grid $h = 1 / 3 2 0$ with the benchmark data of FreeLIFE in [9], the relative errors are below 1% for $Y _ { c } ( t = 3 )$ , around 2% for $t ( C _ { \mathrm { m i n } } ) , V _ { \mathrm { m a x } 1 } , t ( V _ { \mathrm { m a x } 1 } )$ and $t ( V _ { \operatorname* { m a x } 2 } )$ , around 4% for $V _ { \mathrm { m a x } 2 }$ and around 6% for $C _ { \mathrm { m i n } }$ . 

From the findings of the zero gravity simulations for TC2, we estimate the magnitude of the parasitic currents to be in the order of u $_ P \sim 0 . 0 0 5$ . Therefore, we expect a relative uncertainty in the flow field near the bubble interface in the order of $u _ { P } / V | \sim 0 . 0 0 5 / 0 . 2 5 = 2 \%$ . We see that this value is again an upper limit of most of the relative errors between our data and the benchmark results. Similar to TC1, we conclude that the performance of our VOF model for TC2 is also dominated by the existence of the parasitic currents. 

## 5. CONCLUSIONS

We present a numerical model for the rising bubble problem. The numerical model is based on the volume of fluid method. The CSF approach is employed to describe the influence of the surface tension. Surface compression is used to sharpen the resolution of the bubble interface. The model equations are solved with the open-source CFD library OpenFOAM®. The implementation of the surface compression algorithm into $\mathrm { O p e n F O A M ^ { \textregistered } }$ is described in detail. 

Two different realisations of the rising bubble problem are investigated with the numerical model with four different grid resolutions. Both test cases differ in the Reynolds number $R e ,$ , the Etvs number $E o ,$ , the density and the viscosity ratios. Therefore, the final shape of the rising bubbles is ellipsoidal in one and skirted in the other test case. These two test cases have been also investigated by several other groups in the past. Their data for the circularity, the centre of mass and the rising velocity of the rising bubbles, serve as benchmark values for the VOF simulations. 

The comparison of our results with the benchmark data shows that the proposed model can resolve the bubble behaviour in reasonable agreement with the benchmark. In the case of the ellipsoidal bubble (TC1), our results deviate from the reference results of the benchmark data. It is very important to note that TC1 results do not show any grid convergence. This is especially critical because a grid-converged solution is commonly used if there is no benchmark (or experimental) result available. 

In the case of the skirted bubble (TC2), physics and the associated bubble shape are more complex. Here, correct simulation of the formation of the shape and detachment of satellite bubbles is especially critical and is still a challenge. As seen in the paper by Hysing et al. [9], a unique benchmark solution does not exist for TC2. For this test case, we were able to show grid convergence competing with the results given by [9]. 

In both cases, parasitic currents are observed in the flow field. These currents result from the numerical realisation of the CSF approach. We show that the magnitude of the parasitic currents defines an upper measure for the uncertainty of the VOF simulations. The presented VOF approach shows deficiencies with cases that are surface force dominated $( C a < 1 )$ as was the case for TC1. The study of Özkan et al.[19] indicates that such problems occur especially for VOF methods that solve the phase transport equations directly (direct methods). However, we did not experience these problems, such as the mentioned solver divergence problems. Although in our study, the capillary number is one order of magnitude higher, the deviations are much smaller compared with their study. We conclude that the surface compression approach is a more improved method compared with high-resolution interface-capturing scheme approach investigated by Özkan et al. [19]. 



6. Rider WJ, tracking K, DB. Reconstructing volume. Journal of Computational Physics 1998; 141(2):112–152. 



The advantages of the presented VOF approachare are 



7. Scardovelli R, Zaleski S. Direct numerical simulation of free-surface and interfacial flow. Annual Review of Fluid Mechanics Jan 1999; 31(1):567–603. 





8. Sethian J, Smereka P. Level set methods for fluid interfaces. Annual Review of Fluid Mechanics 2003; 35:341–372. 



1. boundedness of the volume fraction $0 \leqslant \gamma \leqslant 1$ , 



9. Hysing S, Turek S, Kuzmin D, Parolini N, Burman E, Ganesan S, Tobiska L. Quantitative benchmark computations of two-dimensional bubble dynamics. International Journal for Numerical Methods in Fluids 2008; 60:1259–1288. 



2. mass conservation without any reinitialisation step (as in LSs), 



10. Ferziger JH, Peric M. ´ Computational Methods for Fluid Dynamics. Springer: Berlin, Germany, 2002. 



3. maintenance of sharp interfaces through the surface compression approach, 



11. Yeoh GH, Barber T. Assessment of interface capturing methods in computational fluid dynamics (CFD) codes—a case study. The Journal of Computational Multiphase Flows Jun 2009; 1(2):201–215. 



4. the surface compression approach gives better results compared to the ones without surface compression. 



12. Štrubelj L, Tiselj I, Mavko B. Simulations of free surface flows with implementation of surface tension and interface sharpening in the two-fluid model. International Journal of Heat and Fluid Flow 2009; 30:741–750. 



Nevertheless, on the basis of the results presented here, we propose that the used surface compression approach should be further improved in future versions of the VOF solver. It remains an open question whether the source of these deviations for surface-tension-dominated flows comes from the direct VOF method, from the CSF method or from the interaction of both methods. For further improvements, we suggest the following: 



13. Ubbink O. Numerical prediction of two fluid systems with sharp interfaces. PhD Thesis, Imperial College, University of London, 1997. 





14. Muzaferija S, Peric M. Computation of free-surface flows using the finite-volume method and moving grids. Numerical Heat Transfer Part B-Fundamentals December 1997; 32(4):369–384. 



1. The implementation of the CSF approach should be refined to minimise the parasitic currents. Or alternative methods, such as the continuous surface stress method, should be tested. 



15. Youngs D. Time-dependent Multi-material flow with Large Fluid Distortion. Academic Press Inc: (London) Ltd, 1982. http://www.scopus.com/inward/record.url?eid=2-s2.0-0019936463&partnerID=40&md5= f1dbae4707c3493b06e31cc4da893bb7, cited By (since 1996) 19. 



2. The improvement of the determination of the curvature (approximation of the interface position) and the correct determination of the pressure jump needs to be addressed. Here, a reconstruction method, for example, moments of fluid method [17, 18] could give some improvements. 



16. Popinet S. An accurate adaptive solver for surface-tension-driven interfacial flows. Journal of Computational Physics 2009; 228(16):5838–5866. http://www.scopus.com/inward/record.url?eid=2-s2.0-67649472402&partnerID= 40&md5=6fcd38f15413e444be2e8a6c13465c45, cited By (since 1996) 35. 



The results presented in this paper should raise the awareness to use the direct VOF method, currently implemented in OpenFOAM®version 1.5.1 through version 2.1.0 [20, 21], cautiously for surface-tension-dominated flows (e.g. free surface liquid metal flows). Here, we recommend careful validation with appropriate reference solutions (results from experiments and published benchmark results). 



17. Dyadechko V, Shashkov M. Moment-of-fluid interface reconstruction. Technical Report, Los Alamos National Laboratoty, P.O. Box 1663, Los Alamos, NM 87545, 2005. 



## ACKNOWLEDGEMENTS



18. Dyadechko V, Shashkov M. Reconstruction of multi-material interfaces from moment data. Journal of Computational Physics 2008; 227(11):5361–5384. http://www.sciencedirect.com/science/article/pii/S0021999107005748. 



This work presents a part of the results obtained under contract number SFB 799 C1 sponsored by the Deutsche Forschungsgemeinschaft (DFG), Germany. 



19. Özkan F, Wörner M, Wenka A, Soyhan H. Critical evaluation of CFD codes for interfacial simulation of bubbletrain flow in a narrow channel. International Journal for Numerical Methods in Fluids 2007; 55(6):537–564. DOI: 10.1002/fld.1468. http://www.scopus.com/inward/record.url?eid=2-s2.0-35248859711&partnerID=40&md5= b1e4b01a78a0402c575caf2d0f551eae, cited By (since 1996) 10. 



## REFERENCES



20. Jasak H, Gschaider B, Nilson H, Beaudoin M, et al. Openfoam-1.5-dev extend project on sourceforge. http: //openfoam-extend.svn.sourceforge.net/viewvc/openfoam-extend/trunk/Core/OpenFOAM-1.5-dev/;2009. 





21. OpenCFD Limited. User guide openfoam 1.6. Technical Report. http://www.opencfd.co.uk/openfoam;2009. 





1. Maiwald A, Scheller P, Brücker C, Schwarze R. Flow-induced emulsification of mold powder slag into liquid steel. STEELSIM2009 3rd International Conference Simulation and Modelling of Metallurgical Processes in Steelmaking, Leoben, 8.-10. September, 2009; 162–167. 





22. Brackbill JU, Kothe DB, Zemach C. A continuum method for modeling surface tension. Journal of Computational Physics 1992; 100(2):335–354. 





2. Harlow F, Welch J. Numerical calculation of time-dependent viscous incompressible flow of fluid with free surface. Physics of Fluids 1965; 8(12):2182–2189. 





23. Rusche H. Computational fluid dynamics of dispersed two-phase flows at high phase fractions. PhD Thesis, Imperial College of Science, Technology and Medicine, London, England, 2002. 





3. Hirt CW, Nichols BD. Volume of fluid (VOF) method for the dynamics of free boundaries. Journal of Computational Physics 1981; 39:201–225. 





24. de M, P B R. Study and numerical simulation of sediment transport in free-surface flow. PhD Thesis, University of Málaga, July 2008. 





4. Osher S, Sethian JA. Fronts propagating with curvature-dependent speed—algorithms based on Hamilton–Jacobi formulations. Journal of Computational Physics 1988; 79(1):12–49. 





25. Rhie CM, Chow LW. Numerical study of the turbulent flow past an airfoil with trailing edge separation. AIAA 1983; 21:1525–1532. http://ci.nii.ac.jp/naid/80001467460/en/. 





5. Sussman M, Smereka P, Osher S. A level set approach for computing solutions to incompressible two-phase flow. Journal of Computational Physics 1994; 114(1):146–159. 





26. Issa RI. Solution of the implicitly discretised fluid flow equations by operator-splitting. Journal of Computational Physics 1986; 62(1):40–65. 





27. LeVeque RL. Finite Volume Methodes for Hyperbolic Problems. Cambridge University Press: Cambridge, UK, 2002. 





28. Hirsch C. Numerical Computation of Internal and External Flows: The Fundamentals of Computational Fluid Dynamics, 2nd ed. Butterworth-Heinemann: Oxford, UK, 2007. 





29. Crank J, Nicolson P. A practical method for numerical evaluation of solutions of partial differential equations of the heat-conduction type. Advances in Computational Mathematics December 1996; 6(1):207–226. http://dx.doi.org/10. 1007/BF02127704. 





30. van L, B. Towards the ultimate conservative difference scheme. IV. monotonicity andconservation combined in a second order scheme. Journal of Computational Physics 1974; 14:361–370. 





31. Jasak H. Error analysis and estimation for the finite volume method with applications to fluid flows. PhD Thesis, Imperial College, University of London, 1996. 





32. Lafaurie B, Nardone C, Scardovelli R, Zaleski S, Zanetti G. Modelling merging and fragmentation in multiphase flows with surfer. Journal of Computational Physics 1994; 113(113):134–147. 





33. Harvie D, Davidson M, Rudman M. An analysis of parasitic current generation in volume of fluid simulations. Applied Mathematics Model October 2006; 30(10):1056–1066. http://www.sciencedirect.com/science/article/ B6TYC-4HC6M40-2/2/aba2563fc39451c85a3adf331f9c1816. 





34. Sussman M, Smith KM, Hussaini MY, Ohta M, Zhi-Wei R. A sharp interface method for incompressible two-phase flows. Journal of Computational Physics 2007; 221(2):469–505. 

