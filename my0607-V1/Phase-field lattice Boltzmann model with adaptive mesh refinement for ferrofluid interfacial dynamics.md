
# Phase-field lattice Boltzmann model with adaptive mesh refinement for ferrofluid interfacial dynamics

**Cite as:** Phys. Fluids 37, 022148 (2025); doi: 10.1063/5.0256574

**Submitted:** 6 January 2025 · **Accepted:** 28 January 2025 · **Published Online:** 25 February 2025

---

**Authors:**

Zhenchao Guo (郭振超),¹ Shiting Zhang (章诗婷),¹ Yuqi Zhu (朱玉麒),¹ Yang Hu (胡洋),¹ᵃ Qiang He (何强),² Xiaolong Yang (杨小龙),³ and Decai Li (李德才)²

**Affiliations:**

1. Beijing key laboratory of Flow and Heat Transfer of Phase Changing in Micro and Small Scale, School of Mechanical, Electronic and Control engineering, Beijing Jiaotong University, Beijing 100044, People's Republic of China
2. State Key Laboratory of Tribology in Advanced Equipment, Tsinghua University, Beijing 100084, People's Republic of China
3. School of Mechanical and Automotive Engineering, Guangxi University of Science and Technology, Liuzhou 545006, People's Republic of China

ᵃAuthor to whom correspondence should be addressed: yanghu@bjtu.edu.cn

---

## ABSTRACT

In this paper, we propose a phase-field model that integrates the lattice Boltzmann method with an adaptive mesh refinement technique to study the interfacial dynamics of ferrofluids. In this model, we employ the second-order conservative Allen–Cahn equation to accurately capture the ferrofluid interface. The velocity-based hydrodynamic equations and a magnetic scalar potential equation with a pseudo-time term are utilized to describe the flow and magnetic fields. All governing equations are solved using a finite difference lattice Boltzmann scheme. To effectively resolve the interfacial dynamics of ferrofluids while reducing computational overhead, the numerical scheme is implemented on a block-structured adaptive mesh. To evaluate the accuracy and efficiency of the proposed model, we conduct simulations on several benchmark problems, including a circular cylinder in a uniform magnetic field, the deformation of a ferrofluid droplet, and the rising of a bubble in ferrofluid. The results obtained show good agreement with exact solutions and well-validated results in the existing literature. Furthermore, three types of ferrofluid instabilities under a uniform magnetic field—namely, the Rosensweig instability, the Rayleigh–Taylor instability, and the Kelvin–Helmholtz instability—are also investigated. Numerical results demonstrate that the magnetic field can significantly promote or suppress the occurrence of flow instabilities.

Published under an exclusive license by AIP Publishing. https://doi.org/10.1063/5.0256574

---

## I. INTRODUCTION

Multiphase phase ferrofluid flows have garnered significant attention due to their inherently rich physical properties and diverse applications across various fields.[1–4] The presence of a magnetic field induces a jump in magnetic stress force across the ferrofluid interface, and this jump, combined with interfacial tension, leads to complex interfacial behaviors. Numerous theoretical and experimental studies have been conducted to investigate the dynamics of ferrofluid interfaces.[3,4] However, theoretical methods are often applicable only in simplified scenarios. Additionally, the opacity of ferrofluids presents challenges for experimental techniques in accurately measuring microscopic phenomena and transient behaviors. With advancements in computer hardware and numerical techniques, numerical simulations have emerged as powerful tools for gaining insights into ferrofluid interfacial dynamics with enhanced precision and efficiency. Various numerical methods, including the finite volume method,[5] finite difference method,[6] and finite element method,[7] have been employed to simulate fluid flows. Among these, the lattice Boltzmann method (LBM) has gained popularity as an effective computational tool. LBM is now widely used for solving complex flow systems, such as multiphase and multicomponent flows,[8–10] turbulent flows,[11,12] micro-flows,[13,14] fluid–solid interactions,[15,16] porous media flows,[17,18] and thermal flows.[19–21]

According to the mesh structure used, the lattice Boltzmann method (LBM) can be categorized into two types: uniform mesh-based LBM and non-uniform mesh-based LBM.[22] The former employs a consistent mesh size throughout the computational domain, simplifying the implementation process. In recent years, numerous researchers have focused on investigating two-phase ferrofluid flows using uniform mesh-based LBM.[23–29] Hu et al.[23] were the first to extend the phase field-lattice Boltzmann model to simulate the dynamics of multiphase ferrofluids. Their work successfully simulated several classical numerical examples, including the deformation of ferrofluid droplets under a uniform magnetic field, the coalescence of bubbles in ferrofluid systems, and the motion and merging of ferrofluid droplets on a flat surface in the presence of a permanent magnet. The research by Hu et al. has sparked a significant amount of subsequent work. For instance, Li et al.[24] developed a multiphase LBM that integrates magnetic fields to simulate the dynamics of multiphase ferrofluid systems, particularly focusing on droplet merging processes. In another significant contribution, Li et al.[25] introduced a lattice Boltzmann flux solver equipped with a self-correcting procedure to explore bubble dynamics in ferrofluids under a uniform magnetic field. Khan et al.[26] employed a simplified multiphase LBM in conjunction with a magnetic field correction method to investigate the dynamics of ferrofluid droplets on solid substrates with varying wettability, analyzing the wetting kinetics involved. Zhang et al.[27] applied a generalized conservative phase-field simplified multiphase lattice Boltzmann model to elucidate the dynamical mechanisms and a general deformation law of a ferrofluid droplet suspended between air and a liquid substrate under the influence of an applied vertical uniform magnetic field. Huang et al.[28] developed an enhanced multicomponent multiphase pseudopotential LBM coupled with a magnetic field solver to study the wetting dynamics of a ferrofluid droplet influenced by non-uniform magnetic fields and gravitational effects. He et al.[29] introduced a hybrid phase-field-LBM for multiphase ferrofluid flow, employing the immersed interface method to solve the Laplace equation for magnetic potential with an interface jump condition. It should be pointed out that, in multiphase ferrofluid problems, the region of interest predominantly lies near the interface, where surface tension, magnetic interfacial forces, and fluid properties exhibit significant variations. The small characteristic thickness of the interface results in large gradients in physical quantities, necessitating a high-resolution mesh for accurate representation. Conversely, in regions farther from the interface, a coarser mesh suffices to maintain numerical accuracy. However, many of the studies referenced above utilize uniform meshes for simulating ferrofluid interfacial dynamics under magnetic field influence, which can lead to computational inefficiencies.

To tackle the inefficiencies associated with uniform meshing, adaptive mesh refinement (AMR) based on the lattice Boltzmann method (LBM) implemented on non-uniform meshes has demonstrated significant improvements in computational efficiency. AMR-LBM techniques intelligently allocate computational resources by refining the mesh only in regions close to interfaces, where high resolution is critical, while utilizing coarser meshes in areas farther away. Several implementation strategies for AMR-LBM have emerged. Tölke et al.[30] introduced an AMR-LBM for multiphase problems using a Rothman–Keller-type model,[31] which employs an unstructured tree-type mesh. Building on a multi-block structured mesh, Yu and Fan[32] developed an enhanced interaction potential model that utilizes AMR-LBM to simulate bubble rising problems. In this approach, the original "collision–propagation" algorithm in LBM is complicated by additional "explosion" and "coalescence" operations to transfer information between coarse and fine levels. Chen et al.[33] proposed a non-uniform mesh-based LBM combined with an image reconstruction technique to create a quadtree mesh for simulating fluid flow in porous media. Liu et al.[34] presented an adaptive mesh refinement method based on an octree structural representation, utilizing individual hash tables for each refinement level to avoid key value conflicts. Hasegawa et al.[35] developed a block-structured AMR-LBM leveraging an octree forest and employed GPU acceleration for aerodynamic simulations. Deiterding and Wood[36] applied block-structured AMR to lattice Boltzmann methods to address complex engineering challenges, highlighting the potential of AMR-LBM to enhance computational efficiency in simulating turbulent wake issues arising from wind turbine operations. It is noteworthy that in the aforementioned AMR-LB models, the standard collision-streaming procedures in LBM necessitate the use of different time steps on meshes with varying lattice sizes, leading to the need for both temporal and spatial interpolations. This results in a significant increase in algorithm complexity and imposes limits on the ratio of the coarsest to the finest mesh spacings. Furthermore, when employing a tree-type data structure in AMR, the memory and time required for tree traversal can be substantial. To mitigate these drawbacks, Fakhari et al.[37] introduced a finite-difference LBM with a block-structured AMR technique, utilizing pointer attributes to determine the neighbors of a specific block through appropriate adjustments of its child identifications. This approach employs the Lax–Wendroff scheme to maintain a uniform time step across meshes with different lattice sizes. Within this framework, both the center-of-cell approach and the apex-of-cell approach[38] were evaluated to assess the accuracy and efficiency of AMR in multiphase fluid simulations. The AMR-LBM has also been applied to three-phase flows and phase change problems.[39,40]

In this paper, we propose a block-structured AMR-LBM for two-phase ferrofluid flows to investigate interfacial dynamics under the influence of a magnetic field. The governing equations in discrete Boltzmann form for the phase field, flow field, and magnetic field are solved using the Lax–Wendroff scheme. Furthermore, to maintain a uniform magnetic field, we implement a specialized treatment of the Neumann boundary condition for the magnetic scalar potential within the Lax–Wendroff framework. Importantly, we note that most existing studies of two-phase ferrofluid flows focus primarily on droplet and bubble dynamics. However, flow instabilities in two-phase ferrofluid systems exhibit a richer variety of interfacial dynamic behaviors. With the aid of the proposed AMR-LBM, we conduct a detailed investigation into Rosensweig instability, Rayleigh–Taylor instability, and Kelvin–Helmholtz instability under different magnetic field configurations.

The remainder of the paper is organized as follows. In Sec. II, we introduce the governing equations for two-phase ferrofluid flows, the lattice Boltzmann scheme, and the implementation details of AMR on block-structured meshes. Section III presents simulations of several typical problems, including the deformation of ferrofluid droplets, bubble rising in ferrofluid, Rosensweig instability, Rayleigh–Taylor instability under a uniform vertical magnetic field, and Kelvin–Helmholtz instability under a uniform horizontal magnetic field. The results of these simulations will be discussed in detail. Finally, Sec. IV offers concluding remarks and outlines future perspectives.

---

## II. MATHEMATICAL MODEL AND NUMERICAL METHOD

### A. Governing equations for two-phase ferrofluid flows

The magnetic governing equations of ferrofluid flows are described by the Maxwell equations for nonconducting fluids[41]

$$
\nabla \cdot \mathbf{B} = 0, \tag{1}
$$

$$
\nabla \times \mathbf{H} = 0, \tag{2}
$$

where *H* (A/m) and *B* (T) are the magnetic field strength and the magnetic flux density, respectively. Considering a ferrofluid domain Ω_d, surrounded by a nonmagnetizable medium Ω_c, *B* can be expressed as

$$
\mathbf{B} = \begin{cases} \mu_0 (\mathbf{H} + \mathbf{M}), & \text{if } \Omega_d; \\ \mu_0 \mathbf{H}, & \text{if } \Omega_c. \end{cases} \tag{3}
$$

where **M** is the magnetization and it can be expressed as **M** = χ**H**. The vacuum permeability μ₀ is 4π × 10⁻⁷ N/A². μ can be expressed as μ = μ₀(1 + χ). μ and χ are the permeability and the magnetic susceptibility, respectively. According to the irrotational condition, a scalar magnetic potential is defined that satisfies

$$
\mathbf{H} = -\nabla \psi. \tag{4}
$$

Then, by substituting Eq. (3) and Eq. (4) into Eq. (1), we ultimately obtain the magnetic potential equation

$$
\nabla \cdot (\mu \nabla \psi) = 0. \tag{5}
$$

The binary ferrofluid flows involved in this paper are presumed to be incompressible, isothermal, and immiscible. The Navier–Stokes equations governing the motion of the two-phase fluids are

$$
\nabla \cdot \mathbf{u} = 0, \tag{6}
$$

$$
\rho \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\nabla p + \eta \nabla^2 \mathbf{u} + \mathbf{F}_s + \mathbf{F}_b + \mathbf{F}_m, \tag{7}
$$

where ρ is the density of the fluid, **u** the velocity, and *p* the pressure, η the kinematic viscosity, **F**_b and **F**_s the body force and surface tension, respectively; **F**_m denotes the magnetic force. When the magnetic susceptibility is constant corresponding to small magnetic field strengths, according to Hu et al.,[23] **F**_m can be expressed by

$$
\mathbf{F}_m = \mu_0 \chi \nabla \left ( \frac{|\mathbf{H}|^2}{2} \right). \tag{8}
$$

In this study, the second-order conservative Allen–Cahn equation is adopted to capture the phase interfaces and it reads[42]

$$
\frac{\partial \phi}{\partial t} + \nabla \cdot (\phi \mathbf{u}) = \nabla \cdot \left[ M_\phi \nabla \phi - \frac{1 - 4 (\phi - \phi_0)^2}{W} \mathbf{n} \right], \tag{9}
$$

where φ is the order parameter of phase field that indicates each phase. φ₀ = (φ_l + φ_h)/2 implies the interface location. In this study, we set φ_l = 0 and φ_h = 1. M_φ and W are the mobility and interface width parameter, respectively. **n** is the normal vector to the interface, which is calculated by

$$
\mathbf{n} = \frac{\nabla \phi}{|\nabla \phi|}. \tag{10}
$$

At equilibrium state, the phase field profile along the normal direction of the interface reads

$$
\phi (n) = \phi_0 + \frac{\phi_h - \phi_l}{2} \tanh\left (\frac{2 n}{W}\right), \tag{11}
$$

where *n* is the signed distance from the interface. The surface tension **F**_s can be calculated as

$$
\mathbf{F}_s = \lambda_\phi \nabla \phi, \tag{12}
$$

where λ_φ is the chemical potential, which can be computed by

$$
\lambda_\phi = 4\beta (\phi - \phi_l)(\phi - \phi_h)(\phi - \phi_0) - \kappa \nabla^2 \phi, \tag{13}
$$

where β and κ are the parameter relating to the surface tension coefficient σ and the interfacial thickness W:

$$
\beta = \frac{12\sigma}{W}, \quad \kappa = \frac{3 W\sigma}{2}. \tag{14}
$$

### B. Finite difference lattice Boltzmann scheme

The traditional lattice Boltzmann method (LBM) is not well suited for non-uniform meshes. To implement the "collision-streaming" process on the hierarchical mesh structure, the finite difference lattice Boltzmann method (FDLBM) proposed by Fakhari et al.[37] is adopted herein. Using a time-splitting method and Lax-Wendroff scheme, the discrete Boltzmann equation (DBE) with a multiple-relaxation-time collision operator for phase field equation can be discretized as

$$
\tilde{f}_\alpha (\mathbf{x}) = f_\alpha (\mathbf{x}) - M_9^{-1} S_f M_9^{\alpha\beta} [f_\beta (\mathbf{x}) - f_\beta^{eq}(\mathbf{x})] + M_9^{-1} (I - \frac{S_f}{2}) M_9^{\alpha\beta} F_\beta, \tag{15}
$$

$$
f_\alpha (\mathbf{x}, t + dt) = \left (1 - \frac{c^2}{2}\right) \tilde{f}_\alpha (\mathbf{x}, t) + \frac{c (c+1)}{2} \tilde{f}_\alpha (\mathbf{x} - d\mathbf{x}_\alpha, t) + \frac{c (c-1)}{2} \tilde{f}_\alpha (\mathbf{x} + d\mathbf{x}_\alpha, t), \tag{16}
$$

where f_α and f̃_α are the distribution function for the order parameter. Here, the D 2 Q 9 lattice velocity model is used, and the corresponding lattice discrete velocity **e**_α is given as

$$
\mathbf{e}_\alpha = \begin{cases} (0, 0), & \alpha = 0; \\ \left (\cos (\alpha-1)\frac{\pi}{2}, \sin (\alpha-1)\frac{\pi}{2}\right) c, & \alpha = 1, 2, 3, 4; \\ \left (\sqrt{2}\cos (2\alpha-1)\frac{\pi}{4}, \sqrt{2}\sin (2\alpha-1)\frac{\pi}{4}\right) c, & \alpha = 5, 6, 7, 8. \end{cases} \tag{17}
$$

where *c* is the lattice speed. It can be expressed as *c* = Δx/dt, where Δx and dt are the lattice spacing and the time step, respectively. The phase-field equilibrium distribution function f_α^eq is given as

$$
f_\alpha^{eq} = \omega_\alpha \phi \left (1 + \frac{\mathbf{e}_\alpha \cdot \mathbf{u}}{c_s^2}\right), \tag{18}
$$

where c_s = c/√3 is the lattice speed of sound and ω_α is the weight coefficient set with ω₀ = 4/9; ω₁₋₄ = 1/9; ω₅₋₈ = 1/36. c = |**e**_α|dt/|dx_α| is the Courant–Friedrichs–Lewy (CFL) number. dx_α is the displacement between two adjacent nodes along the direction of **e**_α. The discrete source term F_β can be calculated by

$$
F_\beta = \omega_\beta \mathbf{e}_\beta \cdot \frac{4\phi (1-\phi)}{W} \mathbf{n}. \tag{19}
$$

The 9×9 transform matrix M₉ reads

$$
M_9 = \begin{pmatrix} 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\ 4 & -1 & -1 & -1 & -1 & 2 & 2 & 2 & 2 \\ 4 & -2 & -2 & -2 & -2 & 1 & 1 & 1 & 1 \\ 0 & 1 & 0 & -1 & 0 & 1 & -1 & -1 & 1 \\ 0 & -2 & 0 & 2 & 0 & 1 & -1 & -1 & 1 \\ 0 & 0 & 1 & 0 & -1 & 1 & 1 & -1 & -1 \\ 0 & 0 & -2 & 0 & 2 & 1 & 1 & -1 & -1 \\ 0 & 1 & -1 & 1 & -1 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 1 & -1 & 1 & -1 \end{pmatrix}. \tag{20}
$$

The relaxation matrix S_f is

$$
S_f = \text{diag}[s_0^f, s_1^f, \ldots, s_8^f], \tag{21}
$$

where the diagonal elements are chosen as

$$
\frac{1}{s_3^f} = \frac{1}{s_5^f} = \frac{M_\phi}{c_s^2 dt} + 0.5, \tag{22}
$$

$$
s_0^f = s_1^f = s_2^f = s_4^f = s_6^f = s_7^f = s_8^f = 1. \tag{23}
$$

The phase-field value φ can be determined by taking the zeroth moment of the phase-field distribution function f_α

$$
\phi = \sum_{\alpha=0}^{8} f_\alpha. \tag{24}
$$

The Navier–Stokes equations can be reformulated in a velocity-based form, which is given by[29]

$$
\frac{\partial \mathbf{u}}{\partial t} + \mathbf{u} \cdot \nabla \mathbf{u} = -\nabla \frac{p}{\rho} + \nabla \cdot [\nu (\nabla \mathbf{u} + \mathbf{u}\nabla)] + \frac{\mathbf{F}^{total}}{\rho}, \tag{25}
$$

where ν = η/ρ is the kinematic viscosity, and **F**^total is defined as

$$
\mathbf{F}^{total} = \mathbf{F}_\nu + \mathbf{F}_p + \mathbf{F}_s + \mathbf{F}_b + \mathbf{F}_m, \tag{26}
$$

$$
\mathbf{F}_\nu = \nu (\nabla \mathbf{u} + \mathbf{u}\nabla) \cdot \nabla \rho, \tag{27}
$$

$$
\mathbf{F}_p = -\frac{p}{\rho} \nabla \rho. \tag{28}
$$

The finite difference LB equation for velocity-based Navier–Stokes equations can be written as

$$
\tilde{g}_\alpha (\mathbf{x}) = g_\alpha (\mathbf{x}) - M^{-1} S_g M_{\alpha\beta} [g_\beta (\mathbf{x}) - g_\beta^{eq}(\mathbf{x})] + M^{-1} (I - \frac{S_g}{2}) M_{\alpha\beta} G_\beta, \tag{29}
$$

$$
g_\alpha (\mathbf{x}, t+dt) = \left (1 - \frac{c^2}{2}\right) \tilde{g}_\alpha (\mathbf{x}, t) + \frac{c (c+1)}{2} \tilde{g}_\alpha (\mathbf{x} - d\mathbf{x}_\alpha, t) + \frac{c (c-1)}{2} \tilde{g}_\alpha (\mathbf{x} + d\mathbf{x}_\alpha, t), \tag{30}
$$

where g_α is the density distribution function. The equilibrium distribution function g_α^eq is

$$
g_\alpha^{eq} = \omega_\alpha \left[ \frac{p}{\rho c_s^2} + \frac{\mathbf{e}_\alpha \cdot \mathbf{u}}{c_s^2} + \frac{(\mathbf{e}_\alpha \cdot \mathbf{u})^2}{2 c_s^4} - \frac{\mathbf{u} \cdot \mathbf{u}}{2 c_s^2} \right]. \tag{31}
$$

The discrete force term can be obtained by

$$
G_\beta = \omega_\beta \left[ \frac{\mathbf{e}_\beta \cdot \mathbf{F}^{total}}{c_s^2} + \frac{\mathbf{u} \mathbf{F}^{total} : (\mathbf{e}_\beta \mathbf{e}_\beta - c_s^2 \mathbf{I})}{c_s^4} \right]. \tag{32}
$$

The diagonal elements s_i^g (0 ≤ i ≤ 8) in relaxation matrix S_g is

$$
\frac{1}{s_7^g} = \frac{1}{s_8^g} = \frac{\eta}{\rho c_s^2 \Delta t} + 0.5, \tag{33}
$$

$$
s_0^g = s_1^g = s_2^g = s_3^g = s_4^g = s_5^g = s_6^g = 1. \tag{34}
$$

With the hydrodynamic distribution function, the conserved quantities can be obtained as

$$
\mathbf{u} = \sum_\alpha g_\alpha \mathbf{e}_\alpha + \frac{dt}{2\rho} \mathbf{F}^{total}, \tag{35}
$$

$$
p = \rho c_s^2 \sum_\alpha g_\alpha. \tag{36}
$$

When solving the potential equation for magnetic field, as suggested by Hu et al.,[23] the original Laplace equation needs to be modified by introducing a pseudo-time derivative term as well as a free parameter ε. The resulting equation is expressed as

$$
\frac{1}{\varepsilon} \frac{\partial \psi}{\partial t} = \nabla \cdot (\mu \nabla \psi). \tag{37}
$$

To ensure the solution of Eq. (37) is approaching to that of Eq. (5), a large value of ε should be used. In our study, to ensure the numerical stability, ε = 1/μ₀ is used.

Unlike the phase field equation and the hydrodynamics equations, the magnetic field equation is relatively simple. As indicated by Li et al.,[43] the D 2 Q 5 model is sufficient to obtain the accurate solutions. In this study, the following finite difference LB model is also used to solve the magnetic field

$$
\tilde{h}_\alpha (\mathbf{x}) = h_\alpha (\mathbf{x}) - M_5^{-1} S_h M_5^{\alpha\beta} [h_\beta (\mathbf{x}) - h_\beta^{eq}(\mathbf{x})], \quad \alpha = 0, 1, 2, 3, 4; \tag{38}
$$

$$
h_\alpha (\mathbf{x}, t+dt) = \left (1 - \frac{c^2}{2}\right) \tilde{h}_\alpha (\mathbf{x}, t) + \frac{c (c+1)}{2} \tilde{h}_\alpha (\mathbf{x} - d\mathbf{x}_\alpha, t) + \frac{c (c-1)}{2} \tilde{h}_\alpha (\mathbf{x} + d\mathbf{x}_\alpha, t), \tag{39}
$$

where h_α is the distribution function for the magnetic potential. The local equilibrium distribution function h_α^eq can be expressed as

$$
h_\alpha^{eq} = \omega_\alpha^5 \psi, \tag{40}
$$

where ω_α^5 = 0.2 (0 ≤ α ≤ 4). Here, 5×5 transform matrix M₅ is expressed as

$$
M_5 = \begin{pmatrix} 1 & 1 & 1 & 1 & 1 \\ 0 & 1 & -1 & 0 & 0 \\ 0 & 0 & 0 & 1 & -1 \\ 4 & -1 & -1 & -1 & -1 \\ 0 & 1 & 1 & -1 & -1 \end{pmatrix}. \tag{41}
$$

The diagonal elements s_h^i in relaxation matrix S_h is

$$
\frac{1}{s_1^h} = \frac{1}{s_2^h} = \frac{2.5 \varepsilon \mu}{c^2 dt} + 0.5, \tag{42}
$$

$$
s_0^h = s_3^h = s_4^h = 1. \tag{43}
$$

The magnetic potential ψ is updated by taking the zeroth moment of the distribution function

$$
\psi = \sum_\alpha f_\alpha. \tag{44}
$$

In the actual simulations, the magnetic field H_n on the boundaries are usually given. In other words, the Neumann boundary condition ∇ψ · **n**_b = H_n should be handled, where **n**_b is the unit normal vector on the boundary. As shown in Fig. 1, a layer of virtual nodes is arranged outside and the actual boundary location is half a cell size away from the lattice nodes. Within the time interval [t_n, t_{n+1}], we can get the following relation based on the conservation law:

$$
\frac{\frac{c}{2}(c+1)\tilde{h}_2 (\mathbf{x}_0, t_n) + \frac{c}{2}(c-1)\tilde{h}_4 (\mathbf{x}_0, t_n)}{dx_2} = \frac{\frac{c}{2}(c-1)\tilde{h}_2 (\mathbf{x}_1, t_n) + \frac{c}{2}(c+1)\tilde{h}_4 (\mathbf{x}_1, t_n)}{dx_2} + dt (\varepsilon \mu H_n) dx. \tag{45}
$$

Then using the non-equilibrium extrapolation method, we have

$$
h_2 (\mathbf{x}_0, t_n) - \tilde{h}_4 (\mathbf{x}_0, t_n) = \tilde{h}_2 (\mathbf{x}_1, t_n) - \tilde{h}_4 (\mathbf{x}_1, t_n). \tag{46}
$$

Combining Eqs. (45) and (46) simultaneously yields

$$
h_2 (\mathbf{x}_0, t_n) = \left (1 - \frac{1}{c}\right) h_2 (\mathbf{x}_1, t_n) + \frac{1}{c} h_4 (\mathbf{x}_1, t_n) + \frac{\varepsilon \mu H_n}{c}, \tag{47}
$$

$$
h_4 (\mathbf{x}_0, t_n) = -\frac{1}{c} h_2 (\mathbf{x}_1, t_n) + \left (1 + \frac{1}{c}\right) h_4 (\mathbf{x}_1, t_n) + \frac{\varepsilon \mu H_n}{c}. \tag{48}
$$

Once the post-collision distribution functions on the virtual nodes j = 0 are obtained, the streaming step can be performed.

At last, we notice that the physical property parameters including density, viscosity, and permeability can be updated by a linear interpolation after obtaining the phase-field value

$$
\rho = \rho_l + \phi (\rho_h - \rho_l), \tag{49}
$$

$$
\eta = \eta_l + \phi (\eta_h - \eta_l), \tag{50}
$$

$$
\mu = \mu_l + \phi (\mu_h - \mu_l), \tag{51}
$$

where the subscripts l and h represent the light and heavy fluids, respectively. Moreover, we compute the gradient term ∇φ and Laplacian term ∇²φ using second-order isotropically centered differences[44]

$$
\nabla \phi = \frac{c}{c_s^2 dt} \sum_{\alpha=1}^{8} \frac{1}{2} e_{\alpha x} \omega_\alpha [\phi (\mathbf{x} + \mathbf{e}_\alpha dt) - \phi (\mathbf{x} - \mathbf{e}_\alpha dt)], \tag{52}
$$

$$
\nabla^2 \phi = \frac{2 c^2}{c_s^2 dt} \sum_{\alpha=1}^{8} \omega_\alpha [\phi (\mathbf{x} + \mathbf{e}_\alpha dt, t) - \phi (\mathbf{x} - \mathbf{e}_\alpha dt, t)]. \tag{53}
$$

### C. Block-structured adaptive mesh refinement

In this work, we implement the adaptive mesh refinement (AMR) method using pointers instead of a tree data structure. This approach utilizes pointers to efficiently track all neighboring and child IDs without the need for maintaining or modifying a tree structure. Given that all blocks are structured and exhibit a degree of self-similarity, we employ a data structure to store relevant information. This structure encompasses both intra-block data—such as coordinates, distribution functions, and macroscopic quantities—and inter-block data, which includes neighbor information.

The AMR hierarchy is constructed using structured rectangular blocks, with each block consisting of n_x × n_y cells. The mesh spacing of these cells varies and is characterized by the refinement level l. The initial blocks with l₀ are referred to as root blocks, while the finest blocks with the maximum refinement level l_max are referred to as leaf blocks. The relationship between the mesh spacing Δx, the refinement level l, and the size of the root block (L_x × L_y) is expressed as

$$
\Delta x = \frac{L_x}{n_x \cdot 2^l} = \frac{L_y}{n_y \cdot 2^l}. \tag{54}
$$

As shown in Fig. 2, refining a mesh block once results in the creation of four child mesh blocks, with the mesh resolution doubling. Block information is updated accordingly, including the parent block identifier (ID), refinement level, four pointers to child blocks, eight pointers to orthogonal and diagonal neighboring blocks, and a logical variable indicating whether the block is a leaf block. Refinement continues until l_max or the refinement criteria are satisfied.

The refinement criteria are based on the gradient of phase field parameter,[37] expressed as

$$
\epsilon = |\nabla \phi| = \sqrt{\left (\frac{\partial \phi}{\partial x}\right)^2 + \left (\frac{\partial \phi}{\partial y}\right)^2}, \tag{55}
$$

where ε represents the refinement error estimate. A block is marked for refinement when ε ≥ ε_r = 0.002, while it is marked for coarsening when ε ≤ ε_d = 0.001. If ε_d < ε < ε_r, the block retains its current refinement level, undergoing neither refinement nor coarsening.

To facilitate information transfer between blocks, each block must have access to extract data from its neighboring blocks. This is achieved by using a ghost layer surrounding each block. When filling the ghost layer of a block, there are three possible scenarios, as illustrated in Fig. 3. Here, l_A and l_B denote the refinement levels of blocks A and B, respectively. The bold solid lines represent block interfaces, the thin solid lines represent mesh lines, the solid circles indicate mesh points, and the dashed circles indicate ghost points.

As shown in Fig. 3 (a), the simplest case occurs when two adjacent blocks are at the same refinement level, and both blocks A and B are leaf blocks. The ghost points of block A are assigned by directly copying data points from block B.

The second scenario, as illustrated in Fig. 3 (b), occurs when two adjacent blocks are on the same hierarchical level, but block B contains sub-blocks. To fill the ghost points of the coarser block A using data from the finer block B, the biquadratic interpolation formula for interior nodes, shown in Fig. 4 (a), is applied

$$
f (I, J) = [9 f (i, j) + 18 f (i+1, j) - 3 f (i+2, j) + 18 f (i, j+1) + 36 (i+1, j+1) - 6 f (i+2, j+1) - 3 f (i, j+2) - 6 f (i+1, j+2) - 6 f (i+2, j+1)]/64 + o (\Delta x^3). \tag{56}
$$

The final scenario, shown in Fig. 3 (c), occurs when blocks A and B have different refinement levels, such that l_A = l_B + 1. In this case, the ghost points of the finer block A are filled using data from the coarser block B. As depicted in Fig. 4 (b), the biquadratic interpolation formula for interior nodes is employed

$$
f (i, j) = [25 f (I-1, J-1) + 150 f (I, J-1) - 15 f (I+1, J-1) + 150 f (I-1, J) + 900 f (I, J) - 90 f (I+1, J-1) - 15 f (I-1, I+1) - 90 f (I, I+1) + 9 f (I+1, I+1)]/1024 + o (\Delta x^3). \tag{57}
$$

It should be noted that when filling the ghost points of the fine block, the ghost points of the adjacent coarse block are needed, as indicated by the rectangular block in Fig. 4 (b). Therefore, the ghost points of the coarse block should be filled first, and then the ghost points of the fine block should be filled.

---

## III. RESULTS AND DISCUSSION

In this section, several typical benchmark problems are simulated to evaluate the performance of the adaptive mesh refinement (AMR) model. First, the magnetic field strength within a stationary cylinder is simulated and compared with its analytical solution to validate the accuracy of the magnetic field model. Next, numerical experiments are conducted to investigate the deformation of ferrodroplets, bubble buoyancy, and Rosensweig instability in magnetic fluids under uniform magnetic fields, thereby demonstrating the model's capability to simulate multiphysics coupling in two-phase flows while also assessing its computational efficiency. Finally, to enrich the interfacial dynamics of two-phase ferrofluid flows, the Rayleigh–Taylor instability and Kelvin–Helmholtz instability are examined at both low and high Reynolds numbers in the presence of magnetic fields.

### A. A circular cylinder in a uniform magnetic field

An analytical solution is provided to determine the magnetic field strength inside a fixed cylinder subjected to an external uniform magnetic field, serving as a benchmark for verifying the model's magnetic field computation. In the polar coordinate system, the Laplace equation for the magnetic potential is expressed as follows:[23]

$$
\frac{\partial}{\partial r}\left (r \frac{\partial \psi}{\partial r}\right) + \frac{1}{r}\frac{\partial^2 \psi}{\partial \theta^2} = 0, \tag{58}
$$

where r and θ represent the radial and angular coordinates, respectively. Using the method of separation of variables and the properties of Legendre polynomials, the analytical solution of Eq. (58) is derived as

$$
\psi = \begin{cases} Ar \sin\theta, & r \leq R; \\ \left (Cr + \frac{D}{r}\right) \sin\theta, & r > R. \end{cases} \tag{59}
$$

where R denotes the radius of the cylinder. The corresponding magnetic field is expressed as

$$
\mathbf{H} = -\nabla\psi = \begin{cases} A\sin\theta\,\mathbf{e}_r - A\cos\theta\,\mathbf{e}_\theta, & r \leq R; \\ \left (\frac{D}{r^2} - C\right)\sin\theta\,\mathbf{e}_r - \left (\frac{D}{r^2} + C\right)\cos\theta\,\mathbf{e}_\theta, & r > R. \end{cases} \tag{60}
$$

where **e**_r and **e**_θ are the unit vectors in the radial and circumferential directions, respectively. H₀ represents an external uniform magnetic field, while A, C, and D are constants

$$
A = -\frac{2\mu_2}{\mu_1 + \mu_2} H_0, \quad C = -H_0, \quad D = \frac{\mu_1 - \mu_2}{\mu_1 + \mu_2} R^2 H_0. \tag{61}
$$

In this simulation, a cylinder with a radius of R is placed at the center of the computational domain. The computational domain is discretized into L × L lattice cells, where L = 256 and R = L/10. The boundary conditions at the bottom and top are

$$
\frac{\partial \psi}{\partial y} = H_0. \tag{62}
$$

The left and right boundaries are subjected to magnetic insulation conditions

$$
\frac{\partial \psi}{\partial x} = 0. \tag{63}
$$

The magnetic field lines, magnetic field strength, and mesh block distribution when μ₁/μ₂ = 2 are shown in Fig. 5. Evidently, the magnetic field lines inside the cylinder and near the outer boundary are aligned with the applied external magnetic field. However, due to the jump in permeability at the interface, the magnetic field lines are distorted near the cylinder. Additionally, it can be observed that the magnetic field inside the cylinder remains uniform, as expected in the analytical solution. The mesh block distribution is denser near the interface and sparser further away, which aligns with our expectations. In Fig. 6, a comparison is made between the numerical results of the magnetic field strength inside the cylinder and the analytical solution as the magnetic permeability ratio μ₁/μ₂ varies. Compared with the analytical solution, consistent results observed.

### B. Ferrodroplet deformation

In engineering applications, various problems are related to the deformation of ferrofluid droplets under the influence of magnetic fields. Examples include microfluidic cancer cell sorting based on magnetron control,[45] disease diagnosis,[46] and chemical engineering.[47] Flament et al.[48] experimentally studied the deformation of ferrofluid droplets confined in a narrow gap between two parallel plane layers, while also measuring the surface tension of the ferrofluid. Hu et al.[23] simulated the deformation of ferrofluid droplets using the same physical parameters as Flament's experiment, employing an LBM phase-field model. He et al.[29] also conducted numerical experiments on ferrofluid droplets, utilizing an immersed interface and phase-field lattice Boltzmann model with the same physical parameters as those in Flament's experiment. Their experimental and numerical results thus serve as benchmarks for validating the ferrofluid multiphase flow model in this paper.

As shown in Fig. 7, an organic solvent is filled in the square cavity and a magnetic liquid droplet immiscible with the organic solvent is placed in the center of the square cavity. A uniform magnetic field is applied in the vertical direction. To validate the model, parameters identical to those in Flament's experiment are used. The length of the square domain is set to L = 7.5 mm, and the radius of the ferrofluid droplet is r = 0.78 mm. The densities of the organic solvent and ferrofluid are 0.8 × 10³ and 1.58 × 10³ kg/m³, respectively. The viscosities are 1.0 × 10⁻³ Pa·s for the organic solvent and 4.0 × 10⁻³ Pa·s for the ferrofluid. The surface tension coefficient is set to σ = 3.07 mN/m, and the relative magnetic permeability is fixed at 2.2. In the simulation, the mobility is set to 0.01, and the interface thickness is set to 4. The lattice unit physical parameters are ρ_l = 1, ρ_h = 1.975, η_l = 0.0142, η_h = 0.0568, and σ = 0.0146. Non-slip boundary conditions are applied at the top and bottom of the computational domain, while periodic boundary conditions are applied at the left and right sides.

Figure 8 shows the deformation of the droplet interface under magnetic field strengths of 1.2, 2.4, 2.9, and 3.7 kA/m in the vertical direction, as well as the comparison of aspect ratios with previous studies. The combined effects of surface tension and magnetic force cause the deformation of the ferrofluid droplet. As shown in Fig. 8, when the magnetic field strength is small, the ferrofluid droplet remains nearly spherical. However, as the magnetic field strength increases, the droplet deforms into an elliptical shape. Moreover, the simulated results align well with published experimental results,[23,29] thus confirming the accuracy of the AMR-LB model's computations.

Figure 9 shows the variation of computational time and mesh number for the ferrofluid droplet at H = 3.7 kA/m under adaptive mesh refinement. It can be concluded that, using adaptive mesh refinement requires less computational time and fewer mesh points to achieve the same accuracy compared to uniform meshes. The computational time is reduced by a factor of 3.5, and the mesh number is reduced to 0.27 times that of the uniform mesh. Figure 10 displays the mesh block structure for both uniform and adaptive meshes at a time step of t = 12000 in the computational domain. It can be observed that, in the regions where numerical changes are significant, the adaptive mesh has a higher refinement level, similar to the uniform mesh. However, in areas away from these regions, the mesh is relatively sparse, which greatly enhances computational efficiency.

### C. Bubble rising in ferrofluid

This subsection uses the AMR-LB model to simulate the bubble rising problem under a uniform magnetic field.[49] We consider the motion of a bubble surrounding ferrofluid in a rectangular channel. Initially, a circular bubble with diameter D = 100 is placed in a computational domain with a length of L = 256 and a height of 2 L, where the mesh is divided into 256 × 512. The coordinates of the bubble center are (L/2, L/2). A vertical uniform magnetic field H₀ is applied. The initial position of the bubble is illustrated in Fig. 11. Periodic boundary conditions are applied in the horizontal direction, and no-slip boundary conditions are applied in the vertical direction. The body force is considered **F**_b = ρG_y, with G_y being the magnitude of gravitational acceleration in the vertical direction. Through dimensionless analysis, in addition to the density ratio (ρ_h/ρ_l) and viscosity ratio (η_h/η_l), the Reynolds number (Re), the Eötvös number (Eo), and the magnetic Bond number (Bo_m) are defined as

$$
Re = \frac{\sqrt{|G_y|\rho_h^2 D^3}}{\eta_h}, \tag{64}
$$

$$
Eo = \frac{|G_y|\rho_h D^2}{\sigma}, \tag{65}
$$

$$
Bo_m = \frac{\mu_0 H_0^2 D}{2\sigma}. \tag{66}
$$

The dimensionless time reads

$$
t^* = t\sqrt{\frac{|G_y|}{D}}. \tag{67}
$$

First, the evolution process of the bubble under buoyancy without the influence of a magnetic field is simulated. Three sets of simulations are performed for different values of the Eo (Eo = 5, 20, 125, and 400). The density ratio is set to 1000, the viscosity ratio to 100, and the Reynolds number to 40. Figure 12 shows the bubble shapes at t* = 1, 2, 3, and 4 for Re = 40 and Eo = 125. Figure 13 shows the bubble shapes at t* = 4 for Re = 40 and Eo = 5, 20, 125, and 400. From the simulation results, we observe that at lower Eo values, the bubble evolves into a two-dimensional elliptical shape. At medium Eo values, the bubble evolves into a two-dimensional hemispherical shape, while at higher Eo values, two filament structures behind bubble starts to form. Moreover, the experimental results are consistent with those observed by Fakhari et al.,[37] Hu et al.,[49] and Li et al.,[50] confirming the accuracy of the AMR-LB model.

In order to further verify the accuracy of the AMR-LB model, we quantitatively compared the center of mass Y_c changes during the bubble rising, which can be expressed as

$$
Y_c = \frac{\sum_{i, j} y \, H_d 5 (0.5 - \phi)}{\sum_{i, j} H_d 5 (0.5 - \phi)}, \tag{68}
$$

where H_d is the smeared-out Heaviside function

$$
H_d (x) = \begin{cases} 0, & x < -2; \\ 0.5 + 0.25x + \frac{0.5}{\pi}\sin(0.5\pi x), & -2 \leq x \leq 2; \\ 1, & x > 2. \end{cases} \tag{69}
$$

We recorded the vertical position of the center of mass Y_c during bubble rise. The numerical experiment uses the same parameters as the previous studies. The comparison results between them are shown in Fig. 14, it can be seen that the obtained results are in good agreement with those of Hu et al.[49] and Aland and Voigt,[51] further verifying the accuracy of the AMR-LB model.

Figure 15 shows the variation in computation time and mesh numbers for the bubble rising simulation under an adaptive mesh, with Re = 40, and Eo = 20, density ratio of 1000, and viscosity ratio of 100. The results indicate that, the use of adaptive mesh leads to a shorter computation time and fewer mesh points for achieving the same accuracy in the bubble rising experiment. At the time step of 20 000, the required computation time is reduced by a factor of 5.10, and the number of meshes is reduced to 0.19 times that of the uniform mesh. This significant improvement in computational efficiency further validates the effectiveness of the AMR-LB model.

Next, we simulate the evolution of bubbles in ferrofluid under gravity field with the presence of a magnetic field. The parameters are set as follows: Eo = 20, density ratio 1000, viscosity ratio 100, and Reynolds number Re = 40. The simulation is conducted with a magnetic field applied, which is uniform and vertical (parallel to gravity). This setup allows us to study the effect of the magnetic field on the evolution of bubbles in a magnetic fluid, and observe the influence of different magnetic field strengths on bubble deformation and behavior under buoyancy.

The results show that due to the normal force acting on the interface, the bubble elongates along the direction of the magnetic field. This elongation is influenced by the interaction between the surface tension, buoyant force, upward motion, pressure variation along the bubble interface, and viscous stress. Figure 16 shows the bubble shapes at t* = 1 for different magnetic Bond numbers (Bo_m = 0, 0.13, 0.27, 0.49, 0.76, and 1.94). It can be seen that, for bubbles in the ferrofluid, when the magnetic field strength is relatively low, the magnetic force is weak, and the bubbles remain relatively elongated in the horizontal direction. As the magnetic field strength increases, the magnetic force becomes stronger, causing the originally horizontally elongated bubble to stretch in the vertical direction and a filament structure is formed behind the bubble. Furthermore, as the magnetic field strength increases, the degree of vertical elongation increases too. Figure 17 shows the mesh block distribution under the corresponding bubble shapes for different magnetic Bond numbers (Bo_m = 0, 0.13, 0.27, 0.49, 0.76, and 1.94) at t* = 1, where each block consists of 8 × 8 mesh cells. It can be seen that the mesh is denser in the region near the bubble interface and sparser in the region away from the interface.

To quantify the simulation results, we measured the longitudinal length of the bubbles at different Bo_m numbers, as well as the bubble velocity, and the results are shown in Fig. 18. Here, the bubble rise velocity is defined as

$$
U_y = \frac{\int_{\phi < 0.5} u_y dx}{\int_{\phi < 0.5} dx}, \tag{70}
$$

where u_y is the vertical component of the velocity. The results show that in the absence of a magnetic field, the elongation length is always smaller than the bubble diameter and gradually decreases until it stabilizes, indicating that the bubble undergoes horizontal stretching only. Under the influence of a magnetic field, the elongation length exceeds the value observed without the magnetic field. It first increases, then decreases, and eventually reaches an equilibrium state. If the magnetic field is too strong, the elongation increases, then decreases until the vertical length becomes zero, causing the bubble to split into two sub-bubbles. A larger Bond number is also associated with a higher bubble rise velocity.

### D. Rosensweig instability

This subsection uses the AMR-LB model to simulate the Rosensweig instability of ferrofluids under the influence of a uniform magnetic field. Rosensweig instability refers to the phenomenon exhibited by ferrofluids when exposed to a magnetic field, where various patterns, such as peaks, valleys, or ripples, form on the surface. This instability arises due to the competition between magnetic forces, gravity, and surface tension within the fluid layer. When the magnetic force reaches a balance with gravity and surface tension, the growth of these patterns stops.[51,52]

As shown in Fig. 19, the immiscible binary fluids, an organic solvent and ferrofluid, fill the cavity, with the former occupying 2/3 of the volume and the latter occupying 1/3. The same parameters as in Flament's experiment[53] are used, and a uniform magnetic field is applied in the vertical direction. The length of the computational domain is set to L_x = 21 and L_y = 7 mm. The densities of the organic solvent and the ferrofluid are 0.8 × 10³ and 1.58 × 10³ kg/m³, respectively. The viscosities are 0.8 × 10⁻³ and 3.2 × 10⁻³ Pa·s, respectively. The surface tension coefficient is set to σ = 3.89 mN/m, and the relative magnetic permeability is fixed at 2.2. The mobility is set to 0.01, and the interface thickness is set to 4. The lattice unit physical parameters are ρ_l = 1, ρ_h = 1.975, η_l = 0.0073, η_h = 0.292, and σ = 0.0071, respectively. No-slip boundary conditions are applied at the top and bottom boundaries, and periodic boundary conditions are applied at the left and right boundaries.

Figure 20 shows a comparison of the morphology and structure of the peak formation of ferrofluid under the application of a uniform magnetic field in the vertical direction, with field strengths of 0, 4.9, and 8.2 kA/m. The figure also presents a comparison between experimental results (on the left) and numerical simulations (on the right). It can be seen clearly that when there is no magnetic field, the surface of the ferrofluid remains stationary. As the magnetic field strength exceeds 4.9 kA/m, wave-like disturbances form at the interface. When the field strength increases further to 8.2 kA/m, the disturbances become more pronounced, forming a comb-like structure. No hysteresis was observed in the simulations, and the magnetic susceptibility χ₀ < 2.54, confirming that the peak instability is a second-order transition.[54] The critical magnetic field is H_c = 4.7 kA/m, the critical wavelength is λ_c = 2.82 mm, and the critical magnetization strength is M_c ≈ 10.35 kA/m. The simulation results are in good agreement with the experimental results, further verifying the accuracy of the model.

Figure 21 shows the variation in computational time and mesh number for the Rosensweig instability of ferrofluids under an applied magnetic field of 4.9 kA/m using the adaptive mesh method. The computation time is reduced by a factor of 3.13, and the number of mesh cells is reduced to 0.31 of that required by a uniform mesh.

It should be pointed out that Krakov et al.[55] highlighted that mesh resolution, magnetic field strength, and the thickness of the diffused peak significantly influence the Rosensweig instability. To explore the effect of mesh resolution on wavelength of Rosensweig instability, we varied the number of mesh cells by adjusting the maximum refinement level and the base block size. Cowley and Rosensweig, in their pioneering work,[56] proposed an expression for the relationship between key parameters affecting Rosensweig instability. The critical magnetic field H_c and critical wavelength λ_c are given as follows:

$$
H_c = 2\sqrt{\frac{l_0}{l_0/\mu + 1}} \left (\frac{l_0}{\mu} - 1\right)^2 \left (\frac{1}{2}\right)^{1/2} (\sigma g \Delta\rho)^{1/4}, \quad \lambda_c = 2\pi\sqrt{\frac{\sigma}{g\Delta\rho}}^{1/2}, \tag{71}
$$

where Δρ = ρ_h - ρ_l. Equation (71) was used as a grounded estimate to predict the number of peaks in the container. We set the magnetic susceptibility to χ₀ = 0.9, the ferrofluid density to ρ_h = 1.58 × 10³ kg/m³, the surface tension to σ = 5.0 mN/m and non-magnetic fluid density to ρ_l = 0.8 × 10³ kg/m³. Assuming λ_c = 2.5 mm, the gravitational acceleration g = 40.49 m/s² was derived using the above equations.

In this simulation, the geometry configurations and boundary conditions of physical fields are same as previous test. Simulations are conducted under different maximum refinement levels (with L = 256, maximum refinement levels of 3, 4, and 5, and L = 512, with a maximum refinement level of 5) to investigate the effects on the peak formation. For each mesh resolution, the interface thickness ε spanned four mesh cells. The simulation results considered six complete wavelengths within the computational domain, where the theoretical value for six full wavelengths should be 15 mm.

Figure 22 presents a comparison of steady-state peak shapes under different mesh sizes, with a magnetic field strength of H₀ = 10 kA/m. The results show that the mesh size significantly affects the uniformity and wavelength of the peak. As the mesh was refined, the six wavelength lengths gradually converged from 16.1 mm to the theoretical value of 15 mm. This matches well with the simulation results for L = 256 and a refinement level of 5. Therefore, to further optimize simulation time and reduce computational cost, simulations were performed using a mesh of L = 256 with a refinement level of 5.

Finally, the effect of magnetic field strength on the Rosensweig instability was investigated by conducting simulations at various magnetic field strengths: H₀ = 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, and 18 kA/m. The numerical results indicate that the magnetic field strength has a significant impact on the peak height and trough depth. Specifically, as the magnetic field strength increases, the peak height and trough depth both increase, with the peak height and trough depth gradually approaching each other. Figure 23 presents a comparison of peak heights and trough depths at different magnetic field strengths. However, beyond a certain magnetic field strength, the peaks will no longer continue to grow.

### E. Rayleigh–Taylor instability in a uniform magnetic field

This subsection uses the AMR-LB model to simulate the Rayleigh–Taylor instability of ferrofluid under the influence of a uniform magnetic field. The Rayleigh–Taylor instability typically occurs at the interface between two fluids of different densities. Under the influence of gravity, the denser fluid tends to move downward due to its higher density, while the lighter fluid moves upward. A small initial perturbation at the interface is amplified over time under the continuous influence of gravity. This instability caused by interface deformation is referred to as Rayleigh–Taylor instability. The phenomenon is characterized by interface deformation, progressing from ripples or wrinkles to mushroom-like or spike-shaped structures.

The computational domain for this problem is set as a rectangular region of [0, L] × [0, 4 L]. Within this domain, the top and bottom boundaries are assigned no-slip solid boundary conditions, while periodic boundary conditions are applied in the horizontal direction. The top region contains the heavier fluid (ρ_h, η_h), and the lighter fluid (ρ_l, η_l) occupies the bottom region. The initial interface is a flat line located at y = 2 L, and it is perturbed using a cosine function

$$
y = 2 L + 0.1 L \cdot \cos\left (\frac{2\pi x}{L}\right), \tag{72}
$$

where L = 256.

In order to compare the results of the AMR-LB model with other models in the literature, we define the dimensionless the Atwood number At, the Reynolds number Re, the Capillary number Ca, the magnetic Bond number Bo_m, and dimensionless time t* as follows:

$$
At = \frac{\rho_h - \rho_l}{\rho_h + \rho_l}, \tag{73}
$$

$$
Re = \frac{\rho_l U_0 L}{\eta_l}, \tag{74}
$$

$$
Ca = \frac{\eta_l U_0}{\sigma}, \tag{75}
$$

$$
Bo_m = \frac{\mu_0 H_0^2 L}{2\sigma}, \tag{76}
$$

$$
t^* = t\sqrt{\frac{L}{(g \cdot At)}}, \tag{77}
$$

where U₀ = √(gL) is the reference characteristic velocity.

In this simulation, we define the position change of the light fluid at x = 0 as the advancement of the light fluid front and the position change of the heavy fluid at x = L/2 as the displacement of the heavy fluid tip. Let √(gL) = 0.01, ρ_l = 1.0, At = 0.5, Ca = 1, and Re = 256. The interface thickness W spans 5 mesh cells. Figure 24 shows the interface morphology of the Rayleigh–Taylor instability problem at Re = 256 with dimensionless times at t* = 0, 1, 2, 3, 4, and 5, respectively. Then we let Re = 3000, other parameters remain unchanged. Figure 25 shows the interface morphology of the Rayleigh–Taylor instability problem at Re = 3000 with dimensionless times at t* = 0, 0.5, 1, 1.5, and 2, respectively. Figure 26 (a) quantitatively present the evolution of the light fluid front and the heavy fluid tip positions over time, while Fig. 26 (b) quantitatively display the velocity changes of the light fluid front and the heavy fluid tip over time. These results are compared with previously reported findings by Nourgaliev et al.[57] and Wang et al.[58,59] The comparison reveals that our AMR-LB method aligns closely with the results from the referenced studies, confirming the applicability of our numerical method under high Reynolds number conditions.

Subsequently, we applied a uniform vertical magnetic field to investigate the effect of the magnetic field on interface deformation. Magnetic fields with five different intensities, Bo_m = 0, 2.57, 16.08, 41.18, and 64.34 are applied. Figure 27 illustrates the interface morphology of the Rayleigh–Taylor instability problem at Bo_m = 41.18, and dimensionless times t* = 0, 0.5, 1, 1.5, and 2. It can be observed that, under the influence of the magnetic field, the tail behind the heavy fluid as it penetrates the light fluid is significantly reduced. The magnetic field suppresses the reverse rotational vortices. Compared to the scenario without a magnetic field, the motion speeds of the light fluid front and the heavy fluid tip are faster, enabling them to reach the same positions more quickly than that in the absence of a magnetic field.

Finally, to quantitatively represent the effect of different magnetic field intensities on Rayleigh–Taylor instability, we examined the position changes of the light fluid front and the heavy fluid tip over time under varying magnetic field strengths. Figure 28 show the position and velocity changes of the heavy fluid tip over time for Bo_m = 0, 2.57, 16.08, 41.18, and 64.34. Similarly, Fig. 29 illustrates the position and velocity changes of the light fluid front over time under the same magnetic field intensities. The numerical results indicate that the magnetic field intensity has a significant impact on the positions, velocities, and interface morphology of both the heavy fluid tip and the light fluid front. Specifically, as the magnetic field intensity increases, more pronounced positional changes and higher velocities are observed.

### F. Kelvin–Helmholtz instability in a uniform magnetic field

In this subsection, the AMR-LB model is employed to simulate the Kelvin–Helmholtz instability of a ferrofluid influenced by a horizontally oriented uniform magnetic field.[60] The Kelvin–Helmholtz instability typically arises at the interface between two fluids with different densities and relative velocities. When there is a velocity difference between the two fluid layers, a small perturbation occurs due to a combination of factors, including fluid viscosity, gravity, and interfacial tension. Under certain conditions, this perturbation is amplified, ultimately resulting in tumbling and mixing at the fluid interface, which creates an unstable state.[37,61]

As shown schematically in Fig. 30, the computational domain is set as a square region of [0, 2 L] × [0, 2 L], where L = 256, within which the top and bottom boundaries are set as no-slip boundary conditions, with periodic conditions on the left and right boundaries. The top half of the domain is filled by an organic solvent (ρ_h, η_h), while the ferrofluid (ρ_l, η_l) is located at the bottom half of the domain. The initial interface position is a flat line at y = 2 L, which is subsequently perturbed by a cosine function

$$
y (x) = 2 L + 0.1 L \cos\left (\frac{2\pi x}{L}\right). \tag{78}
$$

The initial order parameter and velocity field can be expressed as

$$
\phi (i, j) = \phi_0 + \frac{\phi_h - \phi_l}{2} \tanh\left (\frac{2 (y - y (x))}{W}\right), \tag{79}
$$

$$
\mathbf{u}(i, j) = \frac{1}{2}(2\phi (i, j) - 1) \cdot u_0. \tag{80}
$$

We then define the dimensionless parameters, the Reynolds number (Re), the Weber number (We), the Froude number (Fr), the magnetic Bond number (Bo_m), and dimensionless time (t*) as follows:

$$
Re = \frac{\rho_h U_0 L}{\eta_h}, \tag{81}
$$

$$
We = \frac{\rho_h L U^2}{\sigma}, \tag{82}
$$

$$
Fr = \frac{U^2}{gL}, \tag{83}
$$

$$
Bo_m = \frac{\mu_0 H_0^2 L}{2\sigma}, \tag{84}
$$

$$
t^* = t\sqrt{\frac{g}{L}}. \tag{85}
$$

Let us consider the scenario without a magnetic field; first, the physical parameters are set as: L = 256, U = 0.02, ρ_l/ρ_h = 0.99, We = 10 000, and Fr = 1.0. Three different cases with Re = 500, 5000, and 50 000 are simulated. Figure 31 shows the interfacial evolution of the Kelvin–Helmholtz instability for Re = 500, 5000, and 50000. At lower Reynolds numbers (Re = 500), the peak of the heavier fluid is elongated by the flow field of the lighter fluid, leading to the formation of comma-like liquid bridges under the combined effects of gravity and surface tension. A similar phenomenon is observed at higher Re values. As Re increases, the interfacial deformation intensifies, and the peak height increases slightly due to the weakening of the viscous effect. For Re = 5000, the comma-shaped liquid bridge curls further into a helical structure after t* = 4, with the bridge becoming less stable. When Re = 50 000, significant interfacial deformation occurs, forming an elongated liquid bridge, which eventually cracks after t* = 6. The deformation intensifies as the peak height increases slightly. As Re increases, the interfacial deformation becomes progressively more unstable. Before the liquid bridge ruptures, the number of helical turns increases with higher Re. Figure 31 illustrates the time evolution of the Kelvin–Helmholtz unstable interface at Re = 500, 5000, and 50000. At the early stages, shear flow at the interface generates vorticity in the x-direction. The peaks of the heavier and lighter fluids are elongated in the y-direction on the positive and negative sides, respectively, causing the two fluids to rotate and interpenetrate. The interfacial peaks are stretched and wrinkled by the x-direction vorticity, ultimately forming a helical fluid bridge. Subsequently, catenary-like structures appear in the vorticity field, in better agreement with the phenomenon of Fakhari et al.[37] and Li et al.[50]

Then a horizontally oriented uniform magnetic field is applied to investigate the effect of the magnetic field on the Kelvin–Helmholtz instability. Seven different magnitudes of magnetic field strengths of Bo_m = 0, 636, 1434, 2547, and 3980 are used, respectively. Figure 32 shows the interfacial evolution of the Kelvin–Helmholtz instability at Re = 5000, We = 10 000, Fr = 1, and different Bo_m = 636, 1434, and 2547. Figure 33 shows the two-phase interface morphology at t* = 5 under Bo_m = 0, 636, 1434, 2547, and 3980, respectively.

It can be observed that the applied magnetic field exerts an inhibitory effect on Kelvin–Helmholtz instability, specifically manifested in the varying degree of influence on the interface deformation of the two-phase liquid. When the strength of the applied magnetic field is low, the curling of the two-phase liquid experiences some inhibition, although the effect is minimal. As the magnetic field strength increases, the inhibitory effect becomes more pronounced, leading to a significant reduction in the degree of curling at the interface. When the magnetic field strength reaches a certain threshold, the interface of the two-phase liquid tends to flatten. Specifically, the horizontal component of the magnetic field inhibits Kelvin–Helmholtz instability, with greater magnetic field strengths resulting in a stronger inhibitory effect.

To quantify the effect of different magnetic field strengths on the Kelvin–Helmholtz instability, we investigate the height of the two-phase liquid mixing layer over time at different magnetic field strengths. Figure 34 shows the plots of the height of the two-phase liquid mixed layer as a function of time for Re = 5000, We = 10000, Fr = 1, and different Bo_m values. The numerical results indicate that the magnetic field strength has a significant effect on both the height and the interface morphology of the two-phase liquid mixed into. Specifically, the height of the mixed layer of the two-phase liquid decreases by applying a horizontally oriented magnetic field, and the higher the magnetic field strength, the lower the height of the mixed layer.

Figure 35 shows a comparative plot of the height of the mixed layer page after applying a forward magnetic field and a reverse magnetic field. The data shows that the horizontal direction of the magnetic field, whether it is forward or reverse magnetic field, has almost negligible effect on the height of the mixed layer of the two-phase liquid. Finally, the height of the mixed layer is investigated for We₁ = 10 000 and We₂ = 1000 with different Bo_m values (Bo_m = 636, 1434, and 2547), as shown in Fig. 36. The data show that as the We number increases, surface tension decreases, and the height of the mixed layer increases, indicating that surface tension inhibits Kelvin–Helmholtz instability.

---

## IV. CONCLUSIONS

In this paper, an adaptive-mesh-refinement (AMR) lattice Boltzmann method (LBM) is presented for two-phase ferrofluid flows to investigate interfacial dynamics under the influence of a magnetic field. In this study, the two-phase ferrofluid flow system is assumed to be incompressible, Newtonian and immiscible.

The phase interface is captured by a second order diffuse interface model and the hydrodynamics equation in the velocity-based form is used. For the non-conductive fluid, a scalar potential equation is adopted to describe the magnetic field. The discrete Boltzmann equations which are equivalent to three governing equations are solved using the Lax–Wendroff scheme. To reduce the computational overhead, the numerical schemes are implemented on a block-structured adaptive mesh. Moreover, the Neumann boundary conditions for the magnetic scalar potential are specialized within the Lax–Wendroff framework.

In the section of results and discussion, the simulation accuracy and efficiency of the AMR-LB model are tested by a circular cylinder in a uniform magnetic field, the deformation of a ferrofluid droplet, and the bubble rising in ferrofluid. Then, flow instabilities in two-phase ferrofluid systems are simulated, including the Rosensweig instability, the Rayleigh–Taylor instability, and the Kelvin–Helmholtz instability under different magnetic field configurations. The results show that the AMR-LB model can well reflect the rich interfacial dynamics of ferrofluids. For the Rosensweig instability, the adaptive-mesh-refinement method is verified or be very helpful to capture the correct wavelength relation. Moreover, the presence of a magnetic field is a very efficient factor in promoting Rayleigh–Taylor instability and suppresses the Kelvin–Helmholtz instability.

---

## ACKNOWLEDGMENTS

This work was supported by the National Natural Science Foundation of China (Grant Nos. 12172039, 12102228, and 11802159).

---

## AUTHOR DECLARATIONS

### Conflict of Interest

The authors have no conflicts to disclose.

### Author Contributions

**Zhenchao Guo:** Data curation (equal); Formal analysis (equal); Validation (equal); Writing – original draft (equal); Writing – review & editing (equal). **Shiting Zhang:** Formal analysis (equal); Validation (equal). **Yuqi Zhu:** Formal analysis (equal); Validation (equal). **Yang Hu:** Data curation (equal); Formal analysis (equal); Validation (equal); Writing – review & editing (equal). **Qiang He:** Formal analysis (equal); Validation (equal). **Xiaolong Yang:** Formal analysis (equal); Validation (equal). **Decai Li:** Formal analysis (equal); Validation (equal).

---

## DATA AVAILABILITY

The data that support the findings of this study are available within the article.

---

## REFERENCES

1. W. Lee, R. Scardovelli, A. Trubatch, and P. Yecko, "Numerical, experimental and theoretical investigation of bubble aggregation and deformation in magnetic fluids," Phys. Rev. E **82**, 016302 (2010).
2. R. E. Rosensweig, "Stress boundary-conditions in ferrohydrodynamics," Ind. Eng. Chem. Res. **46**(19), 6113–6117 (2007).
3. V. Segal, A. Rabinovich, D. Nattrass, K. Rag, and A. Nunes, "Experimental study of magnetic colloidal fluids behavior in power transformers," J. Magn. Magn. Mater. **215–216**, 513–515 (2000).
4. C. W. Ueberhuber, *Numerical Computation 1: Methods, Software, and Analysis* (Springer, 1997).
5. W. M. Yang, "A finite volume method for ferrohydrodynamic problems coupled with microscopic magnetization dynamics," Appl. Math. Comput. **441**(15), 127704 (2023).
6. A. Khan, X. D. Niu, Q. Z. Li, Y. Li, D. C. Li, and H. Yamaguchi, "Dynamic study of ferrodroplet and bubbles merging in ferrofluid by a simplified multiphase lattice Boltzmann method," J. Magn. Magn. Mater. **495**(1), 165869 (2020).
7. T. Huang, X. Liao, Z. Q. Huang, and R. Y. Wang, "Numerical simulation of ferrofluid flow in heterogeneous and fractured porous media based on finite element method," Front. Earth Sci. **9**, 693531 (2021).
8. Y. Hu, D. C. Li, and Q. He, "Generalized conservative phase field model and its lattice Boltzmann scheme for multicomponent multiphase flows," Int. J. Multiphase Flow **132**, 103432 (2020).
9. Q. Li, K. H. Luo, Q. J. Kang, Y. L. He, Q. Chen, and Q. Liu, "Lattice Boltzmann methods for multiphase flow and phase-change heat transfer," Prog. Energy Combust. Sci. **52**, 62–105 (2016).
10. S. T. Zhang, Y. Hu, Q. Li, D. C. Li, Q. He, and X. D. Niu, "A second-order phase field-lattice Boltzmann model with equation of state inputting for two-phase flow containing soluble surfactants," Phys. Fluids. **36**, 022104 (2024).
11. F. Bösch, S. S. Chikatamarla, and I. V. Karlin, "Entropic multirelaxation lattice Boltzmann models for turbulent flows," Phys. Rev. E **92**(4), 043309 (2015).
12. H. Yu, S. S. Girimaji, and L. S. Luo, "DNS and LES of decaying isotropic turbulence with and without frame rotation using lattice Boltzmann method," J. Comput. Phys. **209**(2), 599–616 (2005).
13. C. Y. Lim, C. Shu, X. D. Niu, and Y. T. Chew, "Application of lattice Boltzmann method to simulate microchannel flows," Phys. Fluids **14**(7), 2299–2308 (2002).
14. S. Ansumali, I. V. Karlin, C. E. Frouzakis, and K. B. Boulouchos, "Entropic lattice Boltzmann method for microflows," Physica A **359**, 289–305 (2006).
15. Y. Hu, D. C. Li, S. Shu, and X. D. Niu, "Modified momentum exchange method for fluid-particle interactions in the lattice Boltzmann method," Phys. Rev. E **91**(3), 033301 (2015).
16. Y. Hu, H. Yuan, S. Shu, X. D. Niu, and M. Li, "An improved momentum exchanged-based immersed boundary–lattice Boltzmann method by using an iterative technique," Comput. Math. Appl. **68**(3), 140–155 (2014).
17. Y. Hu, D. C. Li, S. Shu, and X. D. Niu, "Finite-volume method with lattice Boltzmann flux scheme for incompressible porous media flow at the representative-elementary-volume scale," Phys. Rev. E **93**(2), 023308 (2016).
18. Y. Hu, D. C. Li, S. Shu, and X. D. Niu, "A multiple-relaxation-time lattice Boltzmann model for the flow and heat transfer in a hydrodynamically and thermally anisotropic porous medium," Int. J. Heat Mass Transfer **104**, 544–558 (2017).
19. N. I. Prasianakis and I. V. Karlin, "Lattice Boltzmann method for thermal flow simulation on standard lattices," Phys. Rev. E **76**(1), 016702 (2007).
20. Y. Hu, D. C. Li, S. Shu, and X. D. Niu, "Study of multiple steady solutions for the 2 D natural convection in a concentric horizontal annulus with a constant heat flux wall using immersed boundary-lattice Boltzmann method," Int. J. Heat Mass Transfer **81**, 591–601 (2015).
21. Y. Hu, D. C. Li, S. Shu, and X. D. Niu, "Full Eulerian lattice Boltzmann model for conjugate heat transfer," Phys. Rev. E **92**(6), 063305 (2015).
22. B. An and W. M. Sang, "The numerical study of lattice Boltzmann method based on different grid structure," Chin. J. Theor. Appl. Mech. **45**(5), 699–706 (2013).
23. Y. Hu, D. C. Li, and X. D. Niu, "Phase-field-based lattice Boltzmann model for multiphase ferrofluid flows," Phys. Rev. E **98**(3), 033301 (2018).
24. X. Li, P. Yu, X. D. Niu, D. C. Li, and H. Yamaguchi, "A magnetic field coupling lattice Boltzmann model and its application on the merging process of multiple-ferrofluid-droplet system," J. Appl. Math. Comput. **393**, 125769 (2021).
25. Y. Li, X. D. Niu, A. Khan, D. C. Li, and H. Yamaguchi, "A numerical investigation of dynamics of bubbly flow in a ferrofluid by a self-correcting procedure-based lattice Boltzmann flux solver," Phys. Fluids **31**(8), 082107 (2019).
26. A. Khan, S. T. Zhang, Q. P. Li, H. Zhang, Y. Q. Wang, and X. D. Niu, "Wetting dynamics of a sessile ferrofluid droplet on solid substrates with different wettabilities," Phys. Fluids **33**(4), 042115 (2021).
27. S. T. Zhang, X. D. Niu, Q. P. Li, A. Khan, Y. Hu, and D. C. Li, "A numerical investigation on the deformation of ferrofluid droplets," Phys. Fluids **35**(1), 012102 (2023).
28. Y. Huang, Z. Ke, Z. Li, Y. Gao, Z. Tang, and Y. Zhang, "A non-uniform magnetic field coupled lattice Boltzmann model and its application on the wetting dynamics of a ferrofluid droplet under gravity effects," Comput. Math. Appl. **143**, 73–93 (2023).
29. Q. He, W. F. Huang, J. J. Xu, Y. Hu, and D. C. Li, "A hybrid immersed interface and phase-field-based lattice Boltzmann method for multiphase ferrofluid flow," Comput. Fluids **255**, 105821 (2023).
30. J. Tölke, S. Freudiger, and M. Krafczyk, "An adaptive scheme using hierarchical grids for lattice Boltzmann multi-phase flow simulations," Comput. Fluids **35**, 820–830 (2006).
31. H. B. Huang, M. C. Sukop, and X. Y. Lu, *Multiphase Lattice Boltzmann Methods: Theory and Application* (John Wiley and Sons, 2015).
32. Z. Yu and L. S. Fan, "An interaction potential based lattice Boltzmann method with adaptive mesh refinement (AMR) for two-phase flow simulation," J. Comput. Phys. **228**(17), 6456–6478 (2009).
33. Y. Chen, Q. Kang, Q. D. Cai, and D. X. Zhang, "Lattice Boltzmann method on quadtree grids," Phys. Rev. E **83**(2), 026707 (2011).
34. Z. L. Liu, F. B. Tian, and X. Y. Feng, "An efficient geometry-adaptive mesh refinement framework and its application in the immersed boundary lattice Boltzmann method," Comput. Methods Appl. Mech. Eng. **392**, 114662 (2022).
35. Y. Hasegawa, T. Aoki, H. Kobayashi, Y. Idomura, and N. Onodera, "Tree cutting approach for domain partitioning on forest-of-octrees-based block-structured static adaptive mesh refinement with lattice Boltzmann method," Parallel Comput. **108**, 102851 (2021).
36. R. Deiterding and S. L. Wood, "Predictive wind turbine simulation with an adaptive lattice Boltzmann method for moving boundaries," J. Phys.: Conf. Ser. **753**, 082005 (2016).
37. A. Fakhari, M. Geier, and T. Lee, "A mass-conserving lattice Boltzmann method with dynamic grid refinement for immiscible two-phase flows," J. Comput. Phys. **315**, 434–457 (2016).
38. A. Fakhari and T. Lee, "Finite-difference lattice Boltzmann method with a block-structured adaptive-mesh-refinement technique," Phys. Rev. E **89**(3), 033310 (2014).
39. Q. He, W. F. Huang, Y. Yin, Y. Hu, and D. C. Li, "A mass-conserving and volume-preserving lattice Boltzmann method with dynamic grid refinement for immiscible ternary flows," Phys. Fluids **34**, 093321 (2022).
40. Y. C. Xia, B. W. Yao, K. Wang, and Z. Y. Li, "A three-dimensional fully threaded tree adaptive mesh phase-field lattice Boltzmann method for gas–liquid phase change problems," Phys. Fluids **35**, 103323 (2023).
41. R. E. Rosensweig, *Ferrohydrodynamics* (Cambridge University Press, 1985).
42. P. H. Chiu and Y. T. Lin, "A conservative phase field method for solving incompressible two-phase flows," J. Comput. Phys. **230**, 185–204 (2011).
43. L. Li, R. Mei, and J. F. Klausner, "Lattice Boltzmann models for the convection-diffusion equation: D 2 Q 5 vs D 2 Q 9," Int. J. Heat Mass Transfer **108**, 41–62 (2017).
44. A. Fakhari, D. Bolster, and L. S. Luo, "A weighted multiple-relaxation-time lattice Boltzmann method for multiphase flows and its application to partial coalescence cascades," J. Comput. Phys. **341**, 22–43 (2017).
45. J. Pipper, Y. Zhang, P. Neuzil, and T. M. Hsieh, "Clockwork pcr including sample preparation," Angew. Chem., Int. Ed. **47**, 3900–3904 (2008).
46. L. Hajba and L. Guttman, "Circulating tumor-cell detection and capture using microfluidic devices," Trends Analyt. Chem. **59**, 9–16 (2014).
47. P. K. Yuen, L. J. Kricka, P. M. Fortina, N. J. Panaro, T. Sakazume, and P. Wilding, "Microchip module for blood sample preparation and nucleic acid amplification reactions," Genome Res. **11**(3), 405–412 (2001).
48. C. Flament, S. Lacis, J. C. Bacri, A. Cebers, S. Neveu, and R. Perzynski, "Measurements of ferrofluid surface tension in confined geometry," Phys. Rev. E **53**(5), 4801–4806 (1996).
49. Y. Hu, D. C. Li, and L. Jin, "Hybrid Allen-Cahn-based lattice Boltzmann model for incompressible two-phase flow: The reduction of numerical dispersion," Phys. Rev. E **99**, 023302 (2019).
50. X. Li, Z. Q. Dong, F. Li, L. P. Wang, X. D. Niu, H. Yamaguchi, D. C. Li, and P. Yu, "A fractional-step lattice Boltzmann method for multiphase flows with complex interfacial behavior and large density contrast," Int. J. Multiphase Flow **149**, 103982 (2022).
51. S. Aland and A. Voigt, "Benchmark computations of diffuse interface models for two-dimensional bubble dynamics," Int. J. Numer. Methods Fluids **69**, 747 (2012).
52. X. D. Niu, A. Khan, Y. Ouyang, M. F. Chen, D. C. Li, and H. Yamaguchi, "A simplified phase-field lattice Boltzmann method with a self-corrected magnetic field for the evolution of spike structures in ferrofluids," Appl. Math. Comput. **436**, 127503 (2023).
53. P. F. Yuan, Q. X. Cheng, Y. Hu, Q. He, W. F. Huang, and D. C. Li, "Phase-field-based finite element model for two-phase ferrofluid flows," Phys. Fluids **36**, 022016 (2024).
54. V. Zaitsev and M. I. Shliomis, "Nature of the instability of the interface between two liquids in a constant field," Dokl. Phys. **14**, 1001 (1970).
55. M. S. Krakov, A. R. Zakinyan, and A. A. Zakinyan, "Instability of the miscible magnetic/non-magnetic fluid interface," J. Fluid Mech **913**, A 30 (2021).
56. M. D. Cowley and R. E. Rosensweig, "The interfacial stability of a ferromagnetic fluid," J. Fluid Mech. **30**, 671–688 (1967).
57. R. R. Nourgaliev, T. N. Dinh, and T. G. Theofanous, "A pseudocompressibility method for the numerical simulation of incompressible multifluid flows," Int. J. Multiphase Flow **30**, 901–937 (2004).
58. Y. Wang, C. Shu, and C. J. Teo, "Development of LBGK and incompressible LBGK-based lattice Boltzmann flux solvers for simulation of incompressible flows," Numer. Methods Fluids. **75**, 344–364 (2014).
59. Y. Wang, C. Shu, J. Y. Shao, J. Wu, and X. D. Niu, "A mass-conserved diffuse interface method and its application for incompressible multiphase flows with large density ratio," J. Comput. Phys. **290**, 336–351 (2015).
60. A. Völkel, A. Kögel, and R. Richter, "Measuring the Kelvin-Helmholtz instability, stabilized by a tangential magnetic field," J. Magn. Magn. Mater. **505**, 166693 (2020).
61. H. G. Lee and J. Kim, "Two-dimensional Kelvin–Helmholtz instabilities of multi-component fluids," Eur. J. Mech. B Fluids **49**, 77–88 (2015).
