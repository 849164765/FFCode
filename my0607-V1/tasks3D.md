首先需要澄清，该文献中使用的是**D3Q27**离散速度模型，而非您提到的D3Q19。D3Q27模型虽然计算成本更高，但能提供更好的数值稳定性和精度，尤其适合模拟自密实混凝土（SCC）这类非牛顿宾汉姆流体的复杂流动。

以下是该论文中3D MRT-LBM的核心理论公式与基本原理总结：

### 🧱 1. 理论基础：从玻尔兹曼方程到格子玻尔兹曼方法 (LBM)

LBM的核心思想是架起一座从微观/介观粒子行为到宏观流体运动的桥梁。它通过对连续**玻尔兹曼方程**进行时间、空间和速度的离散，得到**格子玻尔兹曼方程** (LBE)，最终可通过Chapman-Enskog展开，严格推导出宏观的**纳维-斯托克斯 (Navier-Stokes) 方程**。

### 💎 2. 离散速度模型：D3Q27

文献选择了D3Q27模型，即三维空间内包含27个离散速度方向。其离散速度 $\mathbf{e}\_{\alpha}$ 的具体定义如下：

\[
\mathbf{e}\_{\alpha} = \left{ \begin{array}{ll}
(0,0,0)c, & \alpha = 0, \\
(\pm 1,0,0)c, (0,\pm 1,0)c, (0,0,\pm 1)c, & \alpha = 1,\ldots ,6, \\
(\pm 1,\pm 1,0)c, (\pm 1,0,\pm 1)c, (0,\pm 1,\pm 1)c, & \alpha = 7,\ldots ,18, \\
(\pm 1,\pm 1,\pm 1)c, & \alpha = 19,\ldots ,26,
\end{array} \right.
]

其中，$c = \delta x / \delta t$ 为格子速度，$\delta x$ 和 $\delta t$ 分别为格子间距和时间步长。

### ⚙️ 3. 核心方程：MRT-LBM

多松弛时间 (MRT) 模型的核心是通过一个碰撞矩阵 $\Lambda$，对不同阶的矩（如密度、动量、应力张量分量）采用独立的松弛时间，从而提供了比单松弛时间 (SRT/BGK) 模型更高的数值稳定性和灵活性。

#### a. **演化方程**

流体运动的基本方程为**格子玻尔兹曼方程 (Lattice Boltzmann Equation)**：

\[
f\_{\alpha}(\mathbf{x} + \mathbf{e}*{\alpha}\delta t, t + \delta t) - f*{\alpha}(\mathbf{x}, t) = \Lambda\_{\alpha j}\Big\[f\_{j}(\mathbf{x}, t) - f\_{j}^{\mathrm{eq}}(\mathbf{x}, t)\Big] + F\_{\alpha}\delta t
]

- $f\_{\alpha}$：沿第 $\alpha$ 个方向的粒子分布函数。
- 方程左侧描述了粒子在时间 $\delta t$ 内的**迁移 (Streaming)** 过程。
- 方程右侧描述了粒子间的**碰撞 (Collision)** 过程，在MRT框架下，碰撞矩阵 $\Lambda$ 表示为 $\Lambda = -\mathbf{M}^{-1}\mathbf{S}\mathbf{M}$，其中 $\mathbf{M}$ 是将分布函数转换到矩空间的变换矩阵，$\mathbf{S}$ 是一个对角松弛矩阵【15†L3-L5】。

#### b. **平衡态分布函数**

平衡态分布函数 $f\_{\alpha}^{\mathrm{eq}}$ 是从 **Maxwell-Boltzmann 分布**经 Taylor 级数展开至二阶得到：

\[
f\_{\alpha}^{\mathrm{eq}} = \rho w\_{\alpha}\left\[1 + \frac{\mathbf{e}*{\alpha}\cdot\mathbf{u}}{c*{s}^{2}} + \frac{\left(\mathbf{e}*{\alpha}\cdot\mathbf{u}\right)^{2} - \left(c*{s}|\mathbf{u}|\right)^{2}}{2c\_{s}^{4}}\right]
]

- 其中，声速 $c\_{s} = c / \sqrt{3}$。

#### c. **外力项**

为了模拟重力等外部作用，文献通过在演化方程中增加一个外力项 $F\_{\alpha}$ 来实现：

\[
F\_{\alpha} = w\_{\alpha}\rho \frac{\mathbf{e}*{\alpha}\cdot\mathbf{a}}{c*{s}^{2}}
]

- 其中，$\mathbf{a}$ 代表加速度（如重力加速度）。

#### d. **宏观物理量**

流体的宏观密度 $\rho$ 和动量 $\rho \mathbf{u}$ 由各方向的分布函数统计得出：

- **密度**: $\rho = \sum\_{\alpha = 0}^{26} f\_{\alpha}$
- **动量**: $\rho\mathbf{u} = \sum\_{\alpha = 0}^{26} \mathbf{e}*{\alpha}f*{\alpha} + \frac{1}{2}\mathbf{a}\delta t$
- **压力**: $p = \rho c\_{s}^{2}$

### 🌊 4. 非牛顿流体模型

为精确描述SCC这类非牛顿流体，论文采用了**修正的 Herschel-Bulkley 模型**。该模型能有效处理流体在低剪切率下的粘度发散问题，既能捕捉存在屈服应力（$\sigma\_y$）时的类固体行为，也能模拟屈服后的幂律流体特性（当 $n=1$ 时退化为Bingham模型）。模型的具体形式为：

\[
\mu =
\begin{cases}
\mu\_{y}, & \sigma < \sigma\_{y} \\
\frac{\sigma\_{y}}{\dot{\gamma}} + \frac{k}{\dot{\gamma}}\left\[\dot{\gamma}^{n} - \left(\frac{\sigma\_{y}}{\mu\_{y}}\right)^{n}\right], & \sigma \geq \sigma\_{y}
\end{cases}
]

***

以上便是Qiu和Han (2018)在其论文中构建3D MRT-LBM模型所依据的核心公式与理论基础。希望这份梳理对您有帮助。如果对公式中的符号含义或特定推导过程有疑问，我可以继续为您提供更详细的解释。
