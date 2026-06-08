# 基于VOF方法的同轴气泡聚并数值模拟研究

郝文杰，韩晋玉，刘　瑶，赵陈儒*，薄涵亮 

（清华大学 核能与新能源技术研究院，北京　100084） 

摘要：气液两相流动广泛存在于反应堆堆芯以及蒸汽发生器等设备中。气泡聚并、破碎等气泡动力学行 为直接影响相关设备的热工水力性能。群平衡模型（PBM）是多相流系统中处理多组分粒子平衡问题的重 要方法。气泡聚并核函数是 PBM主要核函数之一，用于描述两个或多个气泡在碰撞过程中聚并为一个较 大气泡的概率，其准确性直接影响欧拉-欧拉框架下的气泡动力学行为研究及热工水力性能的预测。排液 模型是当前应用最为广泛的聚并效率模型，是聚并核函数的重要部分。通过将流体体积法（VOF）与连续 表面力模型（CSF）相结合对两同轴气泡聚并过程开展了 2D数值模拟，并根据模拟结果定量获得气泡聚并 过程的排液时间，通过研究不同黏度下气泡的聚并，定量描述了气泡开始排液的接触时刻，讨论了初始液 膜厚度与气泡直径、形变比以及黏度之间的关系，对文献中已有的 4个排液时间模型的适用性及准确性 进行系统评价。研究结果表明，接触时刻的定义，初始液膜厚度的预测，以及黏性力和惯性力作用的综合 考虑和评估对于提高排液时间模型的准确性，扩展其适用范围，进而提高 PBM方法的计算准确性有重要 的意义。 

关键词：VOF方法；气泡聚并；接触时刻；初始液膜厚度；排液时间 

中图分类号：TL33 

文献标志码：A 

文章编号：1000-6931（2025）06-1262-10 

doi：10.7538/yzk.2024.youxian.0776 

# Numerical Simulation Study of Coaxial Bubble Coalescence Based on VOF Method

HAO Wenjie,  HAN Jinyu,  LIU Yao,  ZHAO Chenru*,  BO Hanliang 

(Institute of Nuclear and New Energy Technology, Tsinghua University, Beijing 100084, China) 

Abstract: Gas-liquid  two-phase  flow  is  widely  found  in  reactor  cores  and  equipments  such  as  steam generators.  The  bubble  kinetic  behaviors,  such  as  bubble  coalescence  and  breakup,  directly  affect  the thermal-hydraulic  performance  of  the  related  equipment.  Bubble  coalescence  can  affect  flow characteristics  such  as  turbulence  intensity,  velocity  distribution,  and  system  pressure  drop.  Bubble coalescence can also lead to the occurrence of flow instability, which must be avoided during reactor operation.  It  is  necessary  to  consider  the  influence  of  bubble  coalescence  and  other  behaviors  in  the process of designing heat exchange equipment and systems. The population balance model (PBM) is an important  method  to  deal  with  the  equilibrium  problem  of  multi-component  particles  in  multi-phase flow  systems.  The  bubble  coalescence  kernel  function  is  one  of  the  main  kernel  functions  of  PBM, which  is  used  to  describe  the  probability  that  two  or  more  bubbles  will  coalesce  into  a  larger  bubble during collision, and its accuracy directly affects the study of bubble kinetic behavior and the prediction of  thermal-hydraulic  performance  under  the  Euler-Euler  frameworks.  The  film  drainage  model  is  the most widely used coalescence efficiency model, which is an important part of the coalescence kernel function.  The  definition  of  bubble  contact  time  and  the  quantification  of  initial  liquid  film  thickness have not been clearly pointed out in previous studies. The influence of bubble diameter, and continuous phase viscosity on the initial liquid film thickness and liquid drainage time is still unclear. In this paper, the 2D numerical simulation of two coaxial bubbles coalescence process was carried out by the volume of fluid method (VOF) combined with the continuous surface force model (CSF) by ANSYS_FLUENT, and the drainage time of the bubble coalescence process under different viscosities was quantitatively obtained according to the simulation results. The contact moment when the bubbles started to discharge was  quantitatively  described,  and  the  relationship  between  the  initial  film  thickness  and  the  bubble diameter, deformation ratio, and viscosity was meticulously discussed. The applicability and accuracy of the four existing drainage time models in the literature were systematically evaluated. The moment corresponding to the turning point of the relative velocity of the two bubble films was defined as the contact moment. It is found that the initial liquid film thickness and drainage time are affected by the bubble diameter and liquid viscosity. The results show that the definition of the contact moment and the comprehensive  consideration  and  evaluation  of  the  effects  of  viscous  and  inertial  forces  are  of  great significance to improve the accuracy of the drainage time model, extend its scope of application, and thus improve the calculation accuracy of the PBM method. 

Key words: VOF method; bubble coalescence; contact moment; initial film thickness; drainage time 

气泡聚并、破碎等气泡动力学行为研究对于 反应堆堆芯、蒸汽发生器等换热设备内两相流动 沸腾换热特性有重要的意义。气泡聚并、破碎直 接影响气相分布和气液两相界面面积，进而影响 流道截面含气率及换热特性。气泡聚并对两相流 动及换热的影响是一个复杂的过程。首先，气泡 聚并会影响湍流强度、速度分布和系统压降等流 动特性，气泡聚并的发生也可能会导致流动不稳 定性的发生[ 1 ]，这在反应堆运行过程中必须避 免。其次，对壁面换热而言，气泡聚并可能导致局 部热流密度的波动，进而引起壁面温度的波动，影 响换热性能的同时可能出现传热恶化[2]；对整个 流道而言，气泡聚并会使气泡数量减少，可能导致 沸腾传热效率的降低。因此，在设计换热设备及 系统的过程中，必须考虑气泡聚并等行为的影响， 研究气泡聚并过程对优化模拟计算精度有一定的 意义。 

群平衡模型（PBM）用于描述离散相在连续相 中的大小、分布，并预测其随时间的变化，包括粒 子破碎、聚并、生长、沉积等过程造成的变化。气 泡聚并模型是 PBM 方法中重要的子模型，提高其 预测的准确性对提高气液两相流数值模拟结果的 准确性具有非常重要的意义。气泡聚并频率 $\Gamma ( d _ { 1 } , d _ { 2 } )$ 用于描述气泡之间发生碰撞并且聚并的 概率，定义为碰撞频率 $h ( d _ { 1 } , d _ { 2 } )$ 与聚并效率 $\lambda ( d _ { 1 } , d _ { 2 } )$ di的乘积， （i=1，2）为气泡直径。能量模型、临界 相对速度模型和排液模型是描述聚并效率的 3 种 常见模型，其中，排液模型应用最广。Coulaloglou 和 T a v l a r i d e s [ 3 ] 给 出 的 排 液 模 型 聚 并 效 率 形 式为： $\lambda ( d _ { \mathrm { 1 } } , d _ { \mathrm { 2 } } ) = \mathrm { e x p } ( - t _ { \mathrm { d } } / \tau _ { \mathrm { c } } )$ τc（ 为接触时间， $t _ { \mathrm { d } }$ 为 排液时间）。该公式目前广泛用于聚并效率的计 算[4-5]，但其合理性及理论依据仍有待进一步验 证[6]。接触时间定义为两个气泡从接触开始（液膜 形成）到脱离（或反弹）的时间间隔，接触时两气泡 间液体的厚度称为初始液膜厚度 $\left( \begin{array} { l } { \boldsymbol { h } _ { 0 } } \end{array} \right)$ ，液膜形成 后的液膜厚度（h）随时间不断减薄。排液时间定 义为两个气泡从接触开始到液膜达到临界厚度 $\left( \boldsymbol { h _ { \mathrm { c r } } } \right)$ 而破裂的时间间隔。 $h _ { 0 } .$ 、初始排液速度（ $U _ { \mathrm { r e l } , 0 } )$ 、 平均薄膜半径 $\left( \boldsymbol { r } _ { \mathrm { f } , 0 } ^ { } \right)$ ）作为初始条件是求解排液时 间的核 $\therefore \dot { \mathrm { L } } ^ { [ 7 - 8 ] }$ ， 其中，接触时刻的准确定义是定量 描述 $h _ { 0 }$ 、 $U _ { \mathrm { r e l , 0 } }$ 和 ${ r } _ { \mathrm { f } , 0 }$ 的前提，也是获得接触时间和 排液时间的关键。但目前对接触时刻的定量讨论 仍非常缺乏：Chesters[9] 将液膜变平或 $F = 2 \sigma / R$ σ（F 为压力差，Pa； 为表面张力系数，N/m；R 为气 泡半径，m）看作是排液开始；李昕晨[10] 将高黏度 甘油溶液中垂直同轴上升气泡聚并过程的接触时 刻定义为下气泡顶点与上气泡尾部处于同一水平 线的时刻；Orvalho 等[11] 将垂直平行的毛细管产生 的气泡聚并过程中的接触时刻定义为两气泡恰好 第 1 次物理接触的时刻，在其实验中定义为气泡 间距离达到相机空间分辨率 1 像素（为 5 µm）的时 刻；Nguyen 等[12] 将天然硅酸盐熔体（极高黏度）中 单气泡与界面的接触时刻定义为在界面附近气泡 速度趋近于零的时刻；Kirkpatrick 和 Lockett[13]、 Sanada 等[14] 认为单气泡与液面接触之前，气泡的 减速是不明显的，因此将气泡速度趋于稳定的时 刻定义为接触时刻。但已有研究均未对其定义的 接触时刻所对应的初始液膜厚度进行研究。 

许多研究者对初始液膜厚度和临界液膜厚度 进行了估计：Coulaloglou 和 Tavlarides [3] 将初始液 膜厚度和临界液膜厚度看作常数；Das 等[15] 认为 不可变形液滴间临界液膜厚度对给定系统是特定 的；Lee 等[8] 认为初始液膜厚度在 0.01～0.1 mm 之 间，临界液膜厚度在 0.01～0.1 μm 之间；Kirkpatrick 和 Lockett[13] 估计空气-水系统中初始液膜厚度 为 0.1 mm，临界液膜厚度为 2.5 nm；龚升高[16] 测 试了初始液膜厚度对临界速度预测的影响，选取 $h _ { 0 } = 0 . 2 d$ ，并验证了临界液膜厚度为 50 nm 的合理 性；闻昭权[17] 根据龚升高模型预测的结果，选取 $h _ { 0 } = 0 . 0 5 d$ 作为初始液膜厚度，同样，Song 等[18] 也 0.05d将初始液膜厚度选为 。 

综上所述，已有研究对于气泡接触时刻的定 义尚不明确，且对初始液膜厚度多为估计；气泡 直径、连续相黏度对于初始液膜厚度、排液时间 的影响尚不明确。另一方面，VOF、Level-Set（LS）、 CLSVOF 等界面构造方法可捕捉气液两相界面细 节以及捕获和解析详细信息（包括空间和时间信 息）的能力已经被证明适合于描述气泡间的聚并 等相界面复杂变形的问题[19-21]。Abbassi 等[22] 采 用 VOF 方法研究了不同黏度静态液体中单气泡 上升以及连续两个气泡的同轴聚并，获得了单气 泡尾流结构、液相速度场以及同轴气泡的形状演 变、上升速度、聚并特征。不同运动方向气泡之 间的相互作用是一个复杂的过程，而且同轴气泡 聚并在鼓泡塔反应器等设备中是一个典型的现 象。相较于并排、斜向等形式的气泡聚并，同轴 气泡聚并的研究更为简便。选取同轴气泡聚并 这一简单形式，排除了气泡滑移等情况的影响， 也便于获取气泡界面间液体厚度、气泡界面运动 速度等参数。对同轴气泡的聚并过程进行定量 研究一定程度上弥补了当前在气泡聚并过程中 定量研究的短缺，对于优化聚并模型具有重要的 意义。因此，本文拟采用 VOF 方法对同轴双气泡 聚并开展数值模拟研究，获得气泡界面运动速 度、界面顶点间液体厚度变化等细节信息，得到 可用于定量描述接触时刻的参数并对接触时刻 进行明确定义，进而研究气泡直径、连续相黏度 对同轴气泡聚并的影响，建立其与初始液膜厚度 之间的联系。 

## 1 模型及验证

## 1.1 物理模型

本文采用 Ansys_Fluent 数值模拟软件对气泡 槽中的两同轴气泡聚并过程进行 2D 瞬态数值模 g拟。气泡槽几何模型如图 1 所示，重力加速度（ ） 为 $9 . 8 1 \mathrm { m } / \mathrm { s } ^ { 2 } ,$ ，方向为 y 轴负向。为减小壁面效应对 气泡上升行为的影响，竖直两壁面之间的距离设 置为 10 倍初始气泡直径，上气泡下顶点与下气泡 上顶点间初始间距保持 2 mm，下气泡质心距液面 59 mm，物性参数列于表1。 

![](images/d2e0959749c88896b2de230da34323859349699ca2285bb01e23514250655c1d.jpg)



图 1    气泡槽几何模型



Fig. 1    Geometry of bubble trough



表 1 物性参数



Table 1 Physical property parameter


<table><tr><td>物性参数</td><td>数值</td></tr><tr><td>离散相密度, kg/m3</td><td>1.225</td></tr><tr><td>离散相黏度, Pa·s</td><td><eq>1.7894 \times 10^{-5}</eq></td></tr><tr><td>连续相密度, kg/m3</td><td>998.2</td></tr><tr><td>表面张力系数, N/m</td><td>0.072</td></tr></table>

两侧及底部壁面采用壁面黏附无滑移边界条 件，出口为 1 个标准大气压（101.325 kPa）的压力 出口边界。保持表面张力系数为 0.072 N/m，研究 黏度对同轴气泡聚并的影响。初始时刻放置两同 轴等直径气泡，直径分别为 4、5、6、7、8 mm，连 续相黏度分别为 0.010 03、0.015 045、0.020 06、 0.025 075、0.050 15、0.075 225 Pa·s，共 30 个工 况。黏度分别为 10 倍、15 倍、20 倍、25 倍、 50 倍、75 倍纯净水的黏度。在本研究中，将黏度 0.025 075 Pa·s 看作较低黏度和较高黏度的边界， 黏度小于 0.025 075 Pa·s 看作是较低黏度，黏度大 于 0.025 075 Pa·s 看作是较高黏度。 

## 1.2 控制方程

α α = 1VOF 方法用相体积分数（ ）来表征气（ ）、 α = 0 液（ ）及气液界面 $( 0 { < } \alpha { < } 1 )$ ，通过对流方程 求解出体积函数后，再根据求出的体积函数值构 造和追踪运动相界面的形状。 

连续性方程： 

$$
\sum_ {i = 1} ^ {2} \alpha_ {i} = 1 \tag {1}
$$

$$
\frac {\partial \alpha_ {i}}{\partial t} + \nabla \cdot (\boldsymbol {u} \alpha_ {i}) = \frac {S _ {\mathrm{m} , i}}{\rho_ {i}} \tag {2}
$$

式中： $S _ { \mathrm { m } , i }$ 为体积质量源相； $\rho _ { i }$ 为相密度。界面处 ρ每个单元的物性为体积平均物性（以密度 和黏度 µ为例）： 

$$
\rho = \alpha_ {\mathrm{g}} \rho_ {\mathrm{g}} + \rho_ {\mathrm{l}} (1 - \alpha_ {\mathrm{g}}) \tag {3}
$$

$$
\mu = \alpha_ {\mathrm{g}} \mu_ {\mathrm{g}} + \mu_ {\mathrm{l}} (1 - \alpha_ {\mathrm{g}}) \tag {4}
$$

在 VOF 方法中，求解单个动量方程和能量方 程，得到的速度场和温度场在相之间共享。 

动量方程： 

$$
\frac {\partial (\rho \boldsymbol {u})}{\partial t} + \nabla \cdot (\rho \boldsymbol {u u}) = - \nabla \boldsymbol {p} + \nabla \cdot \mu (\nabla \boldsymbol {u} + \nabla \boldsymbol {u} ^ {\mathrm{T}}) + \rho \boldsymbol {g} + \boldsymbol {F} _ {\mathrm{csf}} \tag {5}
$$

u p式中： 为速度； 为压力；Brackbill[23] 提出的连续 表面张力模型（CSF）下的表面张力 $\pmb { F } _ { \mathrm { c s f } }$ 为： 

$$
\boldsymbol {F} _ {\mathrm{csf}} = \sigma \frac {\alpha_ {1} \rho_ {1} \kappa_ {1} \nabla \alpha_ {\mathrm{g}} + \alpha_ {\mathrm{g}} \rho_ {\mathrm{g}} \kappa_ {\mathrm{g}} \nabla \alpha_ {1}}{0 . 5 (\rho_ {\mathrm{g}} + \rho_ {1})} \tag {6}
$$

界面曲率： 

$$
\kappa = \nabla \cdot \frac {\nabla \alpha}{| \nabla \alpha |} \tag {7}
$$

## 1.3 求解器设置

瞬态计算采用 PISO 压力速度耦合方案，采用 PRESTO！空间离散压力插值方法、QUICK 动量 离散插值方法，显式时间积分格式搭配 Geo-Reconstruct 体积分数几何重构方法被用于构造气液 界面，瞬态方程通过一阶隐式时间项离散（first order implicit）方法求解。 

## 1.4 网格无关性验证

为验证网格无关性，以 Liu 等[ 2 4 ] 的直径为 2.01 mm 的单气泡上升实验数据为参考，物性参数 为：液相，水；密度，999.8 kg/m3；黏度，0.001 38 Pa·s； 表面张力系数，0.074 N/m。在此物性条件下，以尺 寸分别为 0.1、0.05、0.025 mm 的网格对单气泡上 升进行数值模拟，网格数量分别为 28 万、112 万 和448万，网格正交质量平均值均为1。 

不同形状气泡的高度方向投影长度（W）及长 度方向投影长度（L）如图 2 所示，形变比（E）定义 为 $E { = } W / L _ { c }$ 。如图 3 所示，不同网格尺寸下气泡上 升过程中气泡质心距离起始点不同高度（Z）处的 形变比、速度（V）与实验数据整体趋势保持一 致。由于实验中气泡由喷嘴排出，气泡在脱离喷 嘴之初会有小幅度震荡，造成气泡初始阶段的形 变比震荡，这是实验与数值模拟在气泡上升初期 相异的原因。 

![](images/0a9134131e93115bd6dac989afa0ed8de326fa47ef8a9e41e67b135700068b95.jpg)



图 2    不同形状气泡高度、长度方向投影



Fig. 2    Projections of bubble with different shapes in direction of height and length


![](images/85202dc6e7661e3c6d67e539835e1babdce73beac20baf6f157df2d33f05d2af.jpg)


![](images/52a1159b72861825e638e226b144b8119d8a582cc52f3ddee6d1ea16b6e8c66a.jpg)



图 3    直径2.01 mm气泡形变比和浮升速度随上升高度的变化



Fig. 3    Deformation ratio and velocity of bubble with a diameter of 2.01 mm vs. rising height


气泡上升 20 mm 之后形变比趋于稳定，以 20 mm 之后形变比平均值为参考 ， 实验值 为 0.429，计算结果列于表 2，网格尺寸为 0.05 mm 的 气泡平均形变比与实验数据的相对误差为 5.35%， 网格尺寸为 0.025 mm 的气泡平均形变比与实验 数据的相对误差为 0.17%；气泡上升 15 mm 之后 速度趋于稳定，以 15 mm 之后速度平均值为参考， 实验值为 0.302 8 m/s，计算结果列于表 2，网格尺 寸为 0.05 mm 的气泡平均速度与实验数据的相 对误差为 9.42%，网格尺寸为 0.025 mm 的气泡平 均速度与实验数据的相对误差为 9.74%。综合考 虑计算时间成本及计算准确性，本文采用尺寸为 0.05 mm 的网格进行后续同轴双气泡聚并的数值 模拟。 


表 2 不同网格尺寸气泡稳定平均形变比、速度与实验数据相对误差



Table 2 Relative errors of average deformation ratio and velocity at different grid sizes combined with experiment data


<table><tr><td>网格尺寸/mm</td><td>网格数量/万</td><td>形变比</td><td>形变比相对误差/%</td><td>速度/(m/s)</td><td>速度相对误差/%</td></tr><tr><td>0.1</td><td>28</td><td>0.534</td><td>24.44</td><td>0.264</td><td>12.82</td></tr><tr><td>0.05</td><td>112</td><td>0.452</td><td>5.35</td><td>0.2743</td><td>9.42</td></tr><tr><td>0.025</td><td>448</td><td>0.428</td><td>0.17</td><td>0.2733</td><td>9.74</td></tr></table>

## 1.5 模型验证

本文采用 VOF 方法对 Liu 等[24] 直径为 2.01 mm 单气泡上升过程进行数值模拟，获得气泡形变比 E、垂直上升速度 $V \sb { \mathrm { h } } \sb { \textnormal { \scriptsize h } }$ 、水平速度 $V _ { \mathbf { v } } ,$ 并与相应的实 验结果对比，如图 4 所示。气泡垂直上升速度、形 变比在气泡上升 50 mm 之前符合较好，10～50 mm 间相对稳定段的气泡垂直上升速度平均值相对误 差为 9.12%，气泡形变比平均值相对误差为 3.48%。 

数值模拟结果中，气泡上升 30 mm 之后横向 路径 $V _ { \mathrm { v } }$ 出现周期振荡，这一现象主要与气泡上升 路径不稳定性、气泡受力的不稳定性相关，如 $\mathrm { W u }$ 等[25] 表明直径大于 1.5 mm 的气泡，其中保持 球形的气泡沿之字形上升，而保持椭球形的气泡 沿螺旋形上升；Tripathi 等[26] 认为路径不稳定性和 气泡不对称性密切相关且气泡上升无法表现出稳 定的终端速度；Zhang 等[27]、Chang 等[28] 认为气泡 的非对称形状变形由不对称尾涡的不均匀压差引 起，反向漩涡之间的强度平衡的打破产生了非对 称漩涡对之间的角速度，强烈扭曲的尾流造成气 泡的不对称及不稳定的上升轨迹。在本文后续模 拟工况中，为避免气泡路径不稳定性、轨迹振荡 对聚并过程的影响，将运动区域限定在气泡上升 25 mm 以内。 

![](images/55207ee28a3d7892b0f347d1dd50d71fb6f3de55abe8e0483f31d295de28991e.jpg)



图 4    直径 2.01 mm 单气泡上升过程参数对比



Fig. 4    Comparison of parameter of single bubble rise process with a diameter of 2.01 mm


## 2 两同轴气泡聚并数值模拟结果分析

根据数值模拟结果，如图 5 所示，液相黏度为 0.050 15 Pa·s，直径为 4 mm 的两同轴气泡以初速 度 0 m/s 自由释放，前导气泡率先变形并加速上 升，跟随气泡形变滞后于前导气泡，之后跟随气泡 受到前导气泡尾流的影响，在尾流区逐渐形变并 快速上升追赶前导气泡，直到两气泡接触，接触之 后的液膜厚度不断减薄直至发生聚并。 

![](images/642588611252807f7006360e62ca1beb62ff1191569e371fd36bc2449b27a292.jpg)



图 5    两同轴气泡聚并过程（ µ=0.050 15 Pa·s，d=4 mm）



Fig. 5    Coalescence process of two coaxial bubbles µ（ =0.050 15 Pa·s，d=4 mm）


## 2.1 接触时刻定义

跟随气泡上界面顶点速度 $V _ { \mathrm { t , h } }$ 、前导气泡下界 面顶点速度 $V _ { \mathrm { l , l } }$ 、相对速度 $V _ { \mathrm { r e l } } \left( \ { V } _ { \mathrm { r e l } } = V _ { \mathrm { t , h } } - V _ { \mathrm { l , l } } \right)$ 以及 气泡界面间液体厚度 s、初始液膜厚度 $h _ { 0 } .$ 、液膜厚 度 h 的定义如图 6 所示。气泡界面间液体厚度是 指气泡靠近过程中前导气泡下界面顶点与跟随气 泡上界面顶点间液体的厚度，初始液膜厚度是指 两气泡接触时刻的液体厚度；液膜厚度是指接触 时刻之后排液过程的液体厚度。 

不同直径及黏度工况下，跟随气泡加速靠近 前导气泡过程中的 $V _ { \mathrm { t , h } } \setminus V _ { \mathrm { l , l } } \setminus V _ { \mathrm { r e l } }$ 、s 随时间的变化 如图 7 所示。不同工况下，上述参数随时间的变 化趋势基本一致，即 s 随时间减小， $V _ { \mathrm { l , l } }$ 随时间增 大， $V _ { \mathrm { t , h } }$ 随时间先增大后趋于平缓后下降， $V _ { \mathrm { r e l } }$ 随时 间先增大后减小。 

根据数值模拟结果，前导气泡一直加速直到 聚并，跟随气泡在尾流中加速到一个最大值后减 速，两气泡界面顶点相对速度均呈先增大后减小 的趋势。对于连续相密度远大于离散相的系统， 离散相的相互作用受到连续相的影响，气泡-气 泡、液滴-液滴、刚球-刚球的碰撞存在差异。两气 泡界面顶点间的液体厚度由厚变薄的过程中存在 一临界值 $h _ { 0 }$ 使液体厚度对排液过程影响最为明 显，此时认为排液开始，如图 7 所示，以两气泡界 面顶点相对速度转折点对应的时刻为接触时刻， 所对应的液膜厚度即为初始液膜厚度。 

![](images/41b26471fc96240f36a1ab31432a4788bf3f46e554a9a243b354864746415f0e.jpg)



图 6    同轴气泡聚并过程参数示意图



Fig. 6    Schematic diagram of process parameter of coaxial bubble coalescence


## 2.2 初始液膜厚度

由图 7 所示，相同黏度液体中不同直径气泡 接触时的初始液膜厚度不同，相同直径气泡在不 同黏度液体中的初始液膜厚度也不同。因此，初 始液膜厚度与气泡直径、黏度等系统参数有关， 本文数值模拟结果验证了 Das 等[15] 提出的这一观 点，而并非如 Coulaloglou 和 Tavlarides[3]、Kirkpatrick 和 Lockett[13] 等学者提出的初始液膜厚度 为一定值。 

根据本文数值模拟结果，对气泡直径、黏度对 初始液膜厚度的影响进行了定量分析，结果如图 8 所示。由图 8 可知，当黏度大于 $0 . 0 2 0 0 6 \mathrm { P a \cdot s }$ 时， 一定黏度下的初始液膜厚度随直径的增大而增 大，一定直径下的初始液膜厚度随黏度的增大而 减小。对一定直径的气泡，当黏度低于0.025 075 Pa·s 时，初始液膜厚度受黏度影响最为明显，随着黏度 的增大，初始液膜厚度变化逐渐变缓。 

图 9 为不同初始直径的气泡在接触时刻跟随 气泡长度投影 L 随黏度的变化，可看出，跟随气泡 长度投影随黏度增大而减小。图 10 示出了黏度 为 0.025 075 Pa·s 和 0.075 225 Pa·s 下，初始直径为 5 mm 的跟随气泡形变比随时间的变化，可看出， 黏度越大，跟随气泡的形变比越大，气泡越接近球 形。这与图 8 能够在趋势上很好地对应，即黏度 越大，初始液膜厚度越小。由此考虑到初始液膜 厚度与跟随气泡长度投影之间的联系。 

![](images/90de469f7aba4adbc1d86325507b06423bf3c83e38195a65a7a32469e80b3833.jpg)


![](images/89719888811b0815b9da89663a3e3e51706b7ea529ea931f701cd2c94c8dece5.jpg)


![](images/f7a15b0b9a503e92c278d2b2658afe0046cbfc50a0b1b8f858f65a109e6011eb.jpg)


![](images/6ca3d2e9f9f79014a2eeef8c2b78956e4f7adca90231cbbec1b6d15c8e97d0f6.jpg)



图 7    接触时刻定义



Fig. 7    Definition of contact moment


![](images/0775bed4ddae2de22005dce6a1cfa8f35643bc2aa5890b353e7023395b664d30.jpg)



图 8    初始液膜厚度随黏度的变化



Fig. 8    Initial film thickness vs. viscosity


进一步考察气泡形变比对初始液膜厚度的影 响，图 11 示出了初始液膜厚度与接触时刻跟随气 泡长度投影之比随黏度的变化，可看出，在黏度大 于 0.025 075 Pa·s 时，初始液膜厚度与跟随气泡长 度投影之比趋于稳定在 0.16～0.19 范围内，即对 特定系统而言，当黏度大于某值（在当前系统下认 为是 0.025 075 Pa·s）时，初始液膜厚度与跟随气泡 长度投影 L 有关。本文气泡浮升高度在 25 mm 之 内，且无传热传质，因此忽略气泡体积（A）的变 化。在较高黏度下，气泡可看作椭圆形（x 轴半轴 长 W/2，y 轴半轴长 L/2），有： 

![](images/0f74fae6ac4aa63db25469846da56ede6469c38744496e05883af62e7a233562.jpg)



图 9    跟随气泡长度投影随黏度的变化



Fig. 9    Projection of trailing bubble in length direction vs. viscosity


$$
A = \frac {\pi d ^ {2}}{4} = \frac {\pi L W}{4} \tag {8}
$$

将 $E { = } W / L$ 代入式（8）可得： 

$$
L = d E ^ {- 1 / 2} \tag {9}
$$

若较高黏度下气泡形变可预测，那么初始液 膜厚度就是和初始气泡直径直接相关的函数： $h _ { 0 } =$ $C d E ^ { - 1 / 2 }$ 。C为常数，在较高黏度 $( \mu \geqslant 0 . 0 2 5 0 7 5 \mathrm { P a } \cdot \mathrm { s } )$ 的工况下可认为 C 介于 0.16 与 0.19 之间，以本文 的研究结果，可取 C=0.18。 

![](images/b7a1b021a328b76afc3881a05a2d6ef7e25a7cb2b6180af8c5fe98d4a034047f.jpg)



图 10    不同黏度下跟随气泡形变比（d=5 mm） 随时间的变化



Fig. 10    Bubble deformation ratio (d=5 mm) under different viscosities vs. time


对较低黏度 $( \mu { < } 0 . 0 2 5 \ 0 7 5 \ \mathrm { P a \cdot s } )$ 系统下接触 时刻的跟随气泡形变比预测，以及初始液膜厚度 随黏度、气泡直径的变化规律有待于后续进一步 研究。 

## 2.3 排液时间

目前欧拉-欧拉框架下及欧拉拉格朗日框架 下气液两相流模拟气泡聚并模型多基于可变形、 完全移动界面的排液时间模型，如 Chester 模型[7]、 Lee 模型[8]、Prince 和 Blanch 模型[29]，以上模型均 是在平面膜假设下通过理论推导获得，其具体形 式及适用范围列于表3。 


表 3 排液时间模型



Table 3 Drainage time model


<table><tr><td>作者</td><td>模型</td><td>说明</td></tr><tr><td rowspan="2">Chesters[7]</td><td><eq>t_{\text{d}} = \frac{3\mu_{\text{c}}r}{2\sigma} \ln \frac{h_{0}}{h_{\text{f}}}</eq></td><td>适用于纯黏性力控制下的高黏系统</td></tr><tr><td><eq>t_{\text{d}} = 0.5 \frac{\rho_{\text{c}} U_{\text{rel,0}} d^{2}}{\sigma}</eq></td><td>适用于纯惯性力控制下的低黏系统</td></tr><tr><td>Lee 等[8]</td><td><eq>t_{\text{d}} = \frac{R_{\text{a}}}{4} \left( \frac{\rho_{\text{c}} d}{2\sigma} \right)^{1/2} \ln \frac{h_{0}}{h_{\text{f}}}</eq></td><td>适用于惯性力主导的黏度低于 0.01 Pa·s 的系统</td></tr><tr><td>Prince 和 Blanch[29]</td><td><eq>t_{\text{d}} = \left( \frac{r_{\text{eq}}^{3} \rho_{\text{c}}}{16\sigma} \right)^{1/2} \ln \frac{h_{0}}{h_{\text{f}}}</eq></td><td>适用于低黏系统,如空气-水系统</td></tr></table>

需要说明的是，Chesters 模型[7] 对于其适用范 围只给出定性描述，未给出准确的定量范围，其适 用于低黏、高黏系统下的模型分别是在惯性、黏 性无限大的假设下得出的；Lee 模型[8] 以及 Prince 和 Blanch 模型[29] 的排液时间未包含黏度、界面相 对速度的影响，因而无法体现出惯性力以及黏性 力对排液时间的影响，且这两个预测模型之间的 差异在于初始薄膜半径 ${ r } _ { \mathrm { f } , 0 }$ 的差异，因而预测结果 会相差较小。 

![](images/f35b415c60a458bb2f50eb476d106308c695feedfb807ec94ba2413f53ed9012.jpg)



图 11    初始液膜厚度与跟随气泡长度投影之比 随黏度的变化



Fig. 11    Ratio of initial film thickness to projection of trailing bubble in length direction vs. viscosity


将黏度为 0.010 03 Pa·s 下数值模拟获得的 $R _ { \mathrm { a } } ,$ $r _ { \mathrm { e q } } , U _ { \mathrm { r e l } , 0 \setminus } h _ { 0 }$ 代入 Chesters 模型[7]、Lee 模型[8]、 Prince 和 Blanch 模型[29]，将计算结果与数值模拟 结果进行比较，如表 4 所列， Chester 模型[7] 的预 测结果远小于 Lee 模型[8] 以及 Prince 和 Blanch 模 型[29] 预测的排液时间，且 Prince 和 Blanch 模型[29] 预测的排液时间最大，数值模拟计算结果与模型 预测结果存在较大差异。由此可见，仅考虑惯性 力对气泡排液的影响是不合理的。 


表 4 黏度 0.010 03 Pa·s 的排液时间比较



Table 4 Comparison of drainage time at viscosity of 0.010 03 Pa·s


<table><tr><td rowspan="2">d/mm</td><td colspan="4">排液时间/ms</td></tr><tr><td>Chester 模型</td><td>Lee 模型</td><td>Prince 和 Blanch 模型</td><td>模拟</td></tr><tr><td>4</td><td>5.379</td><td>34.837</td><td>39.970</td><td>28</td></tr><tr><td>5</td><td>16.498</td><td>53.057</td><td>62.126</td><td>28</td></tr><tr><td>6</td><td>21.836</td><td>68.281</td><td>79.544</td><td>24</td></tr><tr><td>7</td><td>31.915</td><td>75.842</td><td>83.181</td><td>27</td></tr><tr><td>8</td><td>27.289</td><td>89.540</td><td>95.416</td><td>37</td></tr></table>

而且，在较低黏度工况下，VOF 方法模拟结果 与 Chester 模型[7]、Lee 模型[8]、Prince 和 Blanch 模 型[29] 的预测结果趋势不同，这需在后续研究中综 合考虑较低黏度下不同直径气泡的形状差异、形变 比变化、气泡上升速度以及黏度对气泡聚并的影响。 

在较高黏度下，将 Chester 模型[7] 与数值模拟 

结果相比较，结果列于表 5。Chester 模型[5] 是在 黏性无限大的极限条件下推导得出的，预测模型 与黏度呈线性正相关，而数值模拟结果显示，较高 黏度下的排液时间与黏度正相关，但并非线性关 系，且黏度越大，预测结果与模拟结果越接近。在 较高黏度下，Chester 模型[7] 预测结果与数值模拟 结果趋势相同，即排液时间随黏度、直径的增大 而增大。 



ZHANG  W,  HU  L,  LI  H,  et  al. Numerical  analysis  of［5］ bubble size effect in a gas-liquid two-phase rotodynamic pump  by  using  a  bubble  coalescence  and  collapse model[J]. Chemical  Engineering  Research  and  Design, 2023, 191: 617-629. 





DAS  S  K. Development  of  a  coalescence  model  due  to［6］ turbulence for the population balance equation[J]. Chemical Engineering Science, 2015, 137: 22-30. 




表 5 较高黏度排液时间比较




CHESTERS A K. The applicability of dynamic-similari-［7］ ty  criteria  to  isothermal,  liquid-gas,  two-phase  flows without  mass  transfer[J]. International Journal  of   Multiphase Flow, 1975, 2(2): 191-212. 




Table 5 Comparison of drainage time at higher viscosity




LEE  C  H,  ERICKSON  L  E,  GLASGOW  L  A. Bubble［8］ breakup  and  coalescence  in  turbulent  gas-liquid  dispersions[J]. Chemical  Engineering  Communications, 1987, 59(1/2/3/4/5/6): 65-84. 



<table><tr><td rowspan="3">d/mm</td><td colspan="6">排液时间/ms</td></tr><tr><td colspan="2">0.025 075 Pa·s</td><td colspan="2">0.050 15 Pa·s</td><td colspan="2">0.075 225 Pa·s</td></tr><tr><td>Chester 模型</td><td>模拟</td><td>Chester 模型</td><td>模拟</td><td>Chester 模型</td><td>模拟</td></tr><tr><td>4</td><td>10.241</td><td>23</td><td>20.369</td><td>26</td><td>30.298</td><td>30</td></tr><tr><td>5</td><td>13.254</td><td>32</td><td>26.071</td><td>33</td><td>38.910</td><td>35</td></tr><tr><td>6</td><td>15.975</td><td>37</td><td>31.751</td><td>40</td><td>47.352</td><td>42</td></tr><tr><td>7</td><td>18.931</td><td>40</td><td>37.328</td><td>46</td><td>55.551</td><td>47</td></tr><tr><td>8</td><td>21.747</td><td>43</td><td>42.794</td><td>50</td><td>63.640</td><td>51</td></tr></table>



CHESTERS A A. The modelling of coalescence process-［9］ es in fluid-liquid dispersions: A review of current understanding[J]. Chemical Engineering Research and Design, 1991, 69(4): 259-270. 



综上所述，如 Liao 和 Lucas[30] 所指出的，在气 液两相具有完全自由移动界面的排液状态的纯净 系统中，如沸腾、气泡塔、油水分离器等系统，气 泡的排液过程是由惯性力和黏性力共同控制的。 对于排液时间的准确预估需综合考虑黏性力和惯 性力的影响。 



李昕晨. 双气泡聚并的流体力学行为研究[D]. 北京: 北［10］ 京化工大学, 2015. 



## 3 结论



ORVALHO S, RUZICKA M C, OLIVIERI G, et al. Bub-［11］ ble  coalescence:  Effect  of  bubble  approach  velocity  and liquid viscosity[J]. Chemical Engineering Science, 2015, 134: 205-216. 



本文通过 VOF 方法结合 CSF 模型对同轴双 气泡聚并过程进行了数值模拟，基于数值模拟结 果获得气泡上升速度场、气泡形变比 E、跟随气 泡长度投影 L、两气泡界面顶点间液体厚度 s 等 微观信息，揭示了同轴双气泡聚并的排液过程中， 不同直径及黏度工况下跟随气泡速度 $V _ { \mathrm { t , h } }$ 、前导气 泡速度 $V _ { \mathrm { l , l } }$ 、气泡相对速度 $V _ { \mathrm { r e l } }$ 、s 随时间变化的特 性和规律，主要结论如下。 



NGUYEN C T, GONNERMANN H M, CHEN Y, et al.［12］ Film  drainage  and  the  lifetime  of  bubbles[J]. Geochemistry, Geophysics, Geosystems, 2013, 14(9): 3616-3631. 





KIRKPATRICK R D, LOCKETT M J. The influence of［13］ approach  velocity  on  bubble  coalescence[J]. Chemical Engineering Science, 1974, 29(12): 2363-2373. 



1） 根据气泡聚并过程相界面变化过程微观模 拟，提出接触时刻及初始液膜厚度的准确定义，即 两气泡液膜相对速度转折点对应的时刻为接触时 刻，所对应的液膜厚度即为初始液膜厚度。 



SANADA  T,  WATANABE  M,  FUKANO  T. Effects  of［14］ viscosity on coalescence of a bubble upon impact with a free  surface[J]. Chemical  Engineering  Science, 2005, 60(19): 5372-5384. 



2） 初始液膜厚度受到气泡直径、液相黏度的 影响，随液相黏度的增大而增大。当黏度大于 0.025 075 Pa·s 时，初始液膜厚度随气泡直径的增 大而增大。 



DAS P K, KUMAR R, RAMKRISHNA D. Coalescence［15］ of drops in stirred dispersion: A white noise model for coalescence[J]. Chemical Engineering Science, 1987, 42(2): 213-220. 



3） $h _ { 0 }$ 与L有关，即与跟随气泡形变比E有关。 当黏度大于 0.025 075 Pa·s 时， $h _ { 0 } { = } C d E ^ { - 1 / 2 }$ ，C 为一 



龚升高. 湍流条件下液滴或气泡破裂和聚并过程研［16］ 究[D]. 湘潭: 湘潭大学, 2017. 



介于 0.16～0.19 之间的常数，在较高黏度下可认 为 C=0.18。 



闻昭权. 气泡聚并机理的实验及模型研究[D]. 湘潭: 湘［17］ 潭大学, 2019. 





SONG  R,  HAN  L,  ZHANG  L,  et  al. Experiments  and［18］ modeling  of  bubbles  colliding  head-on  in  water[J]. 



4） 黏度为 0.010 03 Pa·s 时，Lee 模型[8] 以及 Prince 和 Blanch 模型[29] 无法体现惯性力和黏性力 对排液时间的影响，对排液时间的预测均大于 Chester 模型[7]；较高黏度下的排液时间与黏度正 相关，但并非 Chester 模型[7] 中的线性关系；排液 时间随黏度、直径的增大而增大。 

5） 排液过程是由惯性力和黏性力共同控制 的，预测排液时间需综合考虑二者的影响。 



AIChE Journal, 2021, 67(6): e17220. 



## 参考文献：



GANESAN P B, ISLAM M T, POKRAJAC D, et al. Co-［19］ alescence and rising behavior of co-axial and lateral bubbles in viscous fluid: A CFD study[J]. Asia-Pacific Journal of Chemical Engineering, 2017, 12(4): 605-619. 





ZHANG Y, CHEN K, YOU Y, et al. Coalescence of two［20］ initially  spherical  bubbles:  Dual  effect  of  liquid viscosity[J]. International Journal of Heat and Fluid Flow, 2018, 72: 61-72. 





OWOEYE  E  J,  SCHUBRING  D. Numerical  simulation［1］ of vapor bubble condensation in turbulent subcooled flow boiling[J]. Nuclear  Engineering  and  Design, 2015, 289: 126-143. 





MERABTENE T, GAROOSI F, MAHDI T F. Numerical［21］ modeling of liquid spills from the damaged container and collision of  two  rising  bubbles  in  partially  filled   enclosure  using  modified  Volume-Of-Fluid  (VOF)  method[J]. Engineering  Analysis  with  Boundary  Elements, 2023, 154: 83-121. 





COULIBALY A, BI J, LIN X, et al. Effect of bubble coa-［2］ lescence on the wall heat transfer during subcooled pool boiling[J]. International  Journal  of  Thermal  Sciences, 2014, 76: 101-109. 





ABBASSI W, BESBES S, ELHAJEM M, et al. Numeri-［22］ cal simulation of free ascension and coaxial coalescence of  air  bubbles  using  the  volume  of  fluid  method (VOF)[J]. Computers & Fluids, 2018, 161: 47-59. 





COULALOGLOU  C  A,  TAVLARIDES  L  L. Descrip-［3］ tion of interaction processes in agitated liquid-liquid dispersions[J]. Chemical Engineering Science, 1977, 32(11): 1289-1297. 





BRACKBILL J U, KOTHE D B, ZEMACH C. A contin-［23］ uum  method  for  modeling  surface  tension[J]. Journal  of Computational Physics, 1992, 100(2): 335-354. 





ZHANG  X  B,  LUO  Z  H. Effects  of  bubble  coalescence［4］ and  breakup  models  on  the  simulation  of  bubble columns[J]. Chemical  Engineering  Science, 2020, 226: 115850. 





LIU  L,  YAN  H,  ZHAO  G. Experimental  studies  on  the［24］ shape and motion of air bubbles in viscous liquids[J]. Experimental Thermal and Fluid Science, 2015, 62: 109-121. 





WU M, GHARIB M. Experimental studies on the shape［25］ and  path  of  small  air  bubbles  rising  in  clean  water[J]. Physics of Fluids, 2002, 14(7): L49-L52. 





TRIPATHI M K, SAHU K C, GOVINDARAJAN R. Dy-［26］ namics of an initially spherical bubble rising in quiescent liquid[J]. Nature Communications, 2015, 6: 6268. 





ZHANG  J,  NI  M  J. What happens  to  the  vortex   struc-［27］ tures when  the  rising  bubble  transits  from  zigzag  to   spiral?[J]. Journal of Fluid Mechanics, 2017, 828: 353-373. 





CHANG  Y,  MÜLLER  C,  KOVÁTS  P,  et  al. Hydrody-［28］ namics and shape reconstruction of single rising air bubbles  in  water  using  high-speed  tomographic  particle tracking velocimetry and 3D geometric reconstruction[J]. Experiments in Fluids, 2023, 65(1): 6. 





PRINCE  M  J,  BLANCH  H  W. Bubble  coalescence  and［29］ break-up in air-sparged bubble columns[J]. AIChE Journal, 1990, 36(10): 1485-1499. 





LIAO Y, LUCAS D. A literature review on mechanisms［30］ and  models  for  the  coalescence  process  of  fluid particles[J]. Chemical Engineering Science, 2010, 65(10): 2851-2864. 

