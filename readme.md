\documentclass[10pt]{article}
\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage{tensor}
\usepackage{ctex}
\usepackage{cite}
\usepackage[a4paper, margin=1in]{geometry}
\usepackage{xcolor}
\setlength{\abovedisplayskip}{5pt}
\setlength{\belowdisplayskip}{5pt}
\title{Dy+ASMDO+PPC+FTSMC}
\author{Yu Dong}
\date{\today}

\begin{document}
\maketitle

\section{前言}


\section{模型}
\subsection{机械臂模型}
\begin{equation}
    M_m(q) \ddot{q} + C_m(q, \dot{q}) \dot{q} + G_m(q) = \tau_{q} + \tau_{lum}
\end{equation}

\subsection{无人机模型}
牛顿欧拉法：
\begin{align}
     & \dot{P}_{b}={v}_{b}                                                                                                                                \\
     & \dot{v}_{b}=-\frac{F_{l}}{m_{s}} I R_{B} e_{3}+ ge_{3}+\frac{F_{cd}}{m_{s}}+d_{lumA}                                                               \\
     & \dot{\Phi_{b}}=T(\Phi_{b})\tensor[^B]{\omega}{_{b}}                                                                                                \\
     & \dot{\tensor[^B]{\omega}{_{b}}}=I^{-1}_{b}(\tau-\tensor[^B]{\omega}{_{b}}\times(I_{b}\tensor[^B]{\omega}{_{b}})+\tensor[^B]{\tau}{_{cd}})+d_{lumB}
\end{align}

\textcolor{red}{\textbf{控制方程1:姿态动力学变量为角速度}}
\[
    \begin{cases}
        \dot{x}_1 = x_2                                                                                                                                       \\
        \dot{x}_2 = \frac{1}{m_s} \left( u_x + F_{cdx} \right) + d_{lum_1}                                                                                    \\
        \dot{x}_3 = x_4                                                                                                                                       \\
        \dot{x}_4 = \frac{1}{m_s} \left( u_y + F_{cdy} \right) + d_{lum_2}                                                                                    \\
        \dot{x}_5 = x_6                                                                                                                                       \\
        \dot{x}_6 = \frac{1}{m_s} \left( u_z + F_{cdz} \right) + g + d_{lum_3}                                                                                \\
        \dot{x}_7 = x_8 + f_{\phi}(x_a)                                                                                                                       \\
        \dot{x}_8 = \frac{1}{J_{\phi}} \left( \tau_{\phi} + \tensor[^{B}]{\tau}{_{cd_{\phi}}} + (J_{\psi} - J_{\theta}) x_{10} x_{12} \right) + d_{lum_4}     \\
        \dot{x}_9 = x_{10} + f_{\theta}(x_a)                                                                                                                  \\
        \dot{x}_{10} = \frac{1}{J_{\theta}} \left( \tau_{\theta} + \tensor[^{B}]{\tau}{_{cd_{\theta}}} + (J_{\phi} - J_{\psi}) x_8 x_{12} \right) + d_{lum_5} \\
        \dot{x}_{11} = x_{12} + f_{\psi}(x_a)                                                                                                                 \\
        \dot{x}_{12} = \frac{1}{J_{\psi}} \left( \tau_{\psi} + \tensor[^{B}]{\tau}{_{cd_{\psi}}} + (J_{\theta} - J_{\phi}) x_8 x_{10} \right) + d_{lum_6}
    \end{cases}\\
\]

\textcolor{red}{\textbf{控制方程1说明：}}

1.基本定义

$\tensor[^I]{R}{_{B}}$:联系无人机位置$P_{b}$和速度$\dot{v_{b}}$，为坐标系$\sum B$对坐标系$\sum I$的旋转矩阵
\begin{align}
    \tensor[^I]{R}{_{B}} & =R_{z}(\psi)R_{y}(\theta)R_{x}(\psi)                                                                             \\
                         & =\begin{bmatrix}
                                c\phi & -s\phi & 0 \\
                                s\phi & c\phi  & 0 \\
                                0     & 0      & 1 \\
                            \end{bmatrix}
    \begin{bmatrix}
        c\theta  & 0 & s\theta \\
        0        & 1 & 0       \\
        -s\theta & 0 & c\theta \\
    \end{bmatrix}
    \begin{bmatrix}
        1 & -0    & 0      \\
        0 & c\psi & -s\psi \\
        0 & s\psi & c\psi  \\
    \end{bmatrix}                                                                                                                      \\
                         & =\begin{bmatrix}
                                cos\phi cos\theta & -sin\psi cos\psi + cos\psi sin\theta sin\phi & sin\phi sin\psi + cos\phi sin\theta cos\psi  \\
                                sin\phi cos\theta & cos\phi cos\psi +sin\phi sin\theta sin\psi   & -cos\phi sin\psi + sin\phi sin\theta cos\psi \\
                                -sin\theta        & cos\theta sin\psi                            & cos\theta cos\psi                            \\
                            \end{bmatrix}
\end{align}
$T(\Phi_{b})$:联系无人机欧拉角$\Psi$和机体角速度$\omega$
\begin{align}
    \tensor[^I]{\omega}{_{B}} & =T(\Phi_{b})\dot{\Phi_{b}}                 \\
    T(\Phi_{b})_{rot-zyx}     & =\begin{bmatrix}
                                     1 & sin\phi tan\theta & cos\phi tan\theta \\
                                     0 & cos\phi           & -sin\theta        \\
                                     0 & sin\phi sec\theta & cos\psi sec\theta \\
                                 \end{bmatrix}
\end{align}
2.状态变量表达式
\begin{align}
    x_1 = x \quad       & \quad x_2 = \dot{x}           \\
    x_3 = y \quad       & \quad x_4 = \dot{x}_3         \\
    x_5 = z \quad       & \quad x_6 = \dot{x}_5         \\
    x_7 = \psi \quad    & \quad x_8 = \omega_{x} = p    \\
    x_9 = \theta \quad  & \quad x_{10} = \omega_{y} = q \\
    x_{11} = \phi \quad & \quad x_{12} = \omega_{z} = r
\end{align}
\[
    \begin{cases}
        u_x=-F_l(cos\phi sin\theta cos\psi + sin\phi sin\psi) \\
        u_y=-F_l(cos\phi sin\theta sin\psi - sin\phi cos\psi) \\
        u_z=-F_l cos\phi cos\theta                            \\
    \end{cases}
\]
\[
    \begin{cases}
        f_{\phi}(x_a) = q \sin \phi \tan \theta + r \cos \phi \tan \theta \\
        f_{\theta}(x_a) = q \cos \phi - r \sin \phi - q                   \\
        f_{\psi}(x_a) = q \sin \phi \sec \theta + r \cos \phi \sec \theta - r
    \end{cases}
\]

$F_{cd_{x}}, F_{cd_{y}}, F_{cd_{z}}, \tensor[^{B}]{\tau}{_{cd_{\theta}}},\tensor[^{B}]{\tau}{_{cd_{\psi}}},\tensor[^{B}]{\tau}{_{cd_{\phi}}}$：机械臂对无人机的扰动(直接计算部分)

$d_{lum_{1-6}}$：不确定性、外部扰动、机械臂对无人机不可直接计算部分扰动

3.耦合扰动
\begin{equation}
    \left\{
    \begin{aligned}
        \tensor[^{B}]{\mathbf{F}}{_{cd}} & = -m_{man} \, \tensor[^{I}]{R}{_{B}} \left( \tensor[^{B}]{\omega}{_{b}} \times \left( \tensor[^{B}]{\omega}{_{b}} \times \tensor[^{B}]{r}{_{omc}} \right) + \tensor[^{B}]{\dot{\omega}}{_{b}} \times \tensor[^{B}]{r}{_{omc}} + 2 \, \tensor[^{B}]{\omega}{_{b}} \times \tensor[^{B}]{\dot{r}}{_{omc}} + \tensor[^{B}]{\ddot{r}}{_{omc}} \right)                                                        \\
        \tensor[^{B}]{\hat{\tau}}{_{cd}} & = m_s \left( \tensor[^{B}]{r}{_{oc}} \times \left( \tensor[^{B}]{R}{_{I}} g \mathbf{e_3} \right) - \tensor[^{B}]{r}{_{oc}} \times \tensor[^{B}]{\ddot{r}}{_{o}} - \tensor[^{B}]{\dot{r}}{_{oc}} \times \tensor[^{B}]{\dot{r}}{_{o}} \right) - m_{man} \left( \tensor[^{B}]{\dot{r}}{_{o}} \times \tensor[^{B}]{\dot{r}}{_{omc}} + \tensor[^{B}]{r}{_{o}} \times \tensor[^{B}]{\ddot{r}}{_{omc}} \right) \\
                                         & + \tensor[^{B}]{\omega}{_{b}} \times \left( \tensor[^{B}]{r}{_{omc}} \times \tensor[^{B}]{\dot{r}}{_{omc}} \right) + \tensor[^{B}]{r}{_{omc}} \times \tensor[^{B}]{\ddot{r}}{_{omc}} - \tensor[^{B}]{I}{_{man}} \tensor[^{B}]{\omega}{_{b}} \times \tensor[^{B}]{\omega}{_{b}} - \tensor[^{B}]{I}{_{man}} \tensor[^{B}]{\dot{\omega}}{_{b}}
    \end{aligned}
    \right.
\end{equation}

$r_{omc}$为机械臂质心,$r_{oc}$为系统质心,机械臂运动导致质心及惯性矩变动

有文献认为$\dot{\omega}_{b},\ddot{r}_{oc}$不可测,但$\dot{r}_{oc}$可测，不可测采用卡尔曼滤波近似；

另有认为$\dot{\omega}_{b},\ddot{r}_{oc},\dot{r}_{oc}$皆不可测，全部放入ESO扰动观测器中；

还有认为$\dot{\omega}_{b},\ddot{r}_{oc}$不可测，但$\dot{r}_{oc}$可测，不可测部分采用滤波近似（纯前馈）；

UAM整体质心：
\begin{equation}
    \left\{
    \begin{aligned}
        \tensor[^{B}]{\mathbf{r}}{_{oc}}       & = \frac{1}{m_s} \sum_{i=1}^{n} m_i \tensor[^{B}]{\mathbf{p}}{_{ci}}, \\
        \tensor[^{B}]{\dot{\mathbf{r}}}{_{oc}} & = \frac{1}{m_s} \sum_{i=1}^{n} m_i \tensor[^{B}]{\mathbf{v}}{_{ci}},
    \end{aligned}
    \right.
\end{equation}

机械臂连杆位置及速度：
\begin{equation}
    \left\{
    \begin{aligned}
         & {}^B \mathbf{p}_{ci} = {}^B \mathbf{T}_i (\mathbf{q}) {}^i \mathbf{r}_{ci}, \\
         & \begin{bmatrix}
               {}^B \mathbf{v}_{ci} \\
               {}^B \boldsymbol{\omega}_i
           \end{bmatrix} = {}^B \mathbf{J}_{ci} (\mathbf{q}) \dot{\mathbf{q}},
    \end{aligned}
    \right.
\end{equation}

机械臂质心：
\begin{equation}
    \tensor[^{B}]{\mathbf{r}}{_{omc}} = \frac{m_s}{m_{man}} \tensor[^{B}]{\mathbf{r}}{_{oc}},
\end{equation}

% Second Set of Equations
机械臂惯性矩：基于平行移轴定理
\begin{equation}
    \left\{
    \begin{aligned}
        \tensor[^{B}]{\mathbf{I}}{_{man}^{o}} & = \sum_{i=1}^{n} \left( \tensor[^{B}]{\mathbf{R}}{_{i}} \tensor[^{ci}]{\mathbf{I}}{_{i}^{B}} \tensor[^{B}]{\mathbf{R}}{_{i}^{-1}} + m_i ( \|\tensor[^{B}]{\mathbf{p}}{_{ci}}\|^2 \mathbf{I}_{3\times3} - \tensor[^{B}]{\mathbf{p}}{_{cn}} (\tensor[^{B}]{\mathbf{p}}{_{cn}})^{T}) \right),                                                     \\
        \tensor[^{B}]{\mathbf{I}}{_{man}^{o}} & = \sum_{i=1}^{n} m_i \left( 2 (\tensor[^{B}]{\mathbf{p}}{_{ci}})^{T} \tensor[^{B}]{\mathbf{v}}{_{ci}} \mathbf{I}_{3\times3} - \tensor[^{B}]{\mathbf{v}}{_{ci}} (\tensor[^{B}]{\mathbf{p}}{_{ci}})^{T} - \tensor[^{B}]{\mathbf{p}}{_{ci}} (\tensor[^{B}]{\mathbf{v}}{_{ci}})^{T} \right)                                                        \\
                                              & \quad + \sum_{i=1}^{n} \left( \text{Skew}(\tensor[^{B}]{\omega}{_{i}}) \tensor[^{B}]{\mathbf{R}}{_{i}} \tensor[^{ci}]{\mathbf{I}}{_{i}^{B}} \tensor[^{B}]{\mathbf{R}}{_{i}^{-1}} - \tensor[^{B}]{\mathbf{R}}{_{i}} \tensor[^{ci}]{\mathbf{I}}{_{i}^{B}} \tensor[^{B}]{\mathbf{R}}{_{i}^{-1}} \text{Skew}(\tensor[^{B}]{\omega}{_{i}}) \right),
    \end{aligned}
    \right.
\end{equation}

\textcolor{red}{\textbf{控制方程2:姿态动力学变量为欧拉角}}

无人机模型：
\[
    \left\{
    \begin{aligned}
        \dot{p}_b        & = v_b                                                      \\
        m_u \dot{v}_b    & = f R_b e_3 - m_u g e_3 + F_{cd} + d_{lumA}                \\
        \dot{\Phi}_b     & = T(\Phi_b) \omega_b                                       \\
        J \dot{\omega}_b & = \tau - \omega_b \times (J \omega_b) + \tau_{cd}+d_{lumB}
    \end{aligned}
    \right.
\]

整体模型：
\begin{align}
    \begin{cases}
        \dot{x}_1 = x_2 \\
        H(x_1) \dot{x}_2 + C(x_1,x_2) x_2 + G(x_1) = u + d_{cp} +d_{lum}
    \end{cases}
\end{align}
\[x_1 = \begin{bmatrix} p^{T}_{b} & \Phi^{T}_{b} & q^{T} \end{bmatrix}
\]

整体模型简化：
\begin{align}
    \begin{cases}
        \dot{x}_1 = x_2                                         \\
        \dot{x}_2 = f(\vec{x}) + h(\vec{x})(u + d_{cp}+d_{lum}) \\
        y = x_1
    \end{cases}
\end{align}

\textcolor{red}{\textbf{控制矩阵2说明：}}

1.无人机姿态动力学化简
姿态动力学方程利用欧拉角化简:
\begin{align}
     & I_{b} T^{-1}(\Phi_b) \ddot{\Phi}_b  - I_{b} T^{-1}(\Phi_b) \dot{T} (\Phi_b) T^{-1} (\Phi_b) \dot{\Phi}_b \notag + (T^{-1} (\Phi_b) \dot{\Phi}_b) \times I_{b} (T^{-1} (\Phi_b) \dot{\Phi}_b) \notag \\
     & =  \tau + \tensor[^B]{\tau}{_{cd}} + d_{lum}
\end{align}

质量阵对角化：两边同乘$(T^{-1}(\Phi_b))^{T}$
\begin{align}
     & (T^{-1}(\Phi_b))^{T} I_{b} T^{-1}(\Phi_b) \ddot{\Phi}_b  - (T^{-1}(\Phi_b))^{T} I_{b} T^{-1}(\Phi_b) \dot{T} (\Phi_b) T^{-1} (\Phi_b) \dot{\Phi}_b \notag \\
     & + (T^{-1}(\Phi_b))^{T} (T^{-1} (\Phi_b) \dot{\Phi}_b) \times I_{b} (T^{-1} (\Phi_b) \dot{\Phi}_b) \notag                                                  \\
     & =  (T^{-1}(\Phi_b))^{T} \tau + (T^{-1}(\Phi_b))^{T} \tensor[^B]{\tau}{_{cd}} + (T^{-1}(\Phi_b))^{T} d_{lum}
\end{align}

2.惯性矩阵和扰动
\[ H(x_1)=
    \begin{bmatrix}
        m_u I_{3 \times 3} &                           &     \\
                           & T^{-1}(\Phi)^{T} J T^{-1} &     \\
                           &                           & M_m \\
    \end{bmatrix}
\]

\begin{align*}
    d_{cp} & = [F_{cd_{3 \times 1}},(T^{-1}(\Phi_b))^{T} \tau_{cd_{3 \times 1}},0_{m \times 1}] \\
    d_2    & = [d_{lum_{1-6}},\tau_{lum}]
\end{align*}

3.广义控制输入$u_{1-(6+m)}$

无人机控制分为四个层次：位置控制、姿态控制、控制分配、电机控制。

多旋翼系统有6个自由度，4个独立控制输入（升力$F_{l}$和三个力矩$M_{xyz}$），只能跟踪4个期望指令（偏航角和三个位置），剩余的期望指令（预期俯仰角和预期横滚角）由前者计算得到；

$u_{1},u_{2}$（位置控制器）根据期望偏航计算期望滚转、俯仰;然后通过$u_{2-6}$完成控制力$F_{l}$（位置控制器）/力矩$\tau$（姿态控制器）;最后力/力矩通过动力分配矩阵转为电机期望转速（电机控制器）

$u_1,u_2$获取期望$\psi_d,\theta_d$ :
\begin{align*}
    u_{1} & = F_{l}(cos \psi sin \theta cos \phi + sin \psi sin \phi) \\
    u_{2} & = F_{l}(cos \psi sin \theta sin \phi - sin \psi cos \phi) \\
    \begin{bmatrix}
        u_{1} \\
        u_{2}
    \end{bmatrix}
          & \approx F_{l}
    \begin{bmatrix}
        \theta_{d} cos \phi_{d} + \psi_{d} sin \phi_{d} \\
        \theta_{d} sin \phi_{d} + \psi_{d} cos \phi_{d} \\
    \end{bmatrix}
    =F_{l}
    \begin{bmatrix}
        sin \phi_{d} cos \phi_{d} \\
        -cos \psi_{d} sin \psi_{d}
    \end{bmatrix}
    \begin{bmatrix}
        \psi_{d} \\
        \theta_{d}
    \end{bmatrix}                                                    \\
    F_{l} &
    \begin{bmatrix}
        sin \phi_{d} cos \phi_{d} \\
        -cos \psi_{d} sin \psi_{d}
    \end{bmatrix}
    \begin{bmatrix}
        \psi_{d} \\
        \theta_{d}
    \end{bmatrix}
    =
    \begin{bmatrix}
        u_{1} \\
        u_{2} \\
    \end{bmatrix}                                                    \\
\end{align*}


综合广义控制输入：
\begin{align*}
    \begin{bmatrix}
        \psi_{d} \\
        \theta_{d}
    \end{bmatrix}
     & = \frac{1}{F_{l}}
    \begin{bmatrix}
        sin \psi_{d} & -cos \psi_{d} \\
        cos \psi_{d} & sin \psi_{d}
    \end{bmatrix}
    \begin{bmatrix}
        u_{1} \\
        u_{2}
    \end{bmatrix}       \\
    \begin{bmatrix}
        f \\
        M \\
        \tau
    \end{bmatrix}
     & =
    \begin{bmatrix}
        R_b(3,3)           & 0                  & 0                  \\
        0_{3 \times 3}     & (T^{-1}(\Phi_b))^T & 0_{3 \times 3}     \\
        0_{n_r \times n_r} & 0_{n_r \times n_r} & I_{n_r \times n_r}
    \end{bmatrix}^{-1}
    \begin{bmatrix}
        u_3    \\
        \vdots \\
        u_{6+n_r}
    \end{bmatrix}
\end{align*}

四旋翼控制分配矩阵：
\[
    \begin{bmatrix}
        f   \\
        M_1 \\
        M_2 \\
        M_3
    \end{bmatrix}
    =
    \begin{bmatrix}
        c_f   & c_f    & c_f    & c_f   \\
        d c_f & 0      & -d c_f & 0     \\
        0     & -d c_f & 0      & d c_f \\
        -c_M  & c_M    & -c_M   & c_M
    \end{bmatrix}
    \begin{bmatrix}
        \omega_1^2 \\
        \omega_2^2 \\
        \omega_3^2 \\
        \omega_4^2
    \end{bmatrix}
\]

4.简化模型

状态向量:(第$i$个自由度的状态$x_1,x_2$)
\[
    \vec{x}=f(x_{1,i},x_{2,i}),i=6+m
\]

重力项和科氏力：
\[f(\vec{x}) = -H^{-1}(x_1)\left(C(x_1, x_2)x_2 + G(x_1)\right)\]
\[h(\vec{x}) = H^{-1}(x_1)\]


\section{自适应固定时间滑模扰动观测器}

优势：固定时间收敛（与初值无关）、收敛速度快、精确度高

数学角度理解：待补充

\subsection*{3.1 扰动观测器滑模面}
1.\textbf{滑模面:}
系统动力学方程：
\begin{align}
    H(x_1) \dot{x}_2 + C(x_1, x_2) x_2 + G(x_1) = U + d_{cp} + d_{um}
\end{align}

辅助动力学方程：
\begin{align}
    H(x) \dot{\eta} + C x_2 + G & = u + d_{cp} + \hat{d}_{um} + v + \frac{1}{N(s)} \left( \lambda_1 \text{sig}^\alpha (s) + \lambda_2 \text{sig}^\beta \left( s \right) \right) \\
    v                           & = k sign(s)
\end{align}

引入滑模变量：
\begin{align}
    s = x_2 - \eta
\end{align}

\begin{align}
    \rightarrow H \dot{S} = \hat{d}_{um} - \frac{1}{N(s)} \left( \lambda_1 \text{sig}^{1 + \sigma} (s) + \lambda_2 \text{sig}^{\frac{\alpha}{\beta}} \left( s \right) \right) - v
\end{align}

系数：
\begin{align}
    k = \bar{d}_{lum} + \| \hat{d}_{lum} \| + \epsilon
    ,\epsilon > 0 \\
    \sigma = \left( \frac{m_1}{2 n_1}\right)\left(1+sgn\left(\| s \| -1\right)\right)
    ,\sigma = \begin{cases}
                  \frac{m}{n} > 1, & \| S \|  > 1 \\
                  0,               & \| S \| < 1
              \end{cases}
    \\
    m > n >0 ,\beta >\alpha>0, m\text{、}n\text{、}\beta\text{、}\alpha= 2K+1
\end{align}

2.\textbf{滑模面固定时间收敛及证明:}
\begin{align}
    \dot{V}_s & = 2S H \dot{S}                                                                                                                                                                                                                                                    \\
              & = 2S \left[ \tilde{d}_{lum} - \frac{1}{N(s)} \left( \lambda_1 \text{sig}^{1+\delta} (s) + \lambda_2 \text{sig}^{\frac{\alpha}{\beta}} \left( s \right) \right) - v \right]                                                                                        \\
              & = 2S \left[ d_{lum} - \hat{d}_{lum} - \frac{1}{N(s)} \left( \lambda_1 \text{sig}^{1+\delta} (s) + \lambda_2 \text{sig}^{\frac{\alpha}{\beta}} \left( \dot{s} \right) \right) - \left( \bar{d}_{um} + \| \hat{d}_{um} \| + \epsilon \right) \text{sign}(s) \right] \\
              & \leq -\frac{2\lambda_1}{N(s)} \|S\|^{2+\delta}- \frac{2\lambda_2}{N(s)} \|S\|^{1+\frac{\alpha}{\delta}}                                                                                                                                                           \\
              & \leq -\frac{2\lambda_1}{N(s)\lambda_m^{\frac{2+\delta}{m}} } \|V_s\|^{\frac{2+\delta}{2}} - \frac{2\lambda_2}{N(s)\lambda_m^{\frac{\alpha + \beta}{2 \beta}}}  \|V_s\|^{\frac{\alpha + 1}{2\beta}}                                                                \\
              & \rightarrow Fix \text{-} time
\end{align}

\begin{flalign*}
                & \textbf{\textcolor{red} {Proof 1:}  \text{滑模面收敛}}                                                                                                                                         & \\
                & \text{design \space} \dot{y} = \left(1 - \frac{\alpha}{2 \beta}\right) V_{s}^{-\frac{\alpha +1 }{2 \beta}} \dot{V}_{s} \space, \text{\space}  y = V_{s}^{1-\frac{\alpha+1}{2 \beta}}          & 
            \end{flalign*}
            \[
            \begin{aligned}
                &\dot{V}_{s}  = \frac{\dot{y}}{(1-\frac{\alpha + 1 }{2 \beta}) V^{-\frac{\alpha + 1}{2 \beta}}}                                                                                                                                                                                                                                     \\
                & \leq -\frac{2 \lambda_{1}}{N\left(s\right) \lambda_{m}^{\frac{2 + \delta}{2}}} \| V_{s} \| ^{\frac{2+ \delta }{2}} - \frac{2 \lambda_2}{N\left(s\right) \lambda_{m}^{\frac{\alpha+\beta}{2 \beta}} } \| V_{s} \| ^{\frac{\alpha+1}{2 \beta}}                                                                                                            \\
                & \rightarrow \text{\space} \dot{y} \leq \left(1-\frac{\alpha + 1 }{2 \beta}\right) \left[-\frac{2 \lambda_{1}}{N\left(s\right) \lambda_{m}^{\frac{2 + \delta}{2}}} \| V_{s} \| ^{1 + \frac{\delta }{2} - \frac{\alpha + 1}{2 \beta}} - \frac{2 \lambda_2}{N\left(s\right) \lambda_{m}^{\frac{\alpha+\beta}{2 \beta}} } \| V_{s} \| ^{0}  \right]         \\
                & = \left( 1-\frac{\alpha + 1 }{2 \beta}  \right) \left(-\frac{2 \lambda_{1}}{N\left(s\right) \lambda_{m}^{\frac{2 + \delta}{2}}} y^{1+ \frac{\delta}{2}\left( \frac{2 \beta}{2\beta - \alpha -1}\right)} - \frac{2 \lambda_2}{N\left(s\right) \lambda_{m}^{\frac{\alpha+\beta}{2 \beta}} }      \right)                                                  \\
                & = \frac{r}{N\left(s\right)} \left(-\acute{\lambda}_{1} y^{1+ \frac{\delta}{2}\left( \frac{2 \beta}{2\beta - \alpha -1}\right)} - \acute{\lambda}_{2} \right)                                                                                                                                                                                            \\
                & \rightarrow \text{\space} dy  \leq \frac{r}{N\left(s\right)} \left(-\acute{\lambda}_{1} y^{1+ \frac{\delta}{2}\left( \frac{2 \beta}{2\beta - \alpha -1}\right)} - \acute{\lambda}_{2} \right)  dt                                                                                                                                                       \\
                & dt \leq -\frac{N\left(s\right)}{r}  \frac{1}{-\acute{\lambda}_{1} y^{1+ \frac{\delta}{2}\left( \frac{2 \beta}{2\beta - \alpha -1}\right)} - \acute{\lambda}_{2}} dy                                                                                                                         \\
                & \rightarrow \int_{0}^{T_{s}}dt \leq  \int_{0}^{y\left(0\right)} \frac{N\left(s\right)}{r}  \frac{1}{\acute{\lambda}_{1} y^{1+ \frac{\delta}{2}\left( \frac{2 \beta}{2\beta - \alpha -1}\right)} + \acute{\lambda}_{2}} dy                                                                                                                               \\
                \\
                T_{s}& \leq \int_{0}^{1}\frac{N\left(s\right)}{r}  \frac{1}{\acute{\lambda}_{1} y^{1+ \frac{\delta}{2}\left( \frac{2 \beta}{2\beta - \alpha -1}\right)} + \acute{\lambda}_{2}} dy + \int_{1}^{T_{s}}\frac{N\left(s\right)}{r}  \frac{1}{\acute{\lambda}_{1} y^{1+ \frac{\delta}{2}\left( \frac{2 \beta}{2\beta - \alpha -1}\right)} + \acute{\lambda}_{2}} dy  \\
                & \leq \int_{0}^{1}\frac{N\left(s\right)}{r}  \frac{1}{\acute{\lambda}_{1} y + \acute{\lambda_{2}}} dy + \int_{1}^{T_{s}}\frac{N\left(s\right)}{r}  \frac{1}{\acute{\lambda}_{1} y^{1+ \frac{m}{n}\left( \frac{ \beta}{2\beta - \alpha -1}\right)} + \acute{\lambda}_{2}} dy                                                                              \\
                & \leq \int_{0}^{1}\frac{1}{r}  \frac{1}{\acute{\lambda}_{1} y + \acute{\lambda_{2}}} dy + \int_{1}^{T_{s}}\frac{1}{r}  \frac{1}{\acute{\lambda}_{1} y^{1+ \frac{m}{n}\left( \frac{ \beta}{2\beta - \alpha -1}\right)} + \acute{\lambda}_{2}} dy                                                                                                          \\
                & = \frac{1}{r} \int_{0}^{1} \left(    \frac{1}{\acute{\lambda}_{1} y + \acute{\lambda_{2}}} dy + \int_{1}^{T_{s}}  \frac{1}{\acute{\lambda}_{1} y^{\rho} + \acute{\lambda}_{2}} dy  \right)                                                                                                                                                              \\
                & \leq \frac{1}{r} \int_{0}^{1} \left(    \frac{1}{\acute{\lambda}_{1} y+\acute{\lambda}_{2} } dy + \int_{1}^{T_{s}}  \frac{1}{\acute{\lambda}_{1} y^{\rho}} dy  \right)                                                                                                                                                                                  \\
                & =\frac{1}{r} \left(  \frac{1}{\acute{\lambda}_{1} \ln ( 1+ \frac{\acute{\lambda}_1}{\acute{\lambda}_2})} + \frac{\acute{\lambda_1}}{1-\rho}\left(y(0)^{-\rho +1}-1\right)\right)                                                                                                                                                                        \\
                & \leq \frac{1}{r} \left(\frac{1}{\acute{\lambda}_{1} \ln ( 1+ \frac{\acute{\lambda}_1}{\acute{\lambda}_2})} + \frac{\acute{\lambda_1}}{\rho -1 } \right)                                                                                                                                                                                                 \\
                & T_{s}\text{收敛与初值无关 ,滑模面收敛,\space Proof end}  & \\
\end{aligned}
\]
系数：
\begin{align}
    r                   & = \left( 1-\frac{\alpha + 1 }{2 \beta}  \right)                                                                                                                               \\
    \acute{\lambda_{1}} & = \frac{2 \lambda_{1}}{N\left(s\right) \lambda_{m}^{\frac{2 + \delta}{2}}} \text{,\space}\acute{\lambda{2}} = \frac{2 \lambda_2}{\lambda_{m}^{\frac{\alpha+\beta}{2 \beta}} } \\
    \rho                & = 1+ \frac{m}{n}\left( \frac{ \beta}{2\beta - \alpha -1}\right)>1
\end{align}

备注：
\begin{align}
    \begin{cases}
        if \text{\space} y & \in \left(0 ,1 \right) \rightarrow V_{s} \in \left(0,1\right) \rightarrow \delta = 0 \\
        if \text{\space} y & \in \left(0 ,1 \right) \rightarrow V_{s} \in \left(0,1\right) \rightarrow \delta = 0 \\
    \end{cases}
\end{align}


\subsection*{3.2 自适应固定时间滑模扰动估计器}
\begin{flalign}
    \hspace{2em}\text{design ASMDO:}                &                                                                                                                                                                                                                                     & \\
                                                    & \begin{cases}
                                                          \hat{d}_{\text{lum}} = \zeta + k_{2} H x_{2} &                                                                                                                                                      \\
                                                          \dot{\zeta} = \frac{1}{N\left(v\right)}\left(\lambda_{3} \text{sig}^{1+\delta}(v) + \lambda_{4} \text{sig}^{\frac{\alpha}{\beta}}(v)\right) + k_{2} \left(C x_{2} + G - \hat{d}_{\text{lum}}\right) \\
                                                      \end{cases} &      \\
    \rightarrow \text{\space} \dot{\tilde{d}}_{lum} & = -\frac{1}{N}\left(\lambda_{3} \text{sig}^{1+\delta}(v) + \lambda_{4} \text{sig}^{\frac{\alpha}{\beta}}(v)\right)  - k_2 \tilde{d}_{lum} + \dot{d}_{lum}                                                                           &
\end{flalign}
\begin{flalign}
    \hspace{2em}\text{ Lemma: } &                                                                                                                      & \\
                                & \text{if \space} \dot{V} \leq -a v^{r_{1}} - bV^{r_2} = \rho ,\text{\space} r_1 > 1,0<r_2<1,\rho > 0                 & \\
                                & \text{system convergence :\space} T_{max} \leq \frac{1}{a\left(r_1 - 1 \right)} + \frac{1 }{b\left(1 - r_{2}\right)} & 
\end{flalign}

\textbf{\textcolor{red} {Proof 2:}  \text{扰动估计误差收敛}}

\hspace{2em}\text{滑模面收敛：\space} $\left( \tilde{d}_{lum}\right)_{eq} = v = k sign(s) $
\[
\begin{aligned}
    \hspace{6em} V & = \tilde{d}_{lum}^{2}                                                                                                                                                                             & \\
    \dot{V}        & = 2 \tilde{d}_{lum} \dot{\tilde{d}}_{lum}                                                                                                                                                         & \\
                   & = -\frac{2}{N} \tilde{d}_{lum}  \left[ \left(\lambda_{3} \text{sig}^{1+\delta}(v) + \lambda_{4} \text{sig}^{\frac{\alpha}{\beta}}(v)\right) +\dot{d}_{lum}\right]                                   \\
                   & = -\frac{2}{N} \lambda_{3} \tilde{d}_{lum}^{2+\delta} -\frac{2}{N} \lambda_{4} \tilde{d}_{lum}^{\frac{\alpha}{\beta}+1} - 2 k_2 \tilde{d}_{lum}^{2} + 2 \tilde{d}_{lum} \dot{d}_{lum}             & \\
                   & \leq -\frac{2}{N} \lambda_{3} \tilde{d}_{lum}^{2+\delta}  -\frac{2}{N} \lambda_{4} \tilde{d}_{lum}^{\frac{\alpha}{\beta}+1} - 2 k_2 \tilde{d}_{lum}^{2} + \title{d}_{lum}^{2} + \dot{d}_{lum}^{2} & \\
                   & =-\frac{2}{N} \lambda_{3} \tilde{d}_{lum}^{2+\delta}  -\frac{2}{N} \lambda_{4} \tilde{d}_{lum}^{\frac{\alpha}{\beta}+1} +(1- 2 k_2) \tilde{d}_{lum}^{2}  + \dot{d}_{lum}^{2}                      & \\
                   & \leq -\frac{2}{N} \lambda_{3} V^{\frac{2+\delta}{2}} - -\frac{2}{N} \lambda_{4} V^{\frac{\alpha+ \beta}{2 \beta}} + c_{d}                                                                         &
\end{aligned}
\]
\hspace{6em}Consider Lemma : system convergence fix time , Proof end \\

Poof 2 说明：
\begin{align}
     & \begin{cases}
            & \frac{\alpha+\beta}{2 \beta} < 1                                                                        \\
            & \text{if \space} |\tilde{d}|>1 \rightarrow \|v \|>1 \text{\space,\space} \frac{2+\delta}{2}>1           \\
            & \text{if \space} |\tilde{d}|\leq 1 \rightarrow \|v \| \leq1 \text{\space,\space} \frac{2+\delta}{2} = 1
       \end{cases} \\
     & \text{$\| \tilde{d}\|$超出1就会迅速收敛至可控范围}
\end{align}

系数：
\begin{align}
    c_{d} = \dot{d}_{lum}^{2}
\end{align}



\section{不受初始误差影响的固定时间PPC}

\subsection*{4.1 基本概念}
预设性能控制：Prescribed Performance Control

核心思想：对误差人为设定性能包络，通过性能函数的收敛特性刻画被控系统的瞬态和稳态；

关键点：性能函数；转换函数；控制器

经典PPC:

\hspace*{2em}(1)性能函数
\begin{align}
    \rho(t) = \left(\rho_{0} - \rho_{\infty}\right)e^{-lt} + \rho_{\infty}
\end{align}

\hspace*{2em}性能函数约束
\begin{align}
                & -\rho(t) < e(t) < \rho(t)                    \\
    \text{精细化：} & \begin{cases}
                      -\delta \rho(t) < e(t) < \rho(t), & e(0) > 0 \\
                      \rho(t) < e(t) < \delta \rho(t),  & e(0) < 0 \\
                  \end{cases}
\end{align}

\hspace*{2em}(2)转换方程(PPF)

\hspace*{2em}\text{目的：将受约束方程利用空间对等映射转为无约束系统}

\hspace*{2em}\text{同胚映射后新系统状态$\varsigma$有界即可实现对原系统约束控制}
\begin{align}
    e(t)                      & = \Xi(\varsigma)\rho(t) \text{,\space} \Xi(\varsigma) = \frac{e(t)}{\rho(t)}                                                                                                                              \\
    \Xi(\varsigma)\text{选择要求} & \begin{cases}
                                    \text{光滑，连续，单调递减}                                                                                                                                                                                         \\
                                    \text{if \ } e(0)> 0 \text{,\ } \Xi(\varsigma)\in (-\delta,1)\text{,\ } \lim\limits_{\varsigma \to -\infty}{\Xi(\varsigma)} = -\delta\text{,\ }  \lim\limits_{\varsigma \to \infty}{\Xi(\varsigma)} = 1   \\
                                    \text{if \ } e(0) < 0 \text{,\ } \Xi(\varsigma)\in (-1,\delta)\text{,\ } \lim\limits_{\varsigma \to -\infty}{\Xi(\varsigma)} = -1\text{,\ }  \lim\limits_{\varsigma \to \infty}{\Xi(\varsigma)} = -\delta \\
                                \end{cases}
\end{align}

\hspace*{2em}\text{常见转换函数$\Xi(\varsigma)$}
\begin{align}
     & \Xi(\varsigma) =
    \begin{cases}
        \frac{e^{\varsigma} - \delta e ^{- \varsigma}}{e^{\varsigma} +  e ^{- \varsigma}}\text{\space if \space} e(0) > 0     \\
        \frac{\delta e^{\varsigma} -  e ^{- \varsigma}}{e^{\varsigma} +  e ^{- \varsigma}}\text{\space if \space} e(0) \leq 0 \\
    \end{cases} \\
     & \varsigma = \Xi^{-1}\frac{e(t)}{\rho(t)}  =  \begin{cases}
                                                        \frac{1}{2} \ln(\frac{e(t)/\rho(t)+\delta}{1-e(t)/\rho(t)})\text{,\space if \space}e(0)>0   \\
                                                        \frac{1}{2} \ln(\frac{e(t)/\rho(t)+1}{\delta - e(t)/\rho(t)})\text{,\space if \space}e(0)<0 \\
                                                    \end{cases}
\end{align}

\hspace*{2em}(3)转换方程
待补充

%4.2 New PPFunction
\subsection*{4.2 New PPF}
优势：无需初始误差，固定时间快速收敛

数学角度：待补充 %补充第一个函数调整瞬态超调，第二个函数调整收敛时间


\[
\begin{aligned}
    &(1)\quad\rho(t) = \begin{cases}
        coth(kt + \iota)e^{\frac{1}{t-T}}+t_{\infty}, & t\leq T \\
        t_{\infty},                                   & t > T
    \end{cases}\\
    &\hspace{2em}\text{$\iota$ 为设定的无穷小参数，$k$为设定的参数，$T$为设定的收敛时间，$t_{\infty}$ 为设定收敛值}\\
    &(2)\quad\rho(t) = \begin{cases}
        \frac{\left(T-t\right)^{a}}{t^b+\iota}+t_{\infty}, & t\leq T \\
        t_{\infty},                                        & t > T
    \end{cases}\\
    &\hspace{2em}\text{$\iota$ 为设定的无穷小参数，$a,b$为设定的幂次，$T$为设定的收敛时间，$t_{\infty}$ 为设定收敛值}\\
    &(3)\quad\rho(t)= \begin{cases}
        \frac{1}{\left(e^{t+\iota}-1\right)e^{\frac{1}{T-t}}}+t_{\infty}, & t\leq T \\
        t_{\infty},   
        &\text{$\iota$ 为设定的无穷小参数，$T$为设定的收敛时间，$t_{\infty}$ 为设定收敛值}\\ & t > T
    \end{cases}\\
    &(4)\quad\rho(t) = \begin{cases}
        \frac{a}{t+\iota},                                                & 0\leq t <T_1     \\
        e^{-\frac{t}{T_2-t}}\left(\rho_0 - \rho_\infty\right)+t_{\infty}, & T_1 \leq t < T_2 \\
        \rho_{\infty},                                                    & t > T_2
    \end{cases}\\
    &\hspace{2em}\text{$\iota$ 为设定的无穷小参数，$a,\rho_0$ 为计算的参数}\\
    &\hspace{2em}\text{$T_1$为设定的超调时间，$T_2$为设定的收敛时间，$t_{\infty}$ 为设定收敛值}\\
    &\text{$a,\rho_0$计算推导:}    \\
    & \text{函数$t = T_1$时连续：} \frac{a}{T_1+\iota} = e^{-\frac{T_1}{T_2-T_1}}\left(\rho_0 - \rho_\infty\right)+T_{\infty}   \\
    & \text{导数$t = T_1$时连续：} -a  \frac{1}{(T_1 + \iota)^{2}} = -  \frac{T_2 -T_1 + T_1}{(T_2 - T_1)^2} e^{-\frac{T_1}{T_2 -  T_1}} (\rho_0 - \rho_\infty)            \\
    \end{aligned}\\
    \]
\[
    \begin{aligned}
         & T_1 + \iota = \frac{\frac{T_2}{(T_2-T_1)^{2}} e^{- \frac{T_1}{T_2-T_1}} (\rho_0 - \rho_{\infty}) }{e^{-\frac{T_1}{T_2-T_1}}(\rho_0 -\rho_{\infty}) + \rho_{\infty}}               \\
         & (T_1 + \iota)(e^{-\frac{T_1}{T_2-T_1}}(\rho_0 -\rho_{\infty}) + \rho_{\infty}) = \frac{T_2}{(T_2-T_1)^{2}} e^{- \frac{T_1}{T_2-T_1}} (\rho_0 - \rho_{\infty})                     \\
         & (T_1 + \iota)((\rho_0 -\rho_{\infty}) + e^{\frac{T_1}{T_2-T_1}} \rho_{\infty}) = \frac{T_2}{(T_2-T_1)^{2}}(\rho_0 - \rho_{\infty})                                                \\
         & (T_1 + \iota)((\rho_0 -\rho_{\infty}))+(T_1 + \iota)e^{\frac{T_1}{T_2-T_1}} \rho_{\infty} = \frac{T_2}{(T_2-T_1)^{2}}(\rho_0 - \rho_{\infty})                                     \\
         & (\rho_0 -\rho_{\infty})(\frac{T_2}{(T_2-T_1)^{2}} - T_1 -\iota) = (T_1 + \iota) e^{\frac{T_1}{T_2-T_1}} \rho_{\infty}                                                             \\
         & \rho_0 = \frac{(T_1 + \iota) e^{\frac{T_1}{T_2-T_1}} \rho_{\infty}}{\frac{T_2}{(T_2-T_1)^{2}} - T_1 -\iota} + \rho_{\infty}                                                       \\
         & a = (T_1+\iota)(e^{-\frac{T_1}{T_2-T_1}}\left(\rho_0 - \rho_\infty\right)+T_{\infty})                                                                                             \\
         & a =(T_1+\iota)(e^{-\frac{T_1}{T_2-T_1}}\left(\frac{(T_1 + \iota) e^{\frac{T_1}{T_2-T_1}} \rho_{\infty}}{\frac{T_2}{(T_2-T_1)^{2}} - T_1 -\iota}- \rho_{\infty}\right)+T_{\infty}) \\
    \end{aligned}
\]



\end{document}
