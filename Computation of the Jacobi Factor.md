[[2026-03-09]]
可执行的被积函数的变量代换：
$$
\int\frac{b^{2}}{r^{2}\sqrt{ 1-u^{2} }}\,dr=\int\frac{-b^{2}}{\sqrt{ 1-u^{2} }}\,d\left( \frac{1}{r} \right)
$$

---
[[2026-03-12]]
1. 角动量守恒式需要引入无穷远处等离子体折射率修正
2. 分母根式下的变量代换不是单调的，因此需要分两部分积分，接下来我们需要讨论RN时空的被积函数形状
	1. 是否可以只看第一部分的贡献

绝对值符号的处理：
```mathematica
Simplify[Abs[r], Assumptions -> r > 0]
```

Schiwartzchild 偏折角
$$
\frac{b}{r^2 \sqrt{\frac{b^2 (2 M-r)+r^2 \left(2 M \epsilon ^2-r \epsilon ^2+r\right)}{r^3}}}
$$
RN偏折角
$$
\frac{b \sqrt{1-\epsilon ^2}}{r \sqrt{r^2-\frac{\left(r^2 \epsilon ^2-b^2 \left(\epsilon ^2-1\right)\right) \left(r (r-2 M)+Q^2\right)}{r^2}}}
$$

---
我们仍希望构造分母为$r^{2}\sqrt{ 1-u^{2} }$的形式，因此先抓出根号下的部分，可以构造关于$r$和`u2`$=u^{2}$的方程 

$$
r^2-\frac{\left(r^2 \epsilon ^2-b^2 \left(\epsilon ^2-1\right)\right) \left(r (r-2 M)+Q^2\right)}{r^2}=r^{2}(1-u^{2})
$$
变量代换$x=\frac{1}{r}$，整理一下，原方程化为
$$
\left(b^2 x^2 \left(\epsilon ^2-1\right)-\epsilon ^2\right) \left(2 M x-Q^2 x^2-1\right)=u^{2}
$$
求解$x(u^{2})$，即用$u^{2}$表示雅可比因子$J(u^{2})=\frac{dx}{du^{2}}$中所有的$x$，即可将RN时空的偏折角积分化为我们希望的如下形式
$$
\int_{u_{1}^{2}}^{u_{0}^{2}}+\int_{u_{0}^{2}}^{u_{2}^{2}} b\sqrt{ 1-\epsilon^{2} }\cdot\frac{J(u^{2})}{\sqrt{  1-u^{2}}} \,du^{2}
$$
```mathematica
Simplify[Expand[Denominator[d\[Phi]odr]/r^2]];
Simplify[1 - %^2 /. r -> 1/x] == u2 (*u2是u^2的整体代换*)
xtou2 = Solve[%, x][[3, 1, 2]]
```

先讨论积分限$u_{1,0,2}$，$u_{0}$ 是$u(r)$的最小值，且$r\neq \infty$。$u_{1}=u(r_{0})=1$,也定义了$r_{0}$和$b$的关系。 $u_{2}=u(\infty)=\epsilon$

---
testparams的选取需要求解以$r_{0}$为未知数的近日点方程，即
$$
u^{2}(r_{0})=1
$$
解得：$r_{0}=r_{0}(b)$
```mathematica
testparamsRN = {Q -> 2/5, M -> 1, r0 -> 10000};
SimPla = {\[Omega]e[r] -> 
   Sqrt[\[Epsilon]^2 \[Omega]\[Infinity]^2 r^-k], \[Omega]e[r0] -> 
   Sqrt[\[Epsilon]^2 \[Omega]\[Infinity]^2 r0^-k]};
kreplace = {k -> 0};
ToAssumptions = {\[HBar] > 0, \[Omega]\[Infinity] > 0, r > 2  M, 
  r0 > 2  M, ri > r0, rf > r0, M > 0, r0 > 0, b > 0, 
  s[0] > 0, \[Omega]\[Infinity] > 0, 1 > u > 0, \[Epsilon] > 0, v > 0,
   v < 1, n\[Infinity] > 0};
dphiodr = 
 Simplify[-(dphiodk/drodk) /. 
      spr /. {pt -> \[HBar]  \[Omega]\[Infinity], 
      pphi -> (\[HBar] \[Omega]\[Infinity] b)/n\[Infinity]}, 
    Assumptions -> {\[HBar] > 0, \[Omega]\[Infinity] > 0}]
   RNMetric = {Am[r] -> 1 - (2 M)/r + Q^2/r^2, 
   Bm[r] -> 1/(1 - (2 M)/r + Q^2/r^2), Cm[r] -> r^2, 
   Am[r0] -> 1 - (2 M)/r0 + Q^2/r0^2, 
   Bm[r0] -> 1/(1 - (2 M)/r0 + Q^2/r0^2), Cm[r0] -> r0^2};
Simplify[
 d\[Phi]odr /. RNMetric /. SimPla /. 
   kreplace /. {n\[Infinity] -> 1/Sqrt[1 - \[Epsilon]^2]}, 
 Assumptions -> ToAssumptions];
Simplify[Denominator[%] /. {\[Omega]\[Infinity] -> 1}];
denodphiodrsq = %^2;
Solve[% == 0 /. {r -> r0}, b];
%[[2, 1, 2]];
Numerator[%]/r0;
-%^2;
% /. Take[testparamsRN, 3];
Solve[(99999/10)^2 == %, \[Epsilon]];
testparamsRN = Join[testparamsRN, %[[2]]];
Solve[denodphiodrsq == 0 /. {r -> r0} /. testparamsRN, b];
testparamsRN = Join[testparamsRN, %[[2]]];
{\[Epsilon]2 -> \[Epsilon]^2 /. testparamsRN};
(testparamsRN = Join[testparamsRN, %]
```

---
Jacobi factor and the minimum of u2[x]
```mathematica
Simplify[Expand[Denominator[d\[Phi]odr]/r^2]];
du2 = D[1 - %^2 //. {r -> 1/x}, x] dx;
J[x] = dx/du2 // Simplify // Expand;
J[u2] = J[x] /. x -> xtou2;
Print["the Jacobian factor J[(\!\(\*SuperscriptBox[\(u\), \(2\)]\))]=\
\!\(\*FractionBox[\(dx\), SuperscriptBox[\(du\), \(2\)]]\) is saved \
in function J[u2]"];
```

```mathematica
xMin = Solve[du2/dx == 0 , x][[2, 1, 2]];
Simplify[Expand[Denominator[d\[Phi]odr]/r^2]];
Print["minimum of (\!\(\*SuperscriptBox[\(u\), \(2\)]\))=" N[
    u2Min = 1 - %^2 //. {r -> 1/x} /. x -> xMin /. testparamsRN, 
    60]];
(*\[Epsilon] is appoximated to 0.00002000390074814339`*)
```

$$
\frac{1}{J(u^{2})}=0\implies x=x_{min}\,\,\,\,u^{2}(x_{min})=u_{0}^{2}
$$

Do Taylor expansion on J[u2] in the form of (u2-u2Min)^n
$$
J(u^{2})=\sum J_{i} (u^{2}-u_{0}^{2})^{i}
$$
before expansion, replace all $u^{2}$ in $J(u^{2})$ with $u_{0}^{2}+u^{2}_{\mu}$, where $u^{2}_{\mu}$ is written as `u2mu2smallRN` in Mathematica.
$$
J(u^{2})=J(u_{\mu}^{2})=\sum J_{i}(u_{\mu}^{2})^{i/2}
$$
```mathematica
J[u2] // 
  Expand /. 
   u2 -> u0 + 
     u2mu2smallRN ;(*here we use u0 to replace uMin temporarily*)
Series[J[u2], {u2mu2smallRN, 0, 1}]
```

$u_{\mu}$ 用$\epsilon$ 展开的结果应是只有偶数次幂 the series of u2 in $\epsilon^{2}$ can be solved from the equation
``` mathematica
(-4 (12 M \[Epsilon]^2 (-b^2 M + b^2 M \[Epsilon]^2) + (b^2 - 
       b^2 \[Epsilon]^2 + Q^2 \[Epsilon]^2)^2 + 
     12 (-u2 + \[Epsilon]^2) (b^2 Q^2 - 
        b^2 Q^2 \[Epsilon]^2))^3 + (108 (-u2 + \[Epsilon]^2) (-b^2 M +
        b^2 M \[Epsilon]^2)^2 + 
    36 M \[Epsilon]^2 (-b^2 M + b^2 M \[Epsilon]^2) (b^2 - 
       b^2 \[Epsilon]^2 + Q^2 \[Epsilon]^2) + 
    2 (b^2 - b^2 \[Epsilon]^2 + Q^2 \[Epsilon]^2)^3 + 
    108 M^2 \[Epsilon]^4 (b^2 Q^2 - b^2 Q^2 \[Epsilon]^2) - 
    72 (-u2 + \[Epsilon]^2) (b^2 - b^2 \[Epsilon]^2 + 
       Q^2 \[Epsilon]^2) (b^2 Q^2 - b^2 Q^2 \[Epsilon]^2))^2)==0
```
$$
\left(108 M^2 \epsilon ^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+36 M \epsilon ^2 \left(b^2 M \epsilon ^2-b^2 M\right) \left(b^2 \left(-\epsilon ^2\right)+b^2+Q^2 \epsilon ^2\right)+108 \left(\epsilon ^2-\text{u2}\right) \left(b^2 M \epsilon ^2-b^2 M\right)^2-72 \left(\epsilon ^2-\text{u2}\right) \left(b^2 \left(-\epsilon ^2\right)+b^2+Q^2 \epsilon ^2\right) \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+2 \left(b^2 \left(-\epsilon ^2\right)+b^2+Q^2 \epsilon ^2\right)^3\right)^2-4 \left(12 M \epsilon ^2 \left(b^2 M \epsilon ^2-b^2 M\right)+12 \left(\epsilon ^2-\text{u2}\right) \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+\left(b^2 \left(-\epsilon ^2\right)+b^2+Q^2 \epsilon ^2\right)^2\right)^3=0
$$
用待定系数法，将u2Min写成关于$\epsilon^{2}$的多项式形式
$$
u_{2\text{guess}}=\sum_{i=0}^{5}\text{uu}_{i}\,\epsilon^{2}i=u_{0}^{2}(\epsilon^{2};b,M,Q) 
$$
构造$u_{2\text{guess}}-u_{0}^{2}$, 令$\epsilon^{2i}$各项系数均为零，联立求解$uu_{i}$

```mathematica
u2Guess = Sum[uu[i]*\[Epsilon]2^i, {i, 0, 2}];
u2Min /. \[Epsilon] -> Sqrt[\[Epsilon]2];
Solve[CoefficientList[u2Guess - %, \[Epsilon]2] == 0, {uu[0], uu[1], 
  uu[2]}]
```

### 组会纪要[[2026-03-12]]
1. 决定$u_{0}^{2}$的方程变复杂，需要继续测试展开的正确性
2. 需要考虑度规不再是多项式的情形
	1. 整体代换？
	2. **待定系数法**求解各种展开<=>用代数方程求解再展开，就比如上文所述的$uu_{i}$的求解
	3. 微分方程也可以用待定级数法求解，因为级数和多项式的微分是能求得（高数上讲过的办法）
3. Charged Hayward 
4. 调研：均匀等离子体+RN时空的偏折角文献
5. 用待定系数法求$x(u^{2})$

---
## [[2026-03-17]]
focus on this equation
$$
\left(b^2 x^2 \left(\epsilon ^2-1\right)-\epsilon ^2\right) \left(2 M x-Q^2 x^2-1\right)=u^{2}
$$
等式左边展开，对$x$合并同类项，可以构造关于$x$的四次方程:
$$
x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x+(\epsilon ^2 - u^{2})=0
$$
>[!bug] [[2026-03-26]] 常数项错误订正
> 经后续帕德逼近展开检查发现该方程有误，化简之后常数项应是$+(\epsilon^{2}-u^{2})$而不是$-(\epsilon^{2}+u^{2})$，物理上，当$u^{2}=\epsilon^{2}$时，$r\to \infty, x=0$，方程应化为齐次方程，拥有平凡解。

从物理上我们知道：
$$
x=\frac{1}{r}, r_{0}<r<+\infty\implies x \in \left[0, \frac{1}{r_{0}} \right]
$$
从数学上我们知道:
$$
u^{2}=\left(b^2 x^2 \left(\epsilon ^2-1\right)-\epsilon ^2\right) \left(2 M x-Q^2 x^2-1\right) \implies u^{2} \in [u_{0}^{2},1]
$$
且$u_{0}$略小于$\epsilon$

所以该四次方程$Ax^{4}+Bx^{3}+Cx^{2}+Dx+E=0$的系数分别为：
$$
\begin{equation}\begin{aligned}
&A=b^{2}Q^{2}(1-\epsilon^{2})\sim r^{4}\\
&B=2b^{2}M(1-\epsilon^{2})\sim r^{3}\\
&C=b^{2}\left( 1-\epsilon^{2}+\frac{Q^{2}}{b^{2}}\epsilon^{2} \right)\sim r^{2}\\
&D=-2M\epsilon^{2} \sim r\\
&E=(\epsilon^{2}-u^{2})

\end{aligned}\end{equation}
$$
$A,B,C>0\,\,D<0$

### $Q\to{0}$时的 解结构分析
在上周我们用代数求解（Ferrari方法），得到了$x(u^{2})$的解析形式：
$$
x = -\frac{B}{4A} \pm \frac{1}{2}\left( \pm \sqrt{2z - p} \pm \sqrt{-2z - p \pm \frac{2q}{\sqrt{2z - p}}} \right)
$$
其中
$$
\begin{aligned} p &= \frac{8AC - 3B^2}{8A^2} \\ q &= \frac{B^3 - 4ABC + 8A^2D}{8A^3} \\ r &= \frac{-3B^4 + 256A^3E - 64A^2BD + 16AB^2C}{256A^4} \end{aligned}
$$
而 $z$ 是伴随三次方程的三个可能的根，不是单值的。这个伴随三次方程可以写为以$p,q,r$为参数，以$z$为未知数的形式：
$$
z^3 - \frac{p}{2}z^2 - rz + \left(\frac{rp}{2} - \frac{q^2}{8}\right)=0
$$
其显式解为：
$$
z = \frac{p}{6} + \sqrt[3]{ \frac{p^3}{216} - \frac{pr}{6} + \frac{q^2}{16} + \sqrt{\left(\frac{p^3}{216} - \frac{pr}{6} + \frac{q^2}{16}\right)^2 - \left(\frac{p^2}{36} + \frac{r}{3}\right)^3} }

+ \sqrt[3]{ \frac{p^3}{216} - \frac{pr}{6} + \frac{q^2}{16} - \sqrt{\left(\frac{p^3}{216} - \frac{pr}{6} + \frac{q^2}{16}\right)^2 - \left(\frac{p^2}{36} + \frac{r}{3}\right)^3} }
$$
亦可将其写作三根分立的形式：
$$
z_k = \frac{p}{6} + \omega^k U + \omega^{-k} V,\quad k=0,1,2
$$
where $\omega = e^{2\pi i/3}$, $U,V$是两个立方根。

当$Q\to0$时，$A\to0$。Ferrari显式解显含$-\frac{B}{4A},\quad p,q,r \sim \frac{1}{A^2},\frac{1}{A^3},\frac{1}{A^4}$，会导致发散。如果对四个根都求极限，应当有三个根退化为Cardano解，另一个出现奇异摄动，发散到无穷大。

可以设计一个mathematica实验，来检测取$A\to0$极限时，Ferrari解能否收敛到三次方程（史瓦西时空）的Cardano解，或者说奇异摄动是否会影响后续计算：

![[截屏2026-03-17 下午9.06.44.png]]

### 待定系数法解$x(u^{2})$[[2026-03-18]]
由于Ferrari解形式复杂，后续积分Jacobian factor时麻烦。遂尝试用待定系数法解这个方程。从而给出级数形式的$x(u^{2})$。

>[!quote] 级数形式的$x(u^{2})$高阶项截断，是否会丢失积分主值？
>

且如果度规函数关于r再高一阶：$A(r)=1-\frac{2M}{r}+\frac{Q^{2}}{r^{2}}+\frac{a^{3}}{r^{3}}$, $x$与$u^{2}$的数量关系将变为5次方程。因此必须找到一个半解析的办法，给出$x$关于$u^{2}$的表达式。

>[!quote] 我想给出x关于u^2的半解析表达式，因为四次方程的Ferrari解实在太复杂了，不利于分析。有没有可用的方法，比如待定系数法，将x写成u^2的级数形式，然后求各次项系数的线性方程组？
>令 $\delta = u^2\ (\ge 0)$，方程重写为  $p(x) - \delta = 0$ ,其中  
$$p(x) = A x^4 + B x^3 + C x^2 + D x - \epsilon^2$$  
>（A、B、C、D 完全不含 $u^2$，仅 $\epsilon^2$ 在常数项）。
>
>假设  
$$x(\delta) = x_0 + x_1 \delta + x_2 \delta^2 + x_3 \delta^3 + \cdots$$  
>（$x_0$ 是 $u=0$ 时的正实根Root4，satisfies $p(x_0)=0$）。
>
>把 $x = x_0 + s$（$s = x_1\delta + x_2\delta^2 + \cdots$）代入 $p(x_0 + s) = \delta$， Taylor 展开$$p(x_0) + p'(x_0)s + \frac{p''(x_0)}{2}s^2 + \frac{p'''(x_0)}{6}s^3 + \frac{p^{(4)}(x_0)}{24}s^4 + \cdots = \delta$$
>
>因为 $p(x_0)=0$，且 $s$ 的最低阶是 $\delta^1$，逐阶令 $\delta^k$ 系数为 0，即得线性递推：$$\delta^1:\quad p'(x_0)x_1 - 1 = 0 \quad \Rightarrow \quad x_1 = \frac{1}{p'(x_0)}$$
>$$\delta^2:\quad p'(x_0)x_2 + \frac{p''(x_0)}{2}x_1^2 = 0 
\quad \Rightarrow \quad x_2 = -\frac{p''(x_0)x_1^2}{2p'(x_0)}$$
>$$\delta^3:\quad p'(x_0)x_3 + p''(x_0)x_1 x_2 + \frac{p'''(x_0)}{6}x_1^3 = 0
\quad \Rightarrow \quad x_3 = -\frac{p''(x_0)x_1 x_2 + \frac{p'''(x_0)}{6}x_1^3}{p'(x_0)}$$
>$$\delta^4:\quad p'(x_0)x_4 + p''(x_0)\left(x_1 x_3 + \frac{x_2^2}{2}\right) + \frac{p'''(x_0)}{2}x_1^2 x_2 + \frac{p^{(4)}(x_0)}{24}x_1^4 = 0$$
>以此类推,  $p(x)$ 是四次多项式，所以：
>$$p'(x) = 4A x^3 + 3B x^2 + 2C x + D$$
>$$p''(x) = 12A x^2 + 6B x + 2C$$
>$$p'''(x) = 24A x + 6B$$
>$$p^{(4)}(x) = 24A \quad（更高阶导数为 0）$$

>[!bug] 在u->1时，x1 u2 ，x2 u4，x3 u6 等项发散得非常快，下面是代入具体数值计算出的结果![[Pasted image 20260318210144.png|300]]这意味着解的行为在u=1附近不良好，怎么办
>尝试`AsyptoticSolve`函数

#### [[2026-03-19]] 解的渐进逼近
MMA的AsymptoticSlove函数有一些特殊的要求，我做了以下两组实验测试：
![[截屏2026-03-19 下午4.46.58.png]]
第二条指令没有输出的原因在于最后一项常数项除了因变量u2外还多了$\epsilon^{2}$这一参数，因此在后续计算时需要注意常数项不能有多项。

基于此，定义$y=\epsilon^{2}+u^{2}$, 分别执行$y=1$（近日点）和$y=0$（无穷远）附近，$x$的渐进展开，观察解的行为。
![[截屏2026-03-19 下午5.13.07.png]]
经过尝试，无论怎样调整x1=Root[-1 - DD #1 + CC #1^2 + BB #1^3 + AA #1^4 &, 4], 均无法有效计算渐进展开，尝试失败。

#### 帕德逼近

### [[2026-03-19]]组会纪要
真分式×根式=雅可比因子，这个式子积分是否可积？
$$
\int_{u_{1}^{2}}^{u_{0}^{2}}+\int_{u_{0}^{2}}^{u_{2}^{2}} b\sqrt{ 1-\epsilon^{2} }\cdot\frac{\,J(u^{2})}{\sqrt{  1-u^{2}}} \,du^{2}
$$
待定系数法假设的xApprox形式未必是u2的整数次幂级数，目前史瓦西时空是**半整数阶的幂级数。**
数理方法告诉我们：u2=u2small展开，在u2small到1之间，在复平面内有奇点u2singular，从而在usmall附近展开的级数注定收敛不到u=1的值。
	因此在u2small-u2singular之间，用泰勒展开可以make sense
	在u2singular所在的环与实轴的交点->1之间，用洛朗展开（可行域是个环）才能收敛
![[截屏2026-03-19 下午6.58.29.png|400]]

$$
\sqrt{ p(u_{2}) }=\sum p_{i} u_{2}^{i}
$$
其中$p(u_{2})$是关于$u_{2}$的三次多项式，且当$u_{2}\in[u_{2}small,1]$时$p(u_{2})<0$, 泰勒展开后出现了指数爆炸。
	p是负数时，根式函数是多值函数，我们需要在某个单值分支下才能展开。
	调研：需要在不同的范围内有不同的展开方式，导致级数不能工作的原因也可能有复平面内奇点之外的影响因素。需要找复分析和复变函数的资料以解决根号下负值函数的级数展开有意义的条件。

---
### Pade Approximant 的可行性分析[[2026-03-26]]
在检查出$p(x)\text{:=}\sum _{k=0}^n a(k) x^k$的系数有误之后，重新数值检验了帕德逼近的收敛性。

>[!bug] [[2026-03-26]] 常数项错误订正
> 经后续帕德逼近展开检查发现该方程有误，化简之后常数项应是$+(\epsilon^{2}-u^{2})$而不是$-(\epsilon^{2}+u^{2})$，物理上，当$u^{2}=\epsilon^{2}$时，$r\to \infty, x=0$，方程应化为齐次方程，拥有平凡解。


1. 修改四次方程的系数，尤其需要注意，待定系数法求解$x_{0}$和解析法求Ferrari解$x$的方程有所不同

关于$x_{0}$的方程是
$$
x_{0}^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x_{0}^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x_{0}^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x_{0}+\underbrace{ \epsilon ^2 }_{没有u2 } =0
$$
修改后的定义如下：
```mathematica
testAK = {a[0] -> \[Epsilon]^2, a[1] -> DD, a[2] -> CC, a[3] -> BB, 
   a[4] -> AA}; 
$Assumptions = 
  a[0] > 0 && a[1] < 0 && a[2] > 0 && a[3] > 0 && a[4] > 0; 
n = 4; 
p[x_] := Sum[a[k]*x^k, {k, 0, n}];
x0sols = Solve[p[x] == 0 && x > 0, x, Reals];
```

关于$x$的方程是
$$
x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x+\underbrace{ (\epsilon ^2 - u^{2}) }_{ 有u2 }=0
$$
修改后的定义如下:
```mathematica
EE + DD  x + CC x^2 + BB x^3 + AA x^4 == 0;
abcde = {AA -> b^2 Q^2 (1 - \[Epsilon]^2), 
   BB -> 2  b^2  M  (\[Epsilon]^2 - 1), 
   CC -> b^2 (1 - \[Epsilon]^2 + (Q/b)^2 \[Epsilon]^2), 
   DD -> -2  M  \[Epsilon]^2, EE -> \[Epsilon]^2 - u2};
testparamsRN = {Q -> 2/5, M -> 1, 
   r0 -> 10000, \[Epsilon] -> Sqrt[199999/624875001]/4, 
   b -> 999990000/Sqrt[9997800017], \[Epsilon]2 -> 199999/
    9998000016};
FerrariRoots = 
  Solve[AA x^4 + BB x^3 + CC x^2 - DD x - EE == 0 && x > 0, x, 
   Assumptions -> AA > 0 && BB < 0 && CC > 0 && DD < 0 && EE > 0];
```

此时对Ferrari解和Cardano解的数值验证的结果也出现变化：

| Ferrari解                        | 退化的Ferrari解                     | Cardano解                        |
| ------------------------------- | ------------------------------- | ------------------------------- |
| ![[截屏2026-03-26 下午3.52.34.png]] | ![[截屏2026-03-26 下午3.53.06.png]] | ![[截屏2026-03-26 下午3.53.21.png]] |
这说明RN时空的Ferrari的第123个解可以在Q->0时退化到史瓦西时空的Cardano的第123个解。并且
- RN时空下，x应取第2根Root[-EE-DD #1+CC #1^2+BB #1^3+AA #1^4&,2]
- Schwartz-child时空下，x应取第2根Root[-EE-DD #1+CC #1^2+BB #1^3&,2]

2. 第二件事是关于$x_{0}$的求解结果，维持求解$x$时同样的参数设置，唯一的区别体现在常数项$u^{2}$的变化，可能引起判别式的变化从而解得root2为虚根（这可能是正常现象）。
对比$x_{0}$和$x$的Ferrari解如下:

| $x$                             | $x_{0}$                         |
| ------------------------------- | ------------------------------- |
| ![[截屏2026-03-26 下午3.52.34.png]] | ![[截屏2026-03-26 下午4.32.47.png]] |
虽然解的数值发生变化，但四个解的结构依旧类似。基于此，我们仍大胆选用$x_{0}=$root2`x0 = x0sols[[2, 1, 2]] // Normal; (*选物理根root2*)` ,进行后续对
xApprox的求解：
$$
\text{xApprox}=\sum _{k=1}^{\text{maxOrder}} \text{u2}^k \text{xk}(k)+\text{x0}
$$

3. 现在我们一直在做的事情是：反解
$$
u^{2}=p(x)=\left(b^2 x^2 \left(\epsilon ^2-1\right)-\epsilon ^2\right) \left(2 M x-Q^2 x^2-1\right):=\sum _{k=0}^4 a(k) x^k
$$
以得到$x(u^{2})$，但是$p(x)$在区间$\left[ 0, \frac{1}{r_{0}} \right]$上并不单调。这是四次函数$p(x)$关于$x$的图像，$r_{0}=10^{4}$
`testAK = {a[0] -> \[Epsilon]^2, a[1] -> DD, a[2] -> CC, a[3] -> BB,   a[4] -> AA}`

| global                               | local                                |
| ------------------------------------ | ------------------------------------ |
| ![[Pasted image 20260326164315.png]] | ![[Pasted image 20260326164321.png]] |
|                                      | *这个图要画清楚一些，重画*                       |

这意味着反解$u^{2}=p(x)$之后得到的$x(u^{2})$是一个关于$u^{2}$的多值函数
- 在$x\in(0,2\times 10^{-11}]$区间，对应$u^{2}\in[\text{u2small},\epsilon^{2})$
- 在$x\in[2\times 10^{-11},1/r_{0}]$区间，对应$u^{2}\in[\text{u2small},1]$
那么对于最终的偏折角积分
$$
\delta \phi=b\sqrt{ 1-\epsilon^{2} }\left(\int_{\epsilon^{2}}^{\text{u2small}}+\int_{\text{u2small}}^{1} \frac{2u\,J(u^{2})}{\sqrt{  1-u^{2}}} \,du^{2}\right)
$$
其中的雅可比因子
$$
J(u^{2})=\frac{2u\,d\left( \frac{1}{r} \right)}{2udu}=\frac{dx}{du^{2}}
$$
也必须分段处理:
$$
\frac{\delta \phi}{b\sqrt{ 1-\epsilon^{2} }}=\int_{\epsilon^{2}}^{\text{u2small}}\frac{2u\,J_{1}(u^{2})}{\sqrt{  1-u^{2}}}\,du^{2}+\int_{\text{u2small}}^{1} \frac{2u\,J_{2}(u^{2})}{\sqrt{  1-u^{2}}} \,du^{2}
$$
其中
$$
\begin{cases}
J_{1}(u)=\frac{dx_{1}}{du^{2}},\,u^{2}\in[\text{u2small},\epsilon^{2})\\
J_{2}(u)=\frac{dx_{2}}{du^{2}}，\,u^{2}\in[\text{u2small},1]
\end{cases}
$$
$x_{1}(u^{2})$对应$p(x)$极小值点左侧的反函数，$x_{2}(u^{2})$对应$p(x)$极小值点右侧的反函数。

4. 那么接下来对于反函数$x_{1}(u^{2}),x_{2}(u^{2})$的求解，也必须在$p(x)$的极小值点$(\text{xsmall,u2small})$两侧分类讨论。
先求$p(x)$极小值点，解析和数值的结果分别为：
数值极小值坐标：$(2.00003\times10^{-13},0.0000200039)$
解析极小值点 `xsmall=Root[a[1]+2 a[2] #1+3 a[3] #1^2+4 a[4] #1^3&,1]`

>[!bug] 待定系数法中xApprox的形式应随分支而变化
> 待定系数法把x写成u2在0处的展开式是不合理的，尝试在u2=$\epsilon^{2}$处展开再用待定系数法
> 那么对于右分支，$u^{2}<\epsilon^{2}$，关于$x_{0}$的方程化为
> $$x_{0}^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x_{0}^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x_{0}^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x_{0}+\underbrace{ 0 }_{\epsilon^{2}-\epsilon^{2} } =0\tag{1.1}$$
> 关于$x$的方程化为
> $$x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x + 0=u^{2}-\epsilon^{2}:=v_{2}$$
> where $v_{2}$ 是负数，取值范围是（-u2small,0]

![[截屏2026-03-26 下午6.42.34.png]]
让常数项$b_{0}$逐渐逼近0，数值求解$x_{0}$，可以看到root1，2共轭，其虚部和实部均会收敛至0，因此当$b_{0}=0$时，方程(1.1)退化为三次方程，root1和2退化为平凡解$x_{0}=0$ 
既然数值验证了x0=0，那么xAppox就没有零阶项。接下来用待定系数法求解xApprox，解得
$$
x=1.00011\times 10^{90} v^6-5.33505\times 10^{72} v^5+3.04924\times 10^{55} v^4-1.95192\times 10^{38} v^3+1.56187\times 10^{21} v^2-24995.1 v
$$

[[2026-04-02]]
各阶系数仍有振荡，并且数值验证的结果显示，收敛半径约为$10^{-18}$, 这就是$\Delta u^{2}=\epsilon^{2}-u2small$ 

![[截屏2026-04-02 下午12.00.10.png]]

所以这也就意味着极小值点u2small有奇性，在区间$(u_{2}small，1]$，需要用另一套系数展开$v$，或者用另一套级数$\sum x_{i}(u-1)^{i}$ .
这自然带来了在近日点$u^{2}=1$处展开$x(u^{2})$的需求，也就是[[2026-04-02#^43b726]]需要讨论的问题。
![[截屏2026-03-26 下午7.43.01.png|100]]
### [[2026-03-26]]组会纪要
1. 画图
2. 理论上如何从u2=1即近日点处展开
3. 旋转坐标轴，取反函数的左右分支的几个点，用$\sqrt{ x },x^{3/2}$等阶的幂级数去拟合，哪个效果好就用哪个当待定系数法的领头阶

---
>[!quote] 系数振荡发散的原因分析
>根据奇点分析（singularity analysis）和 Darboux 方法，反函数在复 y 平面上距离展开点 y₀ 最近的奇点（通常是分支点、极点或代数奇点）决定了收敛半径 R。 系数满足：
>$$\limsup_{k\to\infty} |xk[k]|^{1/k} = \frac{1}{R}$$
>- 如果 R < 1 → |xk[k]| 呈**指数增长**（(1/R)^k 因子主导）
>- 如果 R = 1 → 系数仍可能先增大，后期由奇点类型（如 1/k^{3/2}）缓慢衰减，但早期阶数已出现“越来越大”
>- 高阶导数 `derivs[[k]]`（k 较大）相对 p' 越大，R 就越小，xk 增长越快
>
>取 p(x) = x + x²（经典 Lagrange 反演）： 
反函数 $x(y) = \sum (-1)^{k+1} \frac{(2k)!}{k!(k+1)!} \frac{y^k}{4^{k-1}}$（Catalan 数相关），系数指数增长 + 振荡，R = 1/4

用帕德逼近可以削弱这种振荡：
![[Pasted image 20260326195054.png]]

>[!quote] 为什么Pade Approximation 不同 [m/n] 差别这么大？
>1. **Padé 逼近会“猜测”并放置极点（poles）** Padé 逼近把幂级数（你的 xApprox，即反函数的泰勒级数）转化为有理函数 P_m(v)/Q_n(v)。它最擅长捕捉**最近的极点**（pole），因为分母多项式 Q_n 的根会尽量去逼近函数在复平面上的最近奇点，从而让逼近在更大区域内有效。 但如果函数的最近奇点**不是极点，而是分支点（branch point）**（反函数 x(y) 经常出现代数分支点，例如 sqrt 类型或更复杂的多值情况），Padé 就只能用**虚假极点（spurious poles）**来“勉强”拟合。
    
2. **虚假极点（Spurious Poles）的出现和位置敏感性** 不同 {m,n} 会导致虚假极点的位置和数量差异很大：
    - 当分母阶数 n 较小（如 [5/1]，分母只有 1 次），它只能放很少的极点。如果级数暗示附近有一个“强”奇点，它可能会把这个虚假极点放得**离展开点 v=0 非常近**（甚至极近）。
    - 当你在 **v → -10^{-10}**（一个非常小的负值）处求值时，如果某个虚假极点正好落在 v ≈ -10^{-10} 附近或略微跨过它，函数值就会**爆炸式放大**（像 10¹⁴ 那样），或者符号剧烈翻转。
    - [4/4] 和 [3/5] 分母阶数更高，能更“均匀”地分布虚假极点，避免把单个极点放得太靠近原点，因此在你这个极小的 v 处表现得“平静”（接近 0）。
    
这就是不同mn的帕德逼近具有极端差异的根本原因：**Padé 对奇点结构的“猜测”对阶数选择非常敏感**，尤其在接近收敛半径边缘或有分支点时。
    
- 用 **反函数的幂级数**（Lagrange 反演得到的 xApprox），这类级数经常有**有限收敛半径**，且最近奇点多为**代数分支点**（而非简单极点）。Padé 对分支点的逼近效果通常不如对极点好。
- v = -1/10¹⁰ 是一个极小的扰动（testBK 里 u2 偏离了某个平衡值很小），说明你正工作在**奇点附近**或**收敛半径边界**。此时 Taylor 级数本身已发散（xk 越来越大），Padé 试图“拯救”它，但不同阶数“拯救”的方式完全不同。
- 真实值 ≈ 2×10^{-13}（很小实部 + 微小虚部）暗示真实反函数在这里接近 0 或某个小值，而 [5/1] 完全“失控”了。

---
四次函数在极小值点附近的 Puiseux 级数反函数

给定四次函数 $p(x)$，在点 $x = x_{\rm small}$ 处取得局部极小值，满足 $p'(x_{\rm small}) = 0$ 且 $p''(x_{\rm small}) > 0$。令 $u_{2,{\rm small}} = p(x_{\rm small})$，我们希望求反函数 $x(u_2)$ 在 $u_2$ 略大于 $u_{2,{\rm small}}$ 时的半解析表达式。

令 $t = x - x_{\rm small}$，$\delta = u_2 - u_{2,{\rm small}}$，则在极小值点附近有局部展开：
$$
\delta = A t^2 + B t^3 + C t^4 + O(t^5),
$$
其中 $A = \frac{p''(x_{\rm small})}{2} \approx 1.00018 \times 10^8$。

由于 $p'(x_{\rm small}) = 0$，普通 Taylor 反演失效，需使用带半整数幂的 Puiseux 级数：
$$
t = c_1 \sqrt{\delta} + c_2 \delta + c_3 \delta^{3/2} + c_4 \delta^2 + c_5 \delta^{5/2} + \cdots
$$

将上式代入局部展开并比较同次幂系数，可逐阶确定 $c_k$。其中最主导的一项为：
$$
c_1 = \frac{1}{\sqrt{2A}} \approx 7.0704314 \times 10^{-5}.
$$

因此，左右两侧的反函数前两项近似形式为：
$$
x_{\rm Right}(u_2) \approx x_{\rm small} + c_1 \sqrt{u_2 - u_{2,{\rm small}}},
$$
$$
x_{\rm Left}(u_2) \approx x_{\rm small} - c_1 \sqrt{u_2 - u_{2,{\rm small}}}.
$$

代入具体数值后得到：
$$
x_{\rm Right}(u_2) \approx 2.00003 \times 10^{-13} + 7.0704314 \times 10^{-5} \sqrt{u_2 - 0.0000200039},
$$
$$
x_{\rm Left}(u_2) \approx 2.00003 \times 10^{-13} - 7.0704314 \times 10^{-5} \sqrt{u_2 - 0.0000200039}.
$$

更高阶形式可通过 InverseSeries 自动生成，包含 $\delta$、$\delta^{3/2}$、$\delta^2$ 等项，适用于稍大范围的 $\delta$。该级数在极小值点附近精度最高，随着 $\delta$ 增大，就需更多阶项.

### [[2026-03-26]]组会纪要
1. 画图
2. 理论上如何从u2=1即近日点处展开
3. 旋转坐标轴，取反函数的左右分支的几个点，用$\sqrt{ x },x^{3/2}$等阶的幂级数去拟合，哪个效果好就用哪个当待定系数法的领头阶

---
[[2026-04-02]]
#### 1. 高精度绘图
之所以出现$u(x)$极小值点$10^{-11}$和$10^{-13}$不匹配的问题，是因为testparams在绘图时选用的是史瓦西组而非RN组。
对于史瓦西时空 绘制$u(r)=\sqrt{\left(1-\frac{2 M}{r}\right) \left(\frac{b^2}{r^2}+\epsilon ^2\right)}$关于$x=\frac{1}{r}$的函数图像

严格按照上面的绘图格式，绘制以下函数的总体和局部
```mathematica
Plot[(Sqrt[(1 - 2 M x) (b^2 x^2 + \[Epsilon]^2)]) /. testparams, {x, 
  0, 1/10^4}]
Plot[(Sqrt[(1 - 2 M x) (b^2 x^2 + \[Epsilon]^2)]) /. testparams, {x, 
  0, 1/10^10}]
```
其中`testparams = {M -> 1, r0 -> 10406, \[Epsilon] -> Sqrt[2239]/1020,  b -> 5301857/510};`
![[Pasted image 20260402114527.png|450]]

#### 2. 近日点展开

^43b726

理论上如何从u2=1即近日点处展开？
先试探可行域，再设法处理过拐点u2=u2small后引起的振荡发散，可以想到三种解决方案，可以一一尝试：
1. 左支用u2-$\epsilon^{2}$泰勒级数，右支用u2-1泰勒级数
2. 左支用u2-1的洛朗级数，右支用u2-1的泰勒级数
3. 左支用u2-1的帕德逼近，右支用u2-1的泰勒级数，但需要检证帕德逼近给出的真分式$\frac{m(u_{2})}{n(u_{2})}$对$u_{2}$得到的雅可比因子$J(u^{2})$是**真分式**
$$
J(u^{2})=\frac{2u\,d\left( \frac{1}{r} \right)}{2udu}=\frac{dx}{du^{2}}=\frac{m'(u_2)n(u_2)-m(u_2)n'(u_2)}{[n(u_2)]^2}
$$
对于这个真分式$J(u^{2})$，其与式$\frac{2u}{\sqrt{ 1-u^{2} }}$的积：
$$
\frac{2u\,J(u^{2})}{\sqrt{  1-u^{2}}} 
$$
的**可积性需要讨论**。

##### 右支泰勒展开
以上周为范例
>[!bug] u2=$\epsilon^{2}$处展开
> 待定系数法把x写成u2在0处的展开式是不合理的，尝试在u2=$\epsilon^{2}$处展开再用待定系数法
> 那么对于区间$u2small<u^{2}<\epsilon^{2}$，关于$x_{0}$的方程化为
> $$x_{0}^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x_{0}^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x_{0}^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x_{0}+\underbrace{ 0 }_{\epsilon^{2}-\epsilon^{2} } =0\tag{1.1}$$
> 关于$x$的方程化为
> $$x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x + 0=u^{2}-\epsilon^{2}:=v_{2}$$
> where $v_{2}$ 是负数，取值范围是（u2small-$\epsilon^{2}$,0]

类似地可以复刻代码和计算文档:
>[!bug]  $u_{2}=1$处展开
> 对于**单调区间**$1>u^{2}>u_{2}small$（极小值点右侧） ![[截屏2026-04-02 下午12.37.18.png|100]]
>关于$x$的四次方程:
>$$x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x+\epsilon ^2 =u^{2}$$
>等式两边-1
>$$x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x+(\epsilon ^2-1) =u^{2}-1:=w_{2}$$
>$w_{2}=0$时的方程即为关于$x_{0}$的方程：
>$$x_{0}^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x_{0}^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x_{0}^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x_{0}+(\epsilon ^2-1) =0$$
> 其中 $w_{2}$ 是负数，取值范围是（u2small-1,0]

区别于$u_{2}=\epsilon^{2}$处展开，u2=1时，方程不会退化为齐次方程，因此还需计算x0. 注意到u=1时是近日点，因此$x_{0}=\frac{1}{r_{0}}=10^{-4}$. 我们可以求解x0并验证根是否有物理意义

> [!Properties] `Solve`函数的机制
> 我们的目的是筛选有物理意义的x0解
> if add `&&x>0` into `Solve` function, we get
> ![[截屏2026-04-02 下午3.50.44.png]]
> 这是因为先构造 **Root 对象**，再尝试判断 x>0的逻辑与mma Root[...] 在没有数值化时，**符号上无法判断正负**的属性矛盾，导致Mathematica 返回 Undefined
> 正确的写法是取消逻辑上的x>0判断
![[截屏2026-04-02 下午3.49.19.png]]

用相同的办法求解$x_{0}$,xApprox=![[截屏2026-04-02 下午5.13.26.png|400]]

数值验证其在$u^{2}=u_{2}small$处是否收敛到解析解$x(u_{2})$:
![[截屏2026-04-02 下午5.15.05.png]]

在区间$[u_{2}small-1,0]$绘图：
![[截屏2026-04-02 下午5.16.09.png]]

观察到级数xApprox只在部分区间拟合好，说明 **展开点到最近奇点的距离限制了收敛半径**。也就是存在一个最近奇点 $w_{2,\text{sing}}$，级数只在$w_2 \in (w_{2,\text{sing}},\,0)$内有效。

---

xApprox[w2] 来自：方程 $p_3(x)-w_{2}=0$在 $w_2=0$ 处对隐函数展开. 隐函数级数失效的条件是：
 $$p_3(x)-w_{2}=0,\,\partial_x p_3(x)=0$$
 这样解出来，奇点就是右支区间的下界:![[截屏2026-04-02 下午6.57.18.png]]

##### 方案三：近日点展开式的帕德逼近
测试完毕，效果一般，帕德逼近只是让泰勒级数更弯曲，还是没办法从根本上解决奇点导数发散的问题。
![[截屏2026-04-02 下午6.58.35.png]]

##### 方案四：尝试提取包含奇点的根式$\sqrt{ u^{2}-u_{2}small }$,对剩下的部分$g(u^{2})=\frac{x(u^{2})}{\sqrt{ u^{2}-u_{2}small }}$作展开
>[!bug]  $u_{2}=u2small$处展开
> 对于右支**单调区间**$1>u^{2}>u_{2}small$（极小值点右侧） ![[截屏2026-04-02 下午12.37.18.png|100]]
>关于$x$的四次方程:
>$$x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x+\epsilon ^2 =u^{2}$$
>等式两边-u2small
>$$x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x+(\epsilon ^2-u_{2}small) =u^{2}-u_{2}small:=z_{2}$$
>$z_{2}=0$时的方程即为关于$x_{0}$的方程：
>$$x_{0}^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x_{0}^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x_{0}^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x_{0}+(\epsilon ^2-u_{2}small) =0$$
> 其中 $z_{2}$ 是正数，取值范围是[0,1-u2small]

##### 方案五（已成功，见下文）：待定系数法的试探解也应该是半整数幂级数
于是自然来到上周组会提出的idea：旋转坐标轴，取反函数的左右分支的几个点，用$\sqrt{ x },x^{3/2}$等阶的幂级数去拟合，哪个效果好就用哪个当待定系数法的领头阶

---
[[2026-04-08]]
### 解的多种示意图绘制
#### 复平面内解$x$的分布
求解方程$p(x)=u^{2}$, 在区间$[-5,5]$中扫描$u^{2}$ （这里我们只关心数学，不关心$u^{2}<0$物不物理)，可以得到一系列复根$x$的序列，在复平面中把它们画出来。

| global$\quad \quad$                  | $\mathrm{Re}x\in[-2\times10^{-8},10^{-8}]$ | $x\in[0,3\times 10^{-4}]$            | $x\in[0.5,13]$                       |
| ------------------------------------ | ------------------------------------------ | ------------------------------------ | ------------------------------------ |
| ![[Pasted image 20260408113640.png]] | ![[Pasted image 20260408113858.png\|134]]  | ![[Pasted image 20260408114447.png]] | ![[Pasted image 20260408114503.png]] |
关注$u=-2,-1,1,2$时，四个数值解的行为，可以发现前两个小根对u2敏感，后两个大根对u2不敏感，且第四根（x=11.978）尤甚，u2增加1，其仅变化了10^-10。这导致了复平面内画根的分布时，出现了散点凝聚在实轴0.52，11.98附近的现象。
![[截屏2026-04-08 下午12.07.01.png]]
绘制$u^{2}=p(x),x\in[0,13]$的函数图像，也可以发现$p'(11.98)>p'(0.52)\gg p'(-0.001\to{0}.001)$
![[Pasted image 20260408121510.png|300]]![[Pasted image 20260409080031.png|300]]![[Pasted image 20260408121534.png|300]]![[Pasted image 20260402114527.png|300]] 

这也就解释了为什么不同的根敏感性不同。同时该四次函数$p(x)$有一个极大值，两个极小值点，由此可得出使方程有四个实数根的$u^{2}$取值范围: $u^{2}\in [u_{2}small,3.9\times 10^{6}]$。 如果要观察复根x的行为，则$u^{2}$取值应该在此区间外。
   
#### $x(u^{2})$ 的绘制
结合上文结论，为了方便观察分支行为，可采用x的关键区间放大，无关区间压缩的办法画图。期望效果应该是这样：
![[IMG_AF06994C43F3-1.jpg]]

#### x在复平面中时，u2的图像
##### 幅角图ComplexPlot
当$u^{2}=-1,-2$时，在实数域其与$p(x)$的图像仅有两个交点，因此前两个根必为虚根
```mathematica
In[2166]:=  NSolve[poly == -1, x, WorkingPrecision -> 16]
 NSolve[poly == -2, x, WorkingPrecision -> 16]

Out[2166]= {{x -> -9.99779922132*10^-9 - 
    0.00009998999696072122 I}, {x -> -9.99779922132*10^-9 + 
    0.00009998999696072122 I}, {x -> 0.5217804013388130}, {x -> 
   11.97821961865679}}

Out[2167]= {{x -> -1.999599694519*10^-8 - 
    0.00014140791356836761 I}, {x -> -1.999599694519*10^-8 + 
    0.00014140791356836761 I}, {x -> 0.5217804213732246}, {x -> 
   11.97821961861877}}
```

以此为参考，画幅角图
![[Pasted image 20260409160458.png]]
![[Pasted image 20260409160515.png|250]]![[Pasted image 20260409160538.png|250]]

### u=u2small处展开为Puiseux级数
#### 右支
试探解的形式为准级数：
$$
\text{xApprox}=x_{0}+\sum _{k=1}^{\text{maxOrder}}x_{k} v_{2}^{(2k-1)/2} 
$$
where $v_{2}=u^{2}-u_{2small}$（在分支点一侧）。
这种情况下能用待定系数法求解关于$x$的四次方程
>$$x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x+\epsilon ^2 =u^{2}$$
>等式两边-u2small
>$$x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x+(\epsilon ^2-u_{2}small) =u^{2}-u_{2}small:=v_{2}$$
>$v_{2}=0$时的方程即为关于$x_{0}$的方程：
>$$x_{0}^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x_{0}^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x_{0}^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x_{0}+(\epsilon ^2-u_{2}small) =0$$
> 其中 $v_{2}$ 是正数，取值范围是（0,1-u2small]

通过对试探解xApprox求导并回代方程组的办法，可以得到关于系数$x_{k}$的一系列线性方程组。

---
在 $v_2 = u^2 - u_{2\text{small}}$ 处存在平方根型分支点奇性，我们采用以下**只含奇次半整数幂**的试探解形式：

$$
x_{\rm Approx} = x_0 + \sum_{k=1}^{4} x_k \, v_2^{(2k-1)/2}
$$
即
$$
x \approx x_0 + x_1 \sqrt{v_2} + x_3 \, v_2^{3/2} + x_5 \, v_2^{5/2} + x_7 \, v_2^{7/2} + \cdots
$$
各阶系数通过待定系数法得到（$xk[1]$ 取正号），递推公式如下：
$$
x_1 = \sqrt{\dfrac{2}{p''(x_0)}}
$$
$$
x_3 = -\dfrac{p'''(x_0) \, x_1^3}{6 \, p''(x_0)}
$$
$$
x_5 = -\dfrac{1}{p''(x_0)} \left[ 
      \dfrac{3}{2} p'''(x_0) \, x_1^2 x_3^2 
    + p^{(4)}(x_0) \, x_1^4 x_3 
    + \dfrac{1}{24} p^{(5)}(x_0) \, x_1^5 
    \right]
$$
$$
x_7 = -\dfrac{1}{p''(x_0)} \Bigg[ 
      \dfrac{1}{2} p'''(x_0) \bigl(15 x_1^2 x_5 x_3 + 5 x_1 x_3^3\bigr) \\
    \quad + p^{(4)}(x_0) \bigl(x_1^4 x_5 + 6 x_1^3 x_3^2\bigr) \\
    \quad + \dfrac{5}{2} p^{(5)}(x_0) \, x_1^4 x_3 
    + \dfrac{1}{120} p^{(6)}(x_0) \, x_1^5 
    \Bigg]
$$

其中 $p^{(m)}(x_0)$ 表示多项式 $p_4(x)$ 在 $x = x_0$ 处的第 $m$ 阶导数（即 `derivs[[m]]`）。

 注意$x_1$ 取正号，对应 $\sqrt{v_2}$ 的正分支，这意味着我们求出来的Puiseux级数是对红色部分的拟合，$u^{2}\in[u_{2}small,1]$。使用时需先求出满足 $p_4(x_0) = u_{2\text{small}}$ 的 $x_0$，(其实就是$p(x)$的极小值点$x_{0}=2\times 10^{-13}$）再依次代入以上公式计算各系数。

##### 待定系数法的计算细节
```mathematica
y = 0; 
For[k = 1, k <= maxOrder, k++, y += xk[k] s^k; (* 逐步加入当前项 *)
temp = p4[x0 + y] - u2small - s^2; (* 残差方程 *)

(* 展开到 s^k 阶，提取系数 *) 
tempSeries = Series[temp, {s, 0, k}] // Normal; coeffK = Coefficient[tempSeries, s^k];

(* 求解当前 xk[k]（线性方程，只有一个未知数） *) 
sol = Solve[coeffK == 0, xk[k]]; xk[k] = sol[[1, 1, 2]]; (* 默认取第一个解 *) (* 如果需要另一支，可改成 sol[[2,1,2]] 或手动选择正/负号 *) 
];
xApprox = x0 + y /. s -> Sqrt[v2];
```

为了求解最低阶系数 $x_1$，我们先取**最低阶近似**：
$$
x \approx x_0 + x_1 \sqrt{v_2}
$$
令 $s = \sqrt{v_2}$（当 $v_2 \to 0^+$ 时，$s \to 0^+$），则上式变为：
$$
x \approx x_0 + x_1 s
$$
将其代入原方程 $p_4(x) = u^2$，并注意到 $u^2 = u_{2\text{small}} + v_2 = u_{2\text{small}} + s^2$，得到：
$$
p_4(x_0 + x_1 s) = u_{2\text{small}} + s^2
$$
在 $x_0$ 处进行 Taylor 展开并匹配系数：
$$
p_4(x_0 + x_1 s) = p_4(x_0) + p'(x_0)(x_1 s) + \frac{p''(x_0)}{2}(x_1 s)^2 + \frac{p'''(x_0)}{6}(x_1 s)^3 + \frac{p^{(4)}(x_0)}{24}(x_1 s)^4 + O(s^5)
$$
已知 $p_4(x_0) = u_{2\text{small}}$ 且 $p'(x_0) = 0$，上式简化为：
$$
p_4(x_0 + x_1 s) = u_{2\text{small}} + \frac{p''(x_0)}{2} x_1^2 s^2 + \frac{p'''(x_0)}{6} x_1^3 s^3 + \frac{p^{(4)}(x_0)}{24} x_1^4 s^4 + O(s^5)
$$
令其等于右边：
$$
u_{2\text{small}} + \frac{p''(x_0)}{2} x_1^2 s^2 + \frac{p'''(x_0)}{6} x_1^3 s^3 + \frac{p^{(4)}(x_0)}{24} x_1^4 s^4 + O(s^5) = u_{2\text{small}} + s^2
$$
两边同时减去 $u_{2\text{small}}$，得到：
$$
\frac{p''(x_0)}{2} x_1^2 s^2 + \frac{p'''(x_0)}{6} x_1^3 s^3 + \frac{p^{(4)}(x_0)}{24} x_1^4 s^4 + O(s^5) = s^2
$$
比较 $s^2$ 项的系数，可得：
$$
\frac{p''(x_0)}{2} x_1^2 = 1
$$
解得：

$$
x_1^2 = \frac{2}{p''(x_0)}
$$
由于要求 $x_1$ 取**正号**，因此：

$$
x_1 = \sqrt{\frac{2}{p''(x_0)}}
$$
**说明**：更高阶的 $s^3$、$s^4$ 等项在求 $x_1$ 时可以暂时忽略，因为当 $s \to 0$ 时，它们远小于 $s^2$ 项。它们将在后续推导 $x_3$、$x_5$ 等系数时才被用于匹配更高次幂。

##### 近似拟合效果

解xApprox：$(v_{2}=u^{2}-u_{2}small)$
$$
-\frac{1}{3} \sqrt{2} \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right) \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{5/2} \text{v2}^{3/2}+\text{x0}+\sqrt{2} \sqrt{\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}} \sqrt{\text{v2}}-\frac{\text{v2}^{5/2} \left(\frac{2 \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)^3}{3 \left(2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2\right)^6}-32 \sqrt{2} b^2 Q^2 \left(1-\epsilon ^2\right) \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right) \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{9/2}\right)}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}-\frac{\text{v2}^{7/2} \left(24 b^2 \left(\frac{8}{3} \sqrt{2} \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)^2 \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{13/2}-\frac{4 \left(\frac{2 \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)^3}{3 \left(2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2\right)^6}-32 \sqrt{2} b^2 Q^2 \left(1-\epsilon ^2\right) \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right) \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{9/2}\right)}{\left(2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2\right)^3}\right) \left(1-\epsilon ^2\right) Q^2+\frac{1}{2} \left(10 \sqrt{2} \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right) \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{9/2} \left(\frac{2 \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)^3}{3 \left(2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2\right)^6}-32 \sqrt{2} b^2 Q^2 \left(1-\epsilon ^2\right) \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right) \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{9/2}\right)-\frac{20 \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)^3}{27 \left(2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2\right)^8}\right) \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)\right)}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}
$$

代入testparamsRN有
$$
x_{Approx\,R}=2.00003\times10^{-13} + 0.000099991 \sqrt{ u^{2}-u_{2}small } + 9.9973\times10^{-13} (u^{2}-u_{2}small)^{\frac{3}{2}} - 
 1.91789\times10^{-28} (u^{2}-u_{2}small)^{\frac{5}{2}}- 1.15092\times10^{-35} (u^{2}-u_{2}small)^{\frac{7}{2}}
$$

![[截屏2026-04-09 下午7.40.05.png|424]]

#### 左支
试探解的形式为准级数：
$$
\text{xApprox}=x_{0}+\sum _{k=1}^{\text{maxOrder}}x_{k} w_{2}^{(2k-1)/2} 
$$
where $u^{2}-u_{2small}=w_{2}$（在分支点左侧）。
这种情况下能用待定系数法求解关于$x$的四次方程
>$$x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x+\epsilon ^2 =u^{2}$$
>等式两边-u2small
>$$x^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x+(\epsilon ^2-u_{2}small) =u^{2}-u_{2}small:=w_{2}$$
>$w_{2}=0$时的方程即为关于$x_{0}$的方程：
>$$x_{0}^4 \left(b^2 Q^2-b^2 Q^2 \epsilon ^2\right)+x_{0}^3 \left(2 b^2 M \epsilon ^2-2 b^2 M\right)+x_{0}^2 \left(-b^2 \epsilon ^2+b^2+Q^2 \epsilon ^2\right)-2 M  \epsilon ^2 x_{0}+(\epsilon ^2-u_{2}small) =0$$
> 其中 $w_{2}$ 是正数，取值范围是$[0,\epsilon^{2}-u_{2}small)$
通过对试探解xApprox求导并回代方程组的办法，可以得到关于系数$x_{k}$的一系列线性方程组。

在 $w_2 = u^{2}-u_{2\text{small}}=0$ 处存在平方根型分支点奇性，我们采用以下**只含奇次半整数幂**的试探解形式：
$$
x_{\rm Approx} = x_0 + \sum_{k=1}^{4} x_k \, w_2^{(2k-1)/2}
$$
即
$$
x \approx x_0 + x_1 \sqrt{w_2} + x_3 \, w_2^{3/2} + x_5 \, w_2^{5/2} + x_7 \, w_2^{7/2} + \cdots
$$
各阶系数通过待定系数法得到（**左支$xk[1]$ 取负号**），递推公式如下：
$$
x_1 = -\sqrt{\dfrac{2}{p''(x_0)}}
$$
$$
x_3 = -\dfrac{p'''(x_0) \, x_1^3}{6 \, p''(x_0)}
$$
$$
x_5 = -\dfrac{1}{p''(x_0)} \left[ 
      \dfrac{3}{2} p'''(x_0) \, x_1^2 x_3^2 
    + p^{(4)}(x_0) \, x_1^4 x_3 
    + \dfrac{1}{24} p^{(5)}(x_0) \, x_1^5 
    \right]
$$
$$
x_7 = -\dfrac{1}{p''(x_0)} \Bigg[ 
      \dfrac{1}{2} p'''(x_0) \bigl(15 x_1^2 x_5 x_3 + 5 x_1 x_3^3\bigr) \\
    \quad + p^{(4)}(x_0) \bigl(x_1^4 x_5 + 6 x_1^3 x_3^2\bigr) \\
    \quad + \dfrac{5}{2} p^{(5)}(x_0) \, x_1^4 x_3 
    + \dfrac{1}{120} p^{(6)}(x_0) \, x_1^5 
    \Bigg]
$$

其中 $p^{(m)}(x_0)$ 表示多项式 $p_4(x)$ 在 $x = x_0$ 处的第 $m$ 阶导数（即 `derivs[[m]]`）。

 注意$x_1$ 取负号，对应 $\sqrt{w_2}$ 的左分支，这意味着我们求出来的Puiseux级数是对绿色部分的拟合，$u^{2}\in(0,u_{2}small]$。使用时需先求出满足 $p_4(x_0) = u_{2\text{small}}$ 的 $x_0$，(其实就是$p(x)$的极小值点$x_{0}=2\times 10^{-13}$）再依次代入以上公式计算各系数。

##### 待定系数法的计算细节
为了求解最低阶系数 $x_1$，我们先取**最低阶近似**：
$$
x \approx x_0 + x_1 \sqrt{w_2}
$$
令 $s = \sqrt{w_2}$（当 $v_2 \to 0^+$ 时，$s \to 0^+$），则上式变为：
$$
x \approx x_0 + x_1 s
$$
将其代入原方程 $p_4(x) = u^2$，并注意到 $u^2 = u_{2\text{small}} +w_2 = u_{2\text{small}} + s^2$，得到：
$$
p_4(x_0 + x_1 s) = u_{2\text{small}} +s^2
$$
在 $x_0$ 处进行 Taylor 展开并匹配系数，比较 $s^2$ 项的系数，仍可得与右支一致的结果：
$$
x_1^2 = \frac{2}{p''(x_0)}
$$
唯一的区别在于符号，由于要求 $x_1$ 取**负号**，因此：
$$
x_1 = -\sqrt{\frac{2}{p''(x_0)}}
$$
那么高次系数$x_{3},x_{5},\dots$也会随着$x_{1}$而改变
#### 准级数在左支的拟合效果
>[!bug] 数值解绘图bug：NxLeft[v2] 在单点代入时正确，但 Plot 出来的曲线却完全不同
>![[截屏2026-04-10 上午1.29.18.png|300]]
>故障原因和前面回代x的解析解试图回归到0，$\epsilon^{2}-u_{2}small$ 时却失败的原因是类似的
>`Table` + `ListPlot` 强制使用了指定的高精度，而 `Plot` 在自动取样时回退到了标准机器精度，导致了数值不稳定。如果希望直接使用 `Plot` 得到和 `ListPlot` 一样准确的结果，就得手动提升 `Plot` 的计算精度：
>```mathematica
>Plot[NxLeft[v2], {v2, 0, v2max}, WorkingPrecision -> 30]
>```
>这就成功了
>![[Pasted image 20260410013244.png|300]]

解xApproxLeft：
$$
\frac{1}{3} \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{5/2} \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right) \sqrt{2} \text{v2}^{3/2}-\frac{\left(32 b^2 Q^2 \left(1-\epsilon ^2\right) \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right) \sqrt{2} \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{9/2}+\frac{2 \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)^3}{3 \left(2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2\right)^6}\right) \text{v2}^{5/2}}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}+\text{x0}-\sqrt{2} \sqrt{\text{v2}} \sqrt{\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}}-\frac{\text{v2}^{7/2} \left(24 b^2 \left(\frac{1}{3} (-8) \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)^2 \sqrt{2} \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{13/2}-\frac{4 \left(32 b^2 Q^2 \left(1-\epsilon ^2\right) \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right) \sqrt{2} \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{9/2}+\frac{2 \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)^3}{3 \left(2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2\right)^6}\right)}{\left(2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2\right)^3}\right) \left(1-\epsilon ^2\right) Q^2+\frac{1}{2} \left(-10 \sqrt{2} \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right) \left(32 b^2 Q^2 \left(1-\epsilon ^2\right) \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right) \sqrt{2} \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{9/2}+\frac{2 \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)^3}{3 \left(2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2\right)^6}\right) \left(\frac{1}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}\right)^{9/2}-\frac{20 \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)^3}{27 \left(2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2\right)^8}\right) \left(24 Q^2 \text{x0} \left(1-\epsilon ^2\right) b^2+12 M \left(\epsilon ^2-1\right) b^2\right)\right)}{2 \left(\frac{Q^2 \epsilon ^2}{b^2}-\epsilon ^2+1\right) b^2+12 Q^2 \text{x0}^2 \left(1-\epsilon ^2\right) b^2+12 M \text{x0} \left(\epsilon ^2-1\right) b^2}
$$
代入testparamsRN得数值系数的级数：
$$
x_{\mathrm{Approx}\,L}(u)

= 2.00003\times10^{-13}

- 0.000099991\,\sqrt{u^2-u_{2\mathrm{small}}}

- 9.9973\times10^{-13}(u^2-u_{2\mathrm{small}})^{3/2}

+ 1.91969\times10^{-28}(u^2-u_{2\mathrm{small}})^{5/2}

+ 1.15122\times10^{-35}(u^2-u_{2\mathrm{small}})^{7/2}.
$$
可以看到高阶项系数不同于$x_{Approx\,R}$, 这也就的确意味着四次函数$p(x)$在极小值点$x_{0}$两侧不对称。
$$
x_{Approx\,R}=2.00003\times10^{-13} + 0.000099991 \sqrt{ u^{2}-u_{2}small } + 9.9973\times10^{-13} (u^{2}-u_{2}small)^{\frac{3}{2}} - 
 1.91789\times10^{-28} (u^{2}-u_{2}small)^{\frac{5}{2}}- 1.15092\times10^{-35} (u^{2}-u_{2}small)^{\frac{7}{2}}
$$
![[截屏2026-04-10 上午1.37.24.png|365]]误差的绝对值随着$v2=u2-u2small$增大而增大，这是正常的，因为展开中心是 $u2=u2small$,随着$x$远离极小值点$x0$而靠近$0$，$u^{2}=p(x)$也会越来越远离 $u2small$，截断级数总会带来这样的误差。

[[2026-04-10]]
在[[2026-03-26]]撰写计算文档中有这样一段，这就是我们今天计算的工作思路。
>[!quote] 工作思路
>反解$u^{2}=p(x)$之后得到的$x(u^{2})$是一个关于$u^{2}$的多值函数
>- 在$x\in(0,2\times 10^{-11}]$区间，对应左分支$u^{2}\in[\text{u2small},\epsilon^{2})$
>- 在$x\in[2\times 10^{-11},1/r_{0}]$区间，对应右分支$u^{2}\in[\text{u2small},1]$
>那么对于最终的偏折角积分
>$$\delta \phi=b\sqrt{ 1-\epsilon^{2}}\left(\int_{\epsilon^{2}}^{\text{u2small}}+\int_{\text{u2small}}^{1}\frac{J(u^{2})}{\sqrt{  1-u^{2}}} \,du^{2}\right)$$
>其中的雅可比因子
>$$J(u^{2})=\frac{\,d\left( \frac{1}{r} \right)}{2udu}=\frac{dx}{du^{2}}$$
>也必须分段处理:
>$$\frac{\delta \phi}{b\sqrt{ 1-\epsilon^{2}}}=\int_{\epsilon^{2}}^{\text{u2small}}\frac{\,J_{L}(u^{2})}{\sqrt{  1-u^{2}}}\,du^{2}+\int_{\text{u2small}}^{1} \frac{\,J_{R}(u^{2})}{\sqrt{  1-u^{2}}}\,du^{2}$$
>其中
>$$\begin{cases}J_{L}(u)=\frac{dx_{L}}{du^{2}},\,u^{2}\in[\text{u2small},\epsilon^{2})\\
J_{R}(u)=\frac{dx_{R}}{du^{2}}，\,u^{2}\in[\text{u2small},1]\end{cases}$$
>$x_{L}(u^{2})$对应$p(x)$极小值点左侧的反函数（绿色部分），$x_{R}(u^{2})$对应$p(x)$极小值点右侧的反函数（红色部分）。

四次方程的系数如下，以下所有AA BB CC DD都是指代这里的ABCD
$$
\begin{equation}\begin{aligned}
&A=b^{2}Q^{2}(1-\epsilon^{2})\sim r^{4}\\
&B=2b^{2}M(1-\epsilon^{2})\sim r^{3}\\
&C=b^{2}\left( 1-\epsilon^{2}+\frac{Q^{2}}{b^{2}}\epsilon^{2} \right)\sim r^{2}\\
&D=-2M\epsilon^{2} \sim r\\
&E=(\epsilon^{2}-u^{2})
\end{aligned}\end{equation}
$$
那么半解析形式的左支反函数$x_{L}(u^{2})$如下
$$
\begin{aligned}
x_{L}(v_{2}=u^{2}-u_{2}small)
&= x_0 \\[4pt]
&\quad - \sqrt{2}\,\sqrt{v_2}
\sqrt{\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}} \\[6pt]
&\quad + \frac{\sqrt{2}}{3} v_2^{3/2}
(6 BB + 24 AA x_0)
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{5/2} \\[6pt]
&\quad - \frac{v_2^{5/2}}{2 CC + 6 BB x_0 + 12 AA x_0^2}
\Bigg[
32\sqrt{2} AA (6 BB + 24 AA x_0)
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{9/2} \\
&\qquad\qquad
+ \frac{2(6 BB + 24 AA x_0)^3}
{3(2 CC + 6 BB x_0 + 12 AA x_0^2)^6}
\Bigg] \\[6pt]
&\quad - \frac{v_2^{7/2}}{2 CC + 6 BB x_0 + 12 AA x_0^2}
\Bigg[
\frac{1}{2}(6 BB + 24 AA x_0) \\
&\qquad\qquad \times
\Bigg(
-\frac{20(6 BB + 24 AA x_0)^3}
{27(2 CC + 6 BB x_0 + 12 AA x_0^2)^8} \\
&\qquad\qquad
- 10\sqrt{2}(6 BB + 24 AA x_0)
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{9/2} \\
&\qquad\qquad \times
\Bigg(
32\sqrt{2} AA (6 BB + 24 AA x_0)
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{9/2} \\
&\qquad\qquad\qquad
+ \frac{2(6 BB + 24 AA x_0)^3}
{3(2 CC + 6 BB x_0 + 12 AA x_0^2)^6}
\Bigg)
\Bigg) \\
&\qquad\qquad
+ 24 AA
\Bigg(
-\frac{8}{3}\sqrt{2}(6 BB + 24 AA x_0)^2
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{13/2} \\
&\qquad\qquad
- \frac{
4\Big(
32\sqrt{2} AA (6 BB + 24 AA x_0)
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{9/2}
+ \frac{2(6 BB + 24 AA x_0)^3}
{3(2 CC + 6 BB x_0 + 12 AA x_0^2)^6}
\Big)}
{(2 CC + 6 BB x_0 + 12 AA x_0^2)^3}
\Bigg)
\Bigg]
\end{aligned}
$$
对$u^{2}$也就是$v_{2}$求导，即得左支雅可比因子
$$
\begin{aligned}
dx_{L}/du^{2}
&= -\frac{1}{\sqrt{2}\,\sqrt{v_2}}
\sqrt{\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}} \\[6pt]
&\quad + \frac{\sqrt{v_2}}{\sqrt{2}}
(6 BB + 24 AA x_0)
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{5/2} \\[6pt]
&\quad - \frac{5 v_2^{3/2}}{2(2 CC + 6 BB x_0 + 12 AA x_0^2)}
\Bigg[
32\sqrt{2} AA (6 BB + 24 AA x_0)
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{9/2} \\
&\qquad\qquad
+ \frac{2(6 BB + 24 AA x_0)^3}
{3(2 CC + 6 BB x_0 + 12 AA x_0^2)^6}
\Bigg] \\[6pt]
&\quad - \frac{7 v_2^{5/2}}{2(2 CC + 6 BB x_0 + 12 AA x_0^2)}
\Bigg[
\frac{1}{2}(6 BB + 24 AA x_0) \\
&\qquad\qquad \times
\Bigg(
-\frac{20(6 BB + 24 AA x_0)^3}
{27(2 CC + 6 BB x_0 + 12 AA x_0^2)^8} \\
&\qquad\qquad
- 10\sqrt{2}(6 BB + 24 AA x_0)
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{9/2} \\
&\qquad\qquad \times
\Bigg(
32\sqrt{2} AA (6 BB + 24 AA x_0)
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{9/2} \\
&\qquad\qquad\qquad
+ \frac{2(6 BB + 24 AA x_0)^3}
{3(2 CC + 6 BB x_0 + 12 AA x_0^2)^6}
\Bigg)
\Bigg) \\
&\qquad\qquad
+ 24 AA
\Bigg(
-\frac{8}{3}\sqrt{2}(6 BB + 24 AA x_0)^2
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{13/2} \\
&\qquad\qquad
- \frac{
4\Big(
32\sqrt{2} AA (6 BB + 24 AA x_0)
\left(\frac{1}{2 CC + 6 BB x_0 + 12 AA x_0^2}\right)^{9/2}
+ \frac{2(6 BB + 24 AA x_0)^3}
{3(2 CC + 6 BB x_0 + 12 AA x_0^2)^6}
\Big)}
{(2 CC + 6 BB x_0 + 12 AA x_0^2)^3}
\Bigg)
\Bigg]
\end{aligned}
$$

数值积分这两项，左支数值积分：
$$\text{Int}_{L}=\int_{\epsilon^{2}}^{\text{u2small}}\frac{\,J_{L}(u^{2})}{\sqrt{  1-u^{2}}}\,du^{2}=2.000050004902\times10^{-13}$$
右支数值积分：
$$
\text{Int}_{R}=\int_{\text{u2small}}^{1} \frac{\,J_{R}(u^{2})}{\sqrt{  1-u^{2}}}\,du^{2}=0.000157065497
$$
可以看到右支积分的数量级远大于左支，这也符合近日点附近对偏折角有主要贡献的物理图像。

基于此，可以数值计算有Puiseux级数近似而来的偏折角：
$$
\delta \phi_{\text{Approx}}=2b \sqrt{1-\text{$\epsilon $2}} (\Re(\text{Int}_{L})+\Re(\text{Int}_{R}))-\pi=5.11\times 10^{-8}
$$
$$
\delta \phi=2\int_{10^{4}}^{\infty}\frac{b \sqrt{1-\epsilon ^2}}{r \sqrt{r^2-\frac{\left(r^2 \epsilon ^2-b^2 \left(\epsilon ^2-1\right)\right) \left(r (r-2 M)+Q^2\right)}{r^2}}}\,dr-\pi=0.000400078
$$
$$
\frac{4M}{b}=0.00039996
$$
---
[[2026-04-15]] 
## 数值误差过大分析
排查数值residual由谁贡献，是右支积分的Puiseux近似精度不够导致，还是变量代换x->u本身的局限导致

为此，我们需要直接数值计算变量代换u(x)=u的右支解析解，并数值计算右支雅可比因子，再数值积分右支，得到解析的右支数值。将其与刚才的IntRightNum对比，即可排查是否是Puiseux的精度不够导致residual.李冠鸣检查后发现是由于Puiseux级数抛弃了整数阶

### 完整的Puiseux级数应“有零有整”
补全整数阶后可以恢复，且可以展开为合理的二重级数
*李冠鸣负责的部分*

### 文献调研
《Exploring plasma and dark matter on photon deflection by Reissner–Nordström black hole with scalar hair and its shadow》

![[6351dc24a05d8bb02cf6a86226e11f31.png]]
![[54a16897709596431de88e535dd405e6.png]]
![[901926551bc1973b159ea83c687058f2.png]]