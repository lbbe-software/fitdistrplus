# Starting values used in fitdistrplus

We denote by the raw empirical moment by
$$m_{j} = \frac{1}{n}\sum\limits_{i = 1}^{n}x_{i}^{j},$$ by the centered
empirical moment by
$$\mu_{j} = \frac{1}{n}\sum\limits_{i = 1}^{n}\left( x_{i}^{j} - m_{1} \right).$$
Starting values are computed in `R/util-startarg.R`. We give below the
starting values for discrete and continuous distributions and refer to
the bibliograhy sections for details.

## 1. Discrete distributions

### 1.1. Base R distribution

#### 1.1.1. Geometric distribution

The MME is used $\widehat{p} = 1/\left( 1 + m_{1} \right)$.

#### 1.1.2. Negative binomial distribution

The MME is used
$\widehat{n} = m_{1}^{2}/\left( \mu_{2} - m_{1} \right)$.

#### 1.1.3. Poisson distribution

Both the MME and the MLE is $\widehat{\lambda} = m_{1}$.

#### 1.1.4. Binomial distribution

The MME is used
$$\left. Var\lbrack X\rbrack/E\lbrack X\rbrack = 1 - p\Rightarrow\widehat{p} = 1 - \mu_{2}/m_{1}. \right.$$
the size parameter is
$$\widehat{n} = \lceil\max\left( \max\limits_{i}x_{i},m_{1}/\widehat{p} \right)\rceil.$$

### 1.2. logarithmic distribution

The expectation simplifies for small values of
$p$$$E\lbrack X\rbrack = - \frac{1}{\log(1 - p)}\frac{p}{1 - p} \approx - \frac{1}{- p}\frac{p}{1 - p} = \frac{1}{1 - p}.$$
So the initial estimate is $$\widehat{p} = 1 - 1/m_{1}.$$

### 1.3. Zero truncated distributions

This distribution are the distribution of $X|X > 0$ when $X$ follows a
particular discrete distributions. Hence the initial estimate are the
one used for base R on sample $x - 1$.

### 1.4. Zero modified distributions

The MLE of the probability parameter is the empirical mass at 0
${\widehat{p}}_{0} = \frac{1}{n}\sum_{i}1_{x_{i} = 0}$. For other
estimators we use the classical estimator with probability parameter
$1 - {\widehat{p}}_{0}$.

### 1.5. Poisson inverse Gaussian distribution

The first two moments are
$$E\lbrack X\rbrack = \mu,Var\lbrack X\rbrack = \mu + \phi\mu^{3}.$$ So
the initial estimate are
$$\widehat{\mu} = m_{1},\widehat{\phi} = \left( \mu_{2} - m_{1} \right)/m_{1}^{3}.$$

## 2. Continuous distributions

### 2.1. Normal distribution

The MLE is the MME so we use the empirical mean and variance.

### 2.2. Lognormal distribution

The log sample follows a normal distribution, so same as normal on the
log sample.

### 2.3. Beta distribution (of the first kind)

The density function for a beta $\mathcal{B}e(a,b)$ is
$$f_{X}(x) = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a + b)}x^{a - 1}(1 - x)^{b - 1}.$$
The initial estimate is the MME
$$\widehat{a} = m_{1}\delta,\widehat{b} = \left( 1 - m_{1} \right)\delta,\delta = \frac{m_{1}\left( 1 - m_{1} \right)}{\mu_{2}} - 1,(\# eq:betaguessestimator)$$

### 2.4. Other continuous distribution in `actuar`

#### 2.4.1. Log-gamma

Use the gamma initial values on the sample $\log(x)$

#### 2.4.2. Gumbel

The distribution function is
$$F(x) = \exp\left( - \exp\left( - \frac{x - \alpha}{\theta} \right) \right).$$
Let $q_{1}$ and $q_{3}$ the first and the third quartiles. \$\$
\left\\\begin{array} -\theta\log(-\log(p_1)) = q_1-\alpha \\
-\theta\log(-\log(p_3)) = q_3-\alpha \end{array}\right. \Leftrightarrow
\left\\\begin{array} -\theta\log(-\log(p_1))+\theta\log(-\log(p_3)) =
q_1-q_3 \\ \alpha= \theta\log(-\log(p_3)) + q_3 \end{array}\right.
\Leftrightarrow \left\\\begin{array} \theta=
\frac{q_1-q_3}{\log(-\log(p_3)) - \log(-\log(p_1))} \\ \alpha=
\theta\log(-\log(p_3)) + q_3 \end{array}\right.. \$\$ Using the median
for the location parameter $\alpha$ yields to initial estimate
$$\widehat{\theta} = \frac{q_{1} - q_{3}}{\log\left( \log(4/3) \right) - \log\left( \log(4) \right)},\widehat{\alpha} = \widehat{\theta}\log\left( \log(2) \right) + q_{2}.$$

#### 2.4.3. Inverse Gaussian distribution

The moments of this distribution are
$$E\lbrack X\rbrack = \mu,Var\lbrack X\rbrack = \mu^{3}\phi.$$ Hence the
initial estimate are $\widehat{\mu} = m_{1}$,
$\widehat{\phi} = \mu_{2}/m_{1}^{3}$.

#### 2.4.4. Generalized beta

This is the distribution of $\theta X^{1/\tau}$ when $X$ is beta
distributed $\mathcal{B}e(a,b)$ The moments are
$$E\lbrack X\rbrack = \theta\beta(a + 1/\tau,b)/\beta(a,b) = \theta\frac{\Gamma(a + 1/\tau)}{\Gamma(a)}\frac{\Gamma(a + b)}{\Gamma(a + b + 1/\tau)},$$$$E\left\lbrack X^{2} \right\rbrack = \theta^{2}\frac{\Gamma(a + 2/\tau)}{\Gamma(a)}\frac{\Gamma(a + b)}{\Gamma(a + b + 2/\tau)}.$$
Hence for large value of $\tau$, we have
$$E\left\lbrack X^{2} \right\rbrack/E\lbrack X\rbrack = \theta\frac{\Gamma(a + 2/\tau)}{\Gamma(a + b + 2/\tau)}\frac{\Gamma(a + b + 1/\tau)}{\Gamma(a + 1/\tau)} \approx \theta.$$
Note that the MLE of $\theta$ is the maximum We use
$$\widehat{\tau} = 3,\widehat{\theta} = \frac{m_{2}}{m_{1}}\max\limits_{i}x_{i}1_{m_{2} > m_{1}} + \frac{m_{1}}{m_{2}}\max\limits_{i}x_{i}1_{m_{2} \geq m_{1}}.$$
then we use beta initial estimate on sample
$\left( \frac{x_{i}}{\widehat{\theta}} \right)^{\widehat{\tau}}$.

### 2.5. Feller-Pareto family

The Feller-Pareto distribution is the distribution
$X = \mu + \theta(1/B - 1)^{1/\gamma}$ when $B$ follows a beta
distribution with shape parameters $\alpha$ and $\tau$. See details at
<https://doi.org/10.18637/jss.v103.i06> Hence let
$Y = (X - \mu)/\theta$, we have
$$\frac{Y}{1 + Y} = \frac{X - \mu}{\theta + X - \mu} = (1 - B)^{1/\gamma}.$$
For $\gamma$ close to 1, $\frac{Y}{1 + Y}$ is approximately beta
distributed $\tau$ and $\alpha$.

The log-likelihood is
$$\mathcal{L}(\mu,\theta,\alpha,\gamma,\tau) = (\tau\gamma - 1)\sum\limits_{i}\log\left( \frac{x_{i} - \mu}{\theta} \right) - (\alpha + \tau)\sum\limits_{i}\log\left( 1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{\gamma} \right) + n\log(\gamma) - n\log(\theta) - n\log\left( \beta(\alpha,\tau) \right).(\# eq:fellerparetologlik).$$
The MLE of $\mu$ is the minimum.

The gradient with respect to $\theta,\alpha,\gamma,\tau$ is
$$\nabla\mathcal{L}(\mu,\theta,\alpha,\gamma,\tau) = \begin{pmatrix}
{- (\tau\gamma - 1)\sum\limits_{i}\frac{x_{i}}{\theta\left( x_{i} - \mu \right)} + (\alpha + \tau)\sum\limits_{i}\frac{x_{i}\gamma\left( \frac{x_{i} - \mu}{\theta} \right)^{\gamma - 1}}{\theta^{2}\left( 1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{\gamma} \right)} - n/\theta} \\
{- \sum\limits_{i}\log\left( 1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{\gamma} \right) - n\left( \psi(\tau) - \psi(\alpha + \tau) \right)} \\
{(\tau - 1)\sum\limits_{i}\log\left( \frac{x_{i} - \mu}{\theta} \right) - (\alpha + \tau)\sum\limits_{i}\frac{\left( \frac{x_{i} - \mu}{\theta} \right)^{\gamma}}{1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{\gamma}}\log\left( \frac{x_{i} - \mu}{\theta} \right) + n/\gamma} \\
{(\gamma - 1)\sum\limits_{i}\log\left( \frac{x_{i} - \mu}{\theta} \right) - \sum\limits_{i}\log\left( 1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{\gamma} \right) - n\left( \psi(\tau) - \psi(\alpha + \tau) \right)}
\end{pmatrix}.(\# eq:fellerparetogradient)$$ Cancelling the first
component of score for $\gamma = \alpha = 2$, we get
$$\left. - (2\tau - 1)\sum\limits_{i}\frac{x_{i}}{\theta\left( x_{i} - \mu \right)} + (2 + \tau)\sum\limits_{i}\frac{x_{i}2\left( x_{i} - \mu \right)}{\theta^{3}\left( 1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{2} \right)} = \frac{n}{\theta}\Leftrightarrow - (2\tau - 1)\theta^{2}\frac{1}{n}\sum\limits_{i}\frac{x_{i}}{x_{i} - \mu} + (2 + \tau)\frac{1}{n}\sum\limits_{i}\frac{x_{i}2\left( x_{i} - \mu \right)}{\left( 1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{2} \right)} = \theta^{2} \right.$$$$\left. \Leftrightarrow(2 + \tau)\frac{1}{n}\sum\limits_{i}\frac{x_{i}2\left( x_{i} - \mu \right)}{1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{2}} = (2\tau - 1)\theta^{2}\left( \frac{1}{n}\sum\limits_{i}\frac{x_{i}}{x_{i} - \mu} - 1 \right)\Leftrightarrow\sqrt{\frac{(2 + \tau)\frac{1}{n}\sum\limits_{i}\frac{x_{i}2\left( x_{i} - \mu \right)}{1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{2}}}{(2\tau - 1)\left( \frac{1}{n}\sum\limits_{i}\frac{x_{i}}{x_{i} - \mu} - 1 \right)}} = \theta. \right.$$
Neglecting unknown value of $\tau$ and the denominator in $\theta$, we
get with $\widehat{\mu}$ set with (@ref(eq:pareto4muinit))
$$\widehat{\theta} = \sqrt{\frac{\frac{1}{n}\sum\limits_{i}\frac{x_{i}2\left( x_{i} - \widehat{\mu} \right)}{1 + \left( x_{i} - \widehat{\mu} \right)^{2}}}{\left( \frac{1}{n}\sum\limits_{i}\frac{x_{i}}{x_{i} - \widehat{\mu}} - 1 \right)}}.(\# eq:fellerparetothetahat)$$
Initial value of $\tau,\alpha$ are obtained on the sample
$\left( z_{i} \right)_{i}$$$z_{i} = y_{i}/\left( 1 + y_{i} \right),y_{i} = \left( x_{i} - \widehat{\mu} \right)/\widehat{\theta},$$
with initial values of a beta distribution which is based on MME
(@ref(eq:betaguessestimator)).

Cancelling the last component of the gradient leads to
$$\left. (\gamma - 1)\frac{1}{n}\sum\limits_{i}\log\left( \frac{x_{i} - \mu}{\theta} \right) - \frac{1}{n}\sum\limits_{i}\log\left( 1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{\gamma} \right) = \psi(\tau) - \psi(\alpha + \tau)\Leftrightarrow(\gamma - 1)\frac{1}{n}\sum\limits_{i}\log\left( \frac{x_{i} - \mu}{\theta} \right) = \psi(\tau) - \psi(\alpha + \tau) + \frac{1}{n}\sum\limits_{i}\log\left( 1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{\gamma} \right). \right.$$
Neglecting the value $\gamma$ on the right-hand side we obtain
$$\widehat{\gamma} = 1 + \frac{\psi(\tau) - \psi(\alpha + \tau) + \frac{1}{n}\sum\limits_{i}\log\left( 1 + \left( \frac{x_{i} - \mu}{\theta} \right) \right)}{\frac{1}{n}\sum\limits_{i}\log\left( \frac{x_{i} - \mu}{\theta} \right)}.(\# eq:fellerparetogammahat)$$

#### 2.5.1. Transformed beta

This is the Feller-Pareto with $\mu = 0$. So the first component of
@ref(eq:fellerparetogradient) simplifies to with
$\gamma = \alpha = 2$$$\left. - (2\tau - 1)\sum\limits_{i}\frac{x_{i}}{\theta\left( x_{i} \right)} + (2 + \tau)\sum\limits_{i}\frac{2x_{i}^{2}}{\theta^{3}\left( 1 + \left( \frac{x_{i}}{\theta} \right)^{2} \right)} = \frac{n}{\theta}\Leftrightarrow - (2\tau - 1)\theta^{2} + (2 + \tau)\frac{1}{n}\sum\limits_{i}\frac{2x_{i}^{2}}{1 + \left( \frac{x_{i}}{\theta} \right)^{2}} = \theta^{2} \right.$$$$\theta^{2} = \frac{2 + \tau}{2\tau}\frac{1}{n}\sum\limits_{i}\frac{2x_{i}^{2}}{1 + \left( \frac{x_{i}}{\theta} \right)^{2}}.$$
Neglecting unknown value of $\tau$ in the denominator in $\theta$, we
get
$$\widehat{\theta} = \sqrt{\frac{1}{n}\sum\limits_{i}\frac{2x_{i}^{2}}{1 + x_{i}^{2}}}.(\# eq:trbetathetahat)$$
Initial value of $\tau,\alpha$ are obtained on the sample
$\left( z_{i} \right)_{i}$$$z_{i} = y_{i}/\left( 1 + y_{i} \right),y_{i} = x_{i}/\widehat{\theta},$$
with initial values of a beta distribution which is based on MME
(@ref(eq:betaguessestimator)). Similar to Feller-Pareto, we set
$$\widehat{\gamma} = 1 + \frac{\psi(\tau) - \psi(\alpha + \tau) + \frac{1}{n}\sum\limits_{i}\log\left( 1 + \frac{x_{i}}{\theta} \right)}{\frac{1}{n}\sum\limits_{i}\log\left( \frac{x_{i}}{\theta} \right)}.(\# eq:fellerparetogammahat)$$

#### 2.5.2. Generalized Pareto

This is the Feller-Pareto with $\mu = 0$$\gamma = 1$. So the first
component of @ref(eq:fellerparetogradient) simplifies to with
$\gamma = 2$$$\left. - (\tau - 1)\frac{n}{\theta} + (2 + \tau)\sum\limits_{i}\frac{x_{i}}{\theta^{2}(1 + \frac{x_{i}}{\theta}} = n/\theta\Leftrightarrow - (\tau - 1)\theta + (2 + \tau)\frac{1}{n}\sum\limits_{i}\frac{x_{i}}{(1 + \frac{x_{i}}{\theta}} = \theta. \right.$$
Neglecting unknown value of $\tau$ leads to
$$\widehat{\theta} = \frac{1}{n}\sum\limits_{i}\frac{x_{i}}{1 + x_{i}}(\# eq:generalizedparetotheta)$$

Initial value of $\tau,\alpha$ are obtained on the sample
$\left( z_{i} \right)_{i}$$$z_{i} = y_{i}/\left( 1 + y_{i} \right),y_{i} = x_{i}/\widehat{\theta},$$
with initial values of a beta distribution which is based on MME
(@ref(eq:betaguessestimator)).

#### 2.5.3. Burr

Burr is a Feller-Pareto distribution with $\mu = 0$, $\tau = 1$.

The survival function is
$$1 - F(x) = \left( 1 + (x/\theta)^{\gamma} \right)^{- \alpha}.$$ Using
the median $q_{2}$, we have
$$\log(1/2) = - \alpha\log\left( 1 + \left( q_{2}/\theta \right)^{\gamma} \right).$$
The initial value is
$$\alpha = \frac{\log(2)}{\log\left( 1 + \left( q_{2}/\theta \right)^{\gamma} \right)},(\# eq:burralpharelation)$$

So the first component of @ref(eq:fellerparetogradient) simplifies to
with $\gamma = \alpha = 2$, $\tau = 1$, $\mu = 0$.
$$\left. - n/\theta + 3\sum\limits_{i}\frac{2x_{i}\left( \frac{x_{i}}{\theta} \right)}{\theta^{2}\left( 1 + \left( \frac{x_{i}}{\theta} \right)^{2} \right)} = n/\theta\Leftrightarrow\theta^{2}\frac{1}{n}\sum\limits_{i}\frac{2x_{i}\left( \frac{x_{i}}{\theta} \right)}{\left( 1 + \left( \frac{x_{i}}{\theta} \right)^{2} \right)} = 2/3. \right.$$
Neglecting unknown value in the denominator in $\theta$, we get
$$\widehat{\theta} = \sqrt{\frac{2}{3\frac{1}{n}\sum\limits_{i}\frac{2x_{i}^{2}}{1 + \left( x_{i} \right)^{2}}}}.(\# eq:trbetathetahat)$$
We use for $\widehat{\gamma}$ @ref(eq:fellerparetogammahat) with
$\tau = 1$ and $\alpha = 2$ and previous $\widehat{\theta}$.

#### 2.5.4. Loglogistic

Loglogistic is a Feller-Pareto distribution with $\mu = 0$, $\tau = 1$,
$\alpha = 1$. The survival function is
$$1 - F(x) = \left( 1 + (x/\theta)^{\gamma} \right)^{- 1}.$$ So
$$\left. \frac{1}{1 - F(x)} - 1 = (x/\theta)^{\gamma}\Leftrightarrow\log\left( \frac{F(x)}{1 - F(x)} \right) = \gamma\log(x/\theta). \right.$$
Let $q_{1}$ and $q_{3}$ be the first and the third quartile.
$$\left. \log\left( \frac{1/3}{2/3} \right) = \gamma\log\left( q_{1}/\theta \right),\log\left( \frac{2/3}{1/3} \right) = \gamma\log\left( q_{3}/\theta \right)\Leftrightarrow - \log(2) = \gamma\log\left( q_{1}/\theta \right),\log(2) = \gamma\log\left( q_{3}/\theta \right). \right.$$
The difference of previous equations simplifies to
$$\widehat{\gamma} = \frac{2\log(2)}{\log\left( q_{3}/q_{1} \right)}.$$
The sum of previous equations
$$0 = \gamma\log\left( q_{1} \right) + \gamma\log\left( q_{3} \right) - 2\gamma\log(\theta).$$$$\widehat{\theta} = \frac{1}{2}e^{\log{(q_{1}q_{3})}}.(\# eq:llogisthetahat)$$

#### 2.5.5. Paralogistic

Paralogistic is a Feller-Pareto distribution with $\mu = 0$, $\tau = 1$,
$\alpha = \gamma$. The survival function is
$$1 - F(x) = \left( 1 + (x/\theta)^{\alpha} \right)^{- \alpha}.$$ So
$$\log\left( 1 - F(x) \right) = - \alpha\log\left( 1 + (x/\theta)^{\alpha} \right).$$
The log-likelihood is
$$\mathcal{L}(\theta,\alpha) = (\alpha - 1)\sum\limits_{i}\log\left( \frac{x_{i}}{\theta} \right) - (\alpha + 1)\sum\limits_{i}\log\left( 1 + \left( \frac{x_{i}}{\theta} \right)^{\alpha} \right) + 2n\log(\alpha) - n\log(\theta).(\# eq:paralogisloglik)$$
The gradient with respect to $\theta$, $\alpha$ is $$\begin{pmatrix}
{(\alpha - 1)\frac{- n}{\theta} - (\alpha + 1)\sum\limits_{i}\frac{- x_{i}\alpha\left( x_{i}/\theta \right)^{\alpha - 1}}{1 + \left( \frac{x_{i}}{\theta} \right)^{\alpha}} - n/\theta} \\
{\sum\limits_{i}\log\left( \frac{\frac{x_{i}}{\theta}}{1 + \left( \frac{x_{i}}{\theta} \right)^{\alpha}} \right) - (\alpha + 1)\sum\limits_{i}\frac{\left( \frac{x_{i}}{\theta} \right)^{\alpha}\log\left( x_{i}/\theta \right)}{1 + \left( \frac{x_{i}}{\theta} \right)^{\alpha}} + 2n/\alpha} \\

\end{pmatrix}.$$ The first component cancels when
$$\left. - (\alpha + 1)\sum\limits_{i}\frac{- x_{i}\alpha\left( x_{i}/\theta \right)^{\alpha - 1}}{1 + \left( \frac{x_{i}}{\theta} \right)^{\alpha}} = \alpha n/\theta\Leftrightarrow(\alpha + 1)\frac{1}{n}\sum\limits_{i}\frac{\left( x_{i} \right)^{\alpha + 1}}{1 + \left( \frac{x_{i}}{\theta} \right)^{\alpha}} = \theta^{\alpha}. \right.$$
The second component cancels when
$$\frac{1}{n}\sum\limits_{i}\log\left( \frac{\frac{x_{i}}{\theta}}{1 + \left( \frac{x_{i}}{\theta} \right)^{\alpha}} \right) = - 2/\alpha + (\alpha + 1)\frac{1}{n}\sum\limits_{i}\frac{\left( \frac{x_{i}}{\theta} \right)^{\alpha}\log\left( x_{i}/\theta \right)}{1 + \left( \frac{x_{i}}{\theta} \right)^{\alpha}}.$$
Choosing $\theta = 1$, $\alpha = 2$ in sums leads to
$$\frac{1}{n}\sum\limits_{i}\log\left( \frac{\frac{x_{i}}{\theta}}{1 + x_{i}^{2}} \right) - \frac{1}{n}\sum\limits_{i}\frac{x_{i}^{2}\log\left( x_{i} \right)}{1 + x_{i}^{2}} = - 2/\alpha + (\alpha)\frac{1}{n}\sum\limits_{i}\frac{x_{i}^{2}\log\left( x_{i} \right)}{1 + x_{i}^{2}}.$$
Initial estimators are
$$\widehat{\alpha} = \frac{\frac{1}{n}\sum\limits_{i}\log\left( \frac{x_{i}}{1 + x_{i}^{2}} \right) - \frac{1}{n}\sum\limits_{i}\frac{x_{i}^{2}\log\left( x_{i} \right)}{1 + x_{i}^{2}}}{\frac{1}{n}\sum\limits_{i}\frac{x_{i}^{2}\log\left( x_{i} \right)}{1 + x_{i}^{2}} - 2},(\# eq:paralogisalphahat)$$$$\widehat{\theta} = \left( \widehat{\alpha} + 1 \right)\frac{1}{n}\sum\limits_{i}\frac{\left( x_{i} \right)^{\widehat{\alpha} + 1}}{1 + \left( x_{i} \right)^{\widehat{\alpha}}}.(\# eq:paralogisthetahat)$$

#### 2.5.6. Inverse Burr

Use Burr estimate on the sample $1/x$ then inverse the scale parameter.

#### 2.5.7. Inverse paralogistic

Use paralogistic estimate on the sample $1/x$ then inverse the scale
parameter.

#### 2.5.8. Inverse pareto

Use pareto estimate on the sample $1/x$ then inverse the scale
parameter.

#### 2.5.9. Pareto IV

The survival function is
$$1 - F(x) = \left( 1 + \left( \frac{x - \mu}{\theta} \right)^{\gamma} \right)^{- \alpha},$$
see [`?Pareto4`](https://rdrr.io/pkg/actuar/man/Pareto4.html) in
`actuar`.

The first and third quartiles $q_{1}$ and $q_{3}$ verify
$$\left( \left( \frac{3}{4} \right)^{- 1/\alpha} - 1 \right)^{1/\gamma} = \frac{q_{1} - \mu}{\theta},\left( \left( \frac{1}{4} \right)^{- 1/\alpha} - 1 \right)^{1/\gamma} = \frac{q_{3} - \mu}{\theta}.$$
Hence we get two useful relations
$$\gamma = \frac{\log\left( \frac{\left( \frac{4}{3} \right)^{1/\alpha} - 1}{(4)^{1/\alpha} - 1} \right)}{\log\left( \frac{q_{1} - \mu}{q_{3} - \mu} \right)},(\# eq:pareto4gammarelation)$$$$\theta = \frac{q_{1} - q_{3}}{\left( \left( \frac{4}{3} \right)^{1/\alpha} - 1 \right)^{1/\gamma} - \left( (4)^{1/\alpha} - 1 \right)^{1/\gamma}}.(\# eq:pareto4thetarelation)$$

The log-likelihood of a Pareto 4 sample (see Equation (5.2.94) of Arnold
(2015) updated with Goulet et al. notation) is
$$\mathcal{L}(\mu,\theta,\gamma,\alpha) = (\gamma - 1)\sum\limits_{i}\log\left( \frac{x_{i} - \mu}{\theta} \right) - (\alpha + 1)\sum\limits_{i}\log\left( 1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{\gamma} \right) + n\log(\gamma) - n\log(\theta) + n\log(\alpha).$$
Cancelling the derivate of $\mathcal{L}(\mu,\theta,\gamma,\alpha)$ with
respect to $\alpha$ leads to
$$\alpha = n/\sum\limits_{i}\log\left( 1 + \left( \frac{x_{i} - \mu}{\theta} \right)^{\gamma} \right).(\# eq:pareto4alpharelation)$$

The MLE of the threshold parameter $\mu$ is the minimum. So the initial
estimate is slightly under the minimum in order that all observations
are strictly above it $$\widehat{\mu} = \begin{cases}
{(1 - \epsilon)\min\limits_{i}x_{i}} & {{\text{if}\mspace{6mu}}\min\limits_{i}x_{i} < 0} \\
{(1 + \epsilon)\min\limits_{i}x_{i}} & {{\text{if}\mspace{6mu}}\min\limits_{i}x_{i} \geq 0} \\
 & 
\end{cases}.(\# eq:pareto4muinit)$$ where $\epsilon = 0.05$.

Initial parameter estimation is $\widehat{\mu}$, $\alpha^{\star} = 2$ ,
$\widehat{\gamma}$ from @ref(eq:pareto4gammarelation) with
$\alpha^{\star}$, $\widehat{\theta}$ from @ref(eq:pareto4thetarelation)
with $\alpha^{\star}$ and $\widehat{\gamma}$, $\widehat{\alpha}$ from
@ref(eq:pareto4alpharelation) with $\widehat{\mu}$, $\widehat{\theta}$
and $\widehat{\gamma}$.

#### 2.5.10. Pareto III

Pareto III corresponds to Pareto IV with $\alpha = 1$.
$$\gamma = \frac{\log\left( \frac{\frac{4}{3} - 1}{4 - 1} \right)}{\log\left( \frac{q_{1} - \mu}{q_{3} - \mu} \right)},$$

$$\theta = \frac{\left( \frac{1}{3} \right)^{1/\gamma} - (3)^{1/\gamma}}{q_{1} - q_{3}}.$$

Initial parameter estimation is $\widehat{\mu}$, $\widehat{\gamma}$ from
, $\widehat{\theta}$ from with $\widehat{\gamma}$.

#### 2.5.11. Pareto II

Pareto II corresponds to Pareto IV with $\gamma = 1$.

$$\theta = \frac{\left( \frac{4}{3} \right)^{1/\alpha} - 4^{1/\alpha}}{q_{1} - q_{3}}.$$

Initial parameter estimation is $\widehat{\mu}$, $\alpha^{\star} = 2$ ,
$\widehat{\theta}$ from with $\alpha^{\star}$ and $\gamma = 1$,
$\widehat{\alpha}$ from with $\widehat{\mu}$, $\widehat{\theta}$ and
$\gamma = 1$,

#### 2.5.12. Pareto I

Pareto I corresponds to Pareto IV with $\gamma = 1$, $\mu = \theta$.

The MLE is
$$\widehat{\mu} = \min\limits_{i}X_{i},\widehat{\alpha} = \left( \frac{1}{n}\sum\limits_{i = 1}^{n}\log\left( X_{i}/\widehat{\mu} \right) \right)^{- 1}.$$

This can be rewritten with the geometric mean of the sample
$G_{n} = \left( \prod_{i = 1}^{n}X_{i} \right)^{1/n}$ as
$$\widehat{\alpha} = \log\left( G_{n}/\widehat{\mu} \right).$$

Initial parameter estimation is $\widehat{\mu}$, $\widehat{\alpha}$ from
.

#### 2.5.13. Pareto

Pareto corresponds to Pareto IV with $\gamma = 1$, $\mu = 0$.
$$\theta = \frac{\left( \frac{4}{3} \right)^{1/\alpha} - 4^{1/\alpha}}{q_{1} - q_{3}}.$$

Initial parameter estimation is
$$\alpha^{\star} = \max\left( 2,2\left( m_{2} - m_{1}^{2} \right)/\left( m_{2} - 2m_{1}^{2} \right) \right),$$
with $m_{i}$ are empirical raw moment of order $i$, $\widehat{\theta}$
from with $\alpha^{\star}$ and $\gamma = 1$, $\widehat{\alpha}$ from
with $\mu = 0$, $\widehat{\theta}$ and $\gamma = 1$.

### 2.6. Transformed gamma family

#### 2.6.1. Transformed gamma distribution

The log-likelihood is given by
$$\mathcal{L}(\alpha,\tau,\theta) = n\log(\tau) + \alpha\tau\sum\limits_{i}\log\left( x_{i}/\theta \right) - \sum\limits_{i}\left( x_{i}/\theta \right)^{\tau} - \sum\limits_{i}\log\left( x_{i} \right) - n\log\left( Gamma(\alpha) \right).$$
The gradient with respect to $\alpha,\tau,\theta$ is given by
$$\begin{pmatrix}
{\tau - n\psi(\alpha))} \\
{n/\tau + \alpha\sum\limits_{i}\log\left( x_{i}/\theta \right) - \sum\limits_{i}\left( x_{i}/\theta \right)^{\tau}\log\left( x_{i}/\theta \right)} \\
{- \alpha\tau/\theta + \sum\limits_{i}\tau\frac{x_{i}}{\theta^{2}}\left( x_{i}/\theta \right)^{\tau - 1}}
\end{pmatrix}.$$ We compute the moment-estimator as in gamma
$$\widehat{\alpha} = m_{2}^{2}/\mu_{2},\widehat{\theta} = \mu_{2}/m_{1}.$$
Then cancelling the first component of the gradient we set
$$\widehat{\tau} = \frac{\psi\left( \widehat{\alpha} \right)}{\frac{1}{n}\sum\limits_{i}\log\left( x_{i}/\widehat{\theta} \right)}.$$

#### 2.6.2. gamma distribution

Transformed gamma with $\tau = 1$

We compute the moment-estimator given by
$$\widehat{\alpha} = m_{2}^{2}/\mu_{2},\widehat{\theta} = \mu_{2}/m_{1}.$$

#### 2.6.3. Weibull distribution

Transformed gamma with $\alpha = 1$

Let $\widetilde{m} = \frac{1}{n}\sum_{i}\log\left( x_{i} \right)$ and
$\widetilde{v} = \frac{1}{n}\sum_{i}\left( \log\left( x_{i} \right) - \widetilde{m} \right)^{2}$.
We use an approximate MME
$$\widehat{\tau} = 1.2/sqrt\left( \widetilde{v} \right),\widehat{\theta} = exp\left( \widetilde{m} + 0.572/\widehat{\tau} \right).$$
Alternatively, we can use the distribution function
$$\left. F(x) = 1 - e^{- {(x/\sigma)}^{\tau}}\Rightarrow\log\left( - \log\left( 1 - F(x) \right) \right) = \tau\log(x) - \tau\log(\theta), \right.$$
Hence the QME for Weibull is
$$\widetilde{\tau} = \frac{\log\left( - \log\left( 1 - p_{1} \right) \right) - \log\left( - \log\left( 1 - p_{2} \right) \right)}{\log\left( x_{1} \right) - \log\left( x_{2} \right)},\widetilde{\tau} = x_{3}/\left( - \log\left( 1 - p_{3} \right) \right)^{1/\widetilde{\tau}}$$
with $p_{1} = 1/4$, $p_{2} = 3/4$, $p_{3} = 1/2$, $x_{i}$ corresponding
empirical quantiles.

Initial parameters are $\widetilde{\tau}$ and $\widetilde{\theta}$
unless the empirical quantiles $x_{1} = x_{2}$, in that case we use
$\widehat{\tau}$, $\widehat{\theta}$.

#### 2.6.4. Exponential distribution

The MLE is the MME $\widehat{\lambda} = 1/m_{1}.$

### 2.7. Inverse transformed gamma family

#### 2.7.1. Inverse transformed gamma distribution

Same as transformed gamma distribution with $\left( 1/x_{i} \right)_{i}$
then inverse the scale parameter.

#### 2.7.2. Inverse gamma distribution

We compute moment-estimator as
$$\widehat{\alpha} = \left( 2m_{2} - m_{1}^{2} \right)/\left( m_{2} - m_{1}^{2} \right),\widehat{\theta} = m_{1}m_{2}/\left( m_{2} - m_{1}^{2} \right).$$

#### 2.7.3. Inverse Weibull distribution

We use the QME.

#### 2.7.4. Inverse exponential

Same as transformed gamma distribution with $\left( 1/x_{i} \right)_{i}$
then inverse the rate parameter.

## 3. Bibliography

### 3.1. General books

- N. L. Johnson, S. Kotz, N. Balakrishnan (1994). Continuous univariate
  distributions, Volume 1, Wiley.
- N. L. Johnson, S. Kotz, N. Balakrishnan (1995). Continuous univariate
  distributions, Volume 2, Wiley.
- N. L. Johnson, A. W. Kemp, S. Kotz (2008). Univariate discrete
  distributions, Wiley.
- G. Wimmer (1999), Thesaurus of univariate discrete probability
  distributions.

### 3.2. Books dedicated to a distribution family

- M. Ahsanullah, B.M. Golam Kibria, M. Shakil (2014). Normal and
  Student’s t Distributions and Their Applications, Springer.
- B. C. Arnold (2010). Pareto Distributions, Chapman and Hall.
- A. Azzalini (2013). The Skew-Normal and Related Families.
- N. Balakrishnan (2014). Handbook of the Logistic Distribution, CRC
  Press.

### 3.3. Books with applications

- C. Forbes, M. Evans, N. Hastings, B. Peacock (2011). Statistical
  Distributions, Wiley.
- Z. A. Karian, E. J. Dudewicz, K. Shimizu (2010). Handbook of Fitting
  Statistical Distributions with R, CRC Press.
- K. Krishnamoorthy (2015). Handbook of Statistical Distributions with
  Applications, Chapman and Hall.
- Klugman, S., Panjer, H. & Willmot, G. (2019). Loss Models: From Data
  to Decisions, 5th ed., John Wiley & Sons.
