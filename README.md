# Estimation of Ornstein-Uhlenbeck SDE Solutions via Monte Carlo Simulation and Variance Reduction Analysis

![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python)
![NumPy](https://img.shields.io/badge/NumPy-Vectorized-013243?logo=numpy)
![SciPy](https://img.shields.io/badge/SciPy-Optimization-8CAAE6?logo=scipy)
![License](https://img.shields.io/badge/License-MIT-green.svg)

Academic project completed for the course **Numerical Methods for Stochastic Differential Equations (MNSRR)**.

**Authors:**
* Julia Mierzyńska
* Jolanta Sułkowska

---

## 📌 Project Overview

The primary objective of this project is the numerical approximation of the stochastic differential equation (SDE) governing the one-factor **Hull-White interest rate model** (an extended **Ornstein-Uhlenbeck process**) using **Monte Carlo simulation**. 

The core tasks include:
1. Pricing a European Call option written on a zero-coupon bond.
2. Calibrating model parameters against market quotes using numerical optimization.
3. Implementing and evaluating the computational performance of **variance reduction techniques**:
   * **Antithetic Variates**
   * **Control Variates**
4. Comparing simulation results against the exact analytical benchmark provided by the **Black-76 bond option formula**.

---

## 📐 Mathematical Formulation

### 1. The Stochastic Process
Under the forward measure, with a long-term mean reversion level $\mu = 0$ and speed of mean reversion $a$, the short-rate process satisfies:

$$dX_t = -a X_t \, dt + \sigma(t) \, dW_t$$

where:
* $a > 0$ is the mean-reversion rate,
* $\sigma(t)$ is a piecewise-constant volatility function,
* $W_t$ is a standard Brownian motion.

The exact analytical integral solution of this Gaussian process is:

$$X(t) = X(0) e^{-at} + \int_0^t \sigma(s) e^{-a(t-s)} \, dW(s)$$

### 2. Zero-Coupon Bond Option Payoff
Let $T$ denote option maturity and $S$ ($S > T$) the maturity of the underlying zero-coupon bond. The price of the underlying bond at time $T$ is given by:

$$P(T,S) = A(T,S) \exp\left(-B(T,S) X_T\right)$$

where:
* $B(T,S) = \frac{1 - e^{-a(S-T)}}{a}$
* $A(T,S)$ aligns the model with the initial term structure of zero-coupon bond prices $P(0, T)$ and $P(0, S)$.

The discounted payoff for a European Call option with strike $K$ is:

$$H = P(0,T) \max\left(P(T,S) - K, 0\right)$$

---

## 🔬 Implementation Pipeline

### Step 1: Analytical Pricing (Benchmark)
* Implemented the closed-form Black formula for European bond options.
* Integrated the piecewise-constant volatility over time intervals to compute forward variance $v(t)$.
* Established an exact ground truth to isolate and measure statistical and discretization errors.

### Step 2: Parameter Calibration
* Constructed an objective function measuring the Sum of Squared Errors (SSE) between model-implied and market prices across strike levels $K_i$:

$$\min_{a} \sum_{i=1}^L \left( Call_{\text{model}}(K_i, a) - Call_{\text{market}}(K_i) \right)^2$$

* Solved using bounded scalar optimization (`scipy.optimize.minimize_scalar`).

### Step 3: Vectorized Euler-Maruyama Discretization
* Discretized the continuous time interval $[0, T]$ into $M$ sub-intervals of length $h = \frac{T}{M}$:

$$Y_n = Y_{n-1} - a Y_{n-1} h + \sigma(t_{n-1}) \sqrt{h} Z_n, \quad Z_n \sim \mathcal{N}(0, 1)$$

* Fully vectorized across paths using NumPy matrices ($M \times N$), eliminating Python-level iteration over Monte Carlo trajectories.

### Step 4: Discretization Error Analysis
* Evaluated systematic time-step bias by sweeping $M \in [10, 200]$ at a fixed high sample size ($N = 100\,000$).
* Observed the transition where deterministic discretization error becomes subdominant relative to Monte Carlo sampling noise.

### Step 5: Variance Reduction Framework
1. **Antithetic Variates:** Exploits the reflective symmetry of the standard Gaussian distribution by evaluating trajectories driven by $+Z$ and $-Z$ simultaneously:

$$\hat{C}_{\text{AV}} = \frac{1}{N/2} \sum_{i=1}^{N/2} \frac{H_i(Z) + H_i(-Z)}{2}$$

2. **Control Variates** 

### Step 6: Convergence and Error Scaling
* Assessed empirical convergence on a log-log scale across $N \in [10^2, 10^5]$ against the theoretical Central Limit Theorem rate $\mathcal{O}(1/\sqrt{N})$.

---

## 📊 Results & Key Findings

Test configuration: $T = 688/360$, $S = 900/360$, $K = 0.77$, $a = 0.04$, $P(0,T) = 0.80$, $P(0,S) = 0.77$.
* **Exact Analytical Price (Black-76):** `0.159489`

| Method | Estimated Price | Estimator Variance ($\text{Var}/N$) | Variance Reduction Ratio |
| :--- | :---: | :---: | :---: |
| **Standard Monte Carlo** | 0.165284 | $1.98 \times 10^{-6}$ | **1.00x** (Baseline) |
| **Antithetic Variates** | 0.164821 | $1.43 \times 10^{-7}$ | **~13.8x** |
| **Control Variates** | 0.164645 | $1.49 \times 10^{-7}$ | **~13.3x** |

### Numerical Insights:
1. **Antithetic Performance:** Yielded an approximate **14-fold variance reduction** with negligible CPU overhead because negative increments $-Z$ require no additional pseudorandom draws. The linear drift structure ($-a Y$) preserves strong negative correlation across paths.
2. **Impact of Payoff Non-Linearity:** While control variates are theoretically superior in linear settings, the compound non-linearity of the bond call payoff ($P(T,S) \propto e^{-B X_T}$ followed by $\max(P - K, 0)$) dampens linear correlation with $X_T$, bringing its efficiency close to the antithetic method.
3. **Discretization Bias vs. Statistical Error:** All Monte Carlo estimates plateau around `0.1646 - 0.1652`, diverging slightly from the Black formula benchmark (`0.159489`). This highlights that variance reduction techniques address estimator variance, not the systematic bias introduced by finite-step Euler-Maruyama approximation. An exact simulation scheme based on conditional Gaussian transitions would eliminate this remaining bias.
