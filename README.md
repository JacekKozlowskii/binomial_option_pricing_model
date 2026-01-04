<a id="readme-top"></a>

# Binomial Option Pricing Model

This project implements a Binomial Tree model for pricing options in Python. It supports both European and American options (calls and puts).

## Table of Contents:
1. [Setup](#setup)
2. [Code Overview](#code-overview)
3. [Optimization](#optimization)
   1. [Node Traversal](#1-node-traversal-the-move-variable)
   2. [Calculating the k-index](#2-calculating-the-k-index)
   3. [American Put Split](#3-american-put-split)
   4. [European and American Call Option Calculations](#4-european-and-american-call-option-calculations)
4. [Disclaimer](#disclaimer)

## Setup

---

### Prerequisites
- Python 3.10+: _recommended_, the code uses union type hints (e.g. ```int | float```).

### Running the Code

#### Installation

1. Clone the repository:
```bash
git clone [https://github.com/JacekKozlowskii/binomial_option_pricing_model.git](https://github.com/JacekKozlowskii/binomial_option_pricing_model.git)
cd binomial_option_pricing_model
```
2. Run the script (this will return calculations for an example option):
```bash
python binomial_option_pricing_model.py
```

3. Import the function:

```python
from binomial_option_pricing_model import price_option
```



## Code Overview

---

The core function, `price_option`, calculates the fair value of an option using a binomial tree model.

### Parameters

* **spot**: Current price of the underlying asset.
* **strike**: Strike price of the option.
* **sigma**: Volatility (standard deviation) of the underlying asset.
* **steps**: Number of steps (height) of the binomial tree.
* **T**: Time to maturity (in years).
* **rate**: Risk-free interest rate.
* **side**: `'call'` (or `1`) / `'put'` (or `-1`).
* **kind**: `'european'` (or `0`) / `'american'` (or `1`).
* **relative**: Determines whether time is measured in absolute or relative terms.

**Returns**

The function returns a tuple containing:

1. **Option Price**: The calculated fair value.
2. **Computation Time**: Time taken to execute the calculation.
3. **Paths Considered**: The number of theoretical paths/nodes processed.

## Optimization

---

The code uses a few optimization techniques to avoid a nested-loop binomial tree implementation.
Below are explanations of the specific optimizations noted in the comments.


### 1. Node Traversal (The move variable)

**Code reference:**
```python 
move = d * d if side == 1 else u * u
```

This approach avoids calculating node prices at terminal time using computationally expensive power operations such as
$S_0 \cdot u^{j} \cdot d^{n-j}$.

With this method, prices are calculated from the most in-the-money (ITM) to the most out-of-the-money (OTM). This ordering is also the reason the `move` variable is determined by the option side.

* **Call (side == 1)**:
  We start at the highest price ($S_0 \cdot u^n$) and move downward.
* **Put (side == -1)**:
  We start at the lowest price ($S_0 \cdot d^n$) and move upward.

### 2. Calculating k-index

We derive the index at which the option payoff becomes zero and skip all non-valuable outcomes.

**Call**
$$S_0e^{(n-j)\sigma\sqrt{t}}e^{-j\sigma\sqrt{t}} > K$$
$$e^{n\sigma\sqrt{t}-2j\sigma\sqrt{t}} > \frac{K}{S_0}$$
$$(n-2j)\sigma\sqrt{t}>\ln(\frac{K}{S_0})$$
$$j < \frac{n}{2} - \frac{\ln(\frac{K}{S_0})}{2\sigma\sqrt{t}}$$

**Put**
$$S_0e^{(i)\sigma\sqrt{t}}e^{-(n-i)\sigma\sqrt{t}} < K$$
$$i < \frac{n}{2} - \frac{\ln(\frac{S_0}{K})}{2\sigma\sqrt{t}}$$

We then take the ceiling of these values to obtain an integer index, ensuring that only non-valuable outcomes are skipped.
If this index is less than or equal to zero, the price will never reach the strike price.

### 3. American Put Split

Code reference: 
```python
if kind and side == -1:
    ...
```

We only consider American puts when performing full tree induction because they can be exercised before expiration.

For American calls (assuming no dividends, constant volatility and risk-neutral world) we can derive the following.
We want to find the index at witch the payoff of current node is greater than the expected value of future prices.
Let $S_{i,j}$ denote $S_0u^jd^{i-j}$, where $i$ is the step such that $i\in[0,n-1]$ and $j$ be a node at step $i$ such that $j \in [0, i]$.

$$S_{i,j}-K>e^{-rt}(p(S_{i,j}u-K) + (1-p)(S_{i,j}d-K))$$
$$S_{i,j}-K>e^{-rt}(pS_{i,j}u-pK + (1-p)S_{i,j}d-(1-p)K)$$
$$S_{i,j}-K>e^{-rt}(pS_{i,j}u + (1-p)S_{i,j}d - (1-p+p)K)$$
$$S_{i,j}-K>e^{-rt}(pS_{i,j}u + (1-p)S_{i,j}d - K)$$
$$S_{i,j}-K>e^{-rt}(S_{i,j}(\frac{e^{rt}-d}{u-d}(u-d) + d) - K)$$
$$S_{i,j}-K>e^{-rt}(S_{i,j}(e^{rt}-d + d) - K)$$
$$S_{i,j}-K>e^{-rt}(S_{i,j}e^{rt} - K)$$
$$S_{i,j}-K>e^{-rt}(S_{i,j}e^{rt} - K)$$
$$S_{i,j} - K > S_{i,j} - Ke^{-rt}$$
$$K < Ke^{-rt}$$
$$1 \not< e^{-rt}$$

This inequality does not hold, meaning that for any state, the intrinsic value is strictly less than the discounted expected value of the option at the next step. This mean that we can show that starting from $i=n-1$,
the expected value of next payoffs is greater than current step payoff and as that we have following inequality:
$$S_{i-1,j}-K < e^{-rt}(p(S_{i-1,j}u-K) + (1-p)(S_{i-1,j}d-K)) \leq f_{i+1}$$

Thus we never exercise the American call and calculate the price the same way as the European call.

For **American puts** we can show the following:
$$K - S_{i,j}>e^{-rt}(p(K - S_{i,j}u) + (1-p)(K - S_{i,j}d))$$
That simplifies to:
$$K - S_{i,j}>e^{-rt}(K-S_{i,j}(pu-(1-p)d))$$
$$K - S_{i,j}>e^{-rt}(K-S_{i,j}(\frac{e^{rt}-d}{u-d}u+(1-\frac{e^{rt}-d}{u-d})d))$$
$$K - S_{i,j}>e^{-rt}(K-S_{i,j}(\frac{e^{rt}-d}{u-d}u-(\frac{u-d}{u-d}-\frac{e^{rt}-d}{u-d})d))$$
$$K - S_{i,j}>e^{-rt}(K-S_{i,j}(\frac{e^{rt}u-ud}{u-d}+\frac{ud-e^{rt}d}{u-d}))$$
$$K - S_{i,j}>e^{-rt}(K-S_{i,j}e^{rt})$$
$$1 >e^{-rt}$$

This means that the inequality can hold for American puts.

### 4. European and American Call option calculations

In cases where full tree induction is unnecessary, the option price is computed as the discounted expected payoff (shown here for calls):


$$f = e^{-rT} \sum_{j>\alpha}^{n} \binom{n}{j} p^{n-j} q^{j} \cdot (S_0u^{n-j}d^{j}-K)$$
where $\alpha$ is the previously calculated k-index.

Direct computation of $\binom{n}{j}$ can cause integer overflow, so logarithms are used:


$$\ln(\binom{n}{j}) = \ln(n!) - \ln(j!(n-j)!)= \ln(n!) - \ln(j!) - \ln((n-j)!)$$

As that we take the logarithm of the important parts of expression:
$$p^*=\ln(\sum_{j>\alpha}^{n} \binom{n}{j} p^{n-j} q^{j})= \sum_{j>\alpha}^{n}\ln(\binom{n}{j}) + (n-j)\ln(p) + j\ln(q)$$

We use the log-gamma function (`math.lgamma`) to compute factorial logarithms efficiently. Since
$\Gamma(x) = (x-1)!$, we use $n+1$ as input.

Because payoffs are precomputed, the final price is:
$$f=e^{-rT}\sum_{j>\alpha}^n\exp(p^*)*(\text{Payoff}_j)$$

This reduces the time complexity from $O(N^2)$ to $O(N)$.

## Disclaimer

Large Language Models (LLMs) were used during the development of this project for the following purposes:
- Correcting grammatical and stylistic errors in the documentation.
- Assisting in the creation of benchmark algorithms used for performance comparisons.

All core implementation logic, modeling decisions, and final code integration were reviewed and validated by the author.
