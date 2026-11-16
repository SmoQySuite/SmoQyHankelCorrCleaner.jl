# Algorithm: Hankel projections with Dykstra's method

This page describes the mathematics behind the denoising algorithm. Every step in
Dykstra's iteration must be an *exact* projection in the Frobenius norm. This page explains
how each step meets that requirement, including the constraints that act only on
anti-diagonal averages.

## 1. Setting

The input is imaginary-time correlation data $G_l = G(\tau_l)$ for $l = 1, \dots, L$, where
$L = 2N - 1$. The data is embedded in an $N \times N$ Hankel matrix

$$
H_{ij} = G_{i+j-1}, \qquad 1 \le i, j \le N .
$$

The $l$-th anti-diagonal is the set of entries with $i + j - 1 = l$. It contains

$$
n_l = \min(l,\; L - l + 1)
$$

entries, and its average is

$$
a_l(H) = \frac{1}{n_l} \sum_{i + j - 1 = l} H_{ij} .
$$

## 2. Dykstra's algorithm

Let $C_1, \dots, C_K$ be closed convex sets of $N \times N$ matrices, and let $P_k$ be the
Frobenius-nearest-point projection onto $C_k$:

$$
P_k(H) = \operatorname*{arg\,min}_{X \in C_k} \lVert X - H \rVert_F^2 .
$$

Starting from $H_0$, Dykstra's algorithm applies the projections cyclically. It keeps one
correction matrix $Q_k$ per set:

$$
\begin{aligned}
X &\leftarrow P_k(H + Q_k), \\
Q_k &\leftarrow H + Q_k - X, \\
H &\leftarrow X .
\end{aligned}
$$

This converges to the nearest point of $C_1 \cap \dots \cap C_K$ to $H_0$:

$$
H^\star = \operatorname*{arg\,min}_{X \in \cap_k C_k} \lVert X - H_0 \rVert_F^2 .
$$

The guarantee requires every $P_k$ to be an exact projection **in the same inner product**,
here the Frobenius inner product. Suppose one step instead minimizes some other quadratic
distance. The iteration may still converge, but possibly to a point of the intersection
that is not the nearest one, or it may stall.

## 3. Constraint sets and their projections

### 3.1 Positive semidefinite (full and shifted)

The physical correlation function gives a PSD Hankel matrix, so we impose

$$
C_{\text{psd}} = \{ X : X = X^{\mathsf T},\; X \succeq 0 \}.
$$

The projection symmetrizes, $S = (H + H^{\mathsf T})/2$, then eigendecomposes
$S = U \operatorname{diag}(\epsilon) U^{\mathsf T}$ and clips the eigenvalues:

$$
P_{\text{psd}}(H) = U \operatorname{diag}\big(\max(\epsilon, 0)\big)\, U^{\mathsf T}.
$$

The *shifted* Hankel matrix $H'_{ij} = G_{i+j}$ must also be PSD. It is the submatrix
`H[1:end-1, 2:end]`. The constraint involves only those entries, and the Frobenius norm
separates over disjoint sets of entries. So the exact projection applies $P_{\text{psd}}$ to
that block and leaves the first column and last row unchanged.

### 3.2 Hankel structure (with optional fixed corners)

$$
C_{\text{hankel}} = \{ X : X_{ij} \text{ depends only on } i + j \},
$$

optionally with $X_{11} = e_1$ and $X_{NN} = e_L$. Within each anti-diagonal, the constant
closest to the entries $h_k$ is their mean, $\operatorname*{arg\,min}_c \sum_k (c - h_k)^2 = a_l(H)$.
So the projection replaces each anti-diagonal by its average and then sets the two corners.
The corner anti-diagonals each have a single entry, so fixing them is still exact.

### 3.3 Constraints on anti-diagonal averages (general recipe)

Several constraints restrict only the vector of averages $a(X) = (a_1(X), \dots, a_L(X))$:

$$
C = \{ X : a(X) \in K \}
$$

for some closed convex $K \subseteq \mathbb{R}^L$. These sets do **not** require $X$ to be
Hankel. To project onto such a set, decompose each anti-diagonal of $X$ and $H$ into its mean
plus zero-sum deviations:

$$
x_k = a_l(X) + \delta_k, \qquad h_k = a_l(H) + \eta_k, \qquad \sum_k \delta_k = \sum_k \eta_k = 0 .
$$

A constant vector is orthogonal to every zero-sum vector. Pythagoras therefore gives, for each
anti-diagonal,

$$
\sum_{k=1}^{n_l} (x_k - h_k)^2 = n_l \big(a_l(X) - a_l(H)\big)^2 + \sum_{k=1}^{n_l} (\delta_k - \eta_k)^2 .
$$

The constraint involves only the means, so the deviations are free. The best choice is
$\delta = \eta$, which keeps the within-anti-diagonal structure of $H$ unchanged. The
projection then reduces to a weighted problem in $L$ variables:

$$
f^\star = \operatorname*{arg\,min}_{f \in K} \; \sum_{l=1}^{L} n_l \,\big(f_l - a_l(H)\big)^2 ,
$$

followed by a uniform shift of each anti-diagonal:

$$
P_C(H)_{ij} = H_{ij} + \big(f^\star_l - a_l(H)\big), \qquad l = i + j - 1 .
$$

The weights $n_l$ appear because moving the mean of anti-diagonal $l$ by $\Delta$ shifts
$n_l$ entries by $\Delta$ each. That costs $n_l \Delta^2$ in squared Frobenius norm, not
$\Delta^2$. Equivalently, $\sqrt{n_l}\, a_l$ are orthonormal coordinates for the Hankel
subspace. Dropping the weights would give an exact projection in a different metric, one
that over-weights the short anti-diagonals near the corners. That projection is incompatible
with the other, Frobenius-exact steps.

### 3.4 Positive curvature (with optional fixed endpoints)

The averages must form a convex sequence:

$$
K_{\text{curv}} = \{ f \in \mathbb{R}^L : f_{l+1} - 2 f_l + f_{l-1} \ge 0,\; l = 2, \dots, L-1 \},
$$

optionally intersected with $\{f_1 = e_1,\; f_L = e_L\}$. Applying §3.3 gives a small convex
quadratic program:

$$
\begin{aligned}
\min_{f \in \mathbb{R}^L} \quad & \sum_{l=1}^{L} n_l \,\big(f_l - a_l(H)\big)^2 \\
\text{s.t.} \quad & f_{l+1} - 2 f_l + f_{l-1} \ge 0, \qquad l = 2, \dots, L-1, \\
& f_1 = e_1, \quad f_L = e_L \quad \text{(if endpoints are fixed)},
\end{aligned}
$$

followed by the uniform per-anti-diagonal shift. The weights do not affect the endpoint
constraints, because $n_1 = n_L = 1$.

### 3.5 Symmetry about $\tau = \beta/2$ (optional)

$$
K_{\text{sym}} = \{ f : f_l = f_{L+1-l} \} .
$$

Since $n_l = n_{L+1-l}$, the weighted problem in §3.3 is solved by the plain pairwise average

$$
m_l = \tfrac{1}{2}\big(a_l(H) + a_{L+1-l}(H)\big),
$$

and the projection shifts anti-diagonal $l$ uniformly by $m_l - a_l(H)$.

### 3.6 Merging Hankel and average constraints (optional optimization)

The set $C_{\text{hankel}} \cap \{a(X) \in K\}$ also has a closed-form projection. For a Hankel
$X$ with values $f_l$, the cost is

$$
\sum_{l} \sum_{k=1}^{n_l} (f_l - h_k)^2 = \sum_l n_l \big(f_l - a_l(H)\big)^2 + \text{const},
$$

so the projection solves the same weighted problem for $f^\star$ and sets $X = \operatorname{hankel}(f^\star)$.
Merging these steps reduces the number of Dykstra sets, which typically speeds up convergence.

## 4. Summary of the algorithm

1. Build $H_0 = \operatorname{hankel}(G^{\text{noisy}})$.
2. Run Dykstra's algorithm, starting from $H_0$, over these sets, each with its exact
   Frobenius projection:

   | Set | Exact projection |
   |:----|:-----------------|
   | Full PSD | Symmetrize, eigendecompose, clip negative eigenvalues to zero |
   | Shifted PSD (`H[1:end-1, 2:end]`) | Same, applied to that block only |
   | Hankel (+ fixed corners) | Replace each anti-diagonal by its average; set $H_{11} = e_1$, $H_{NN} = e_L$ |
   | Convex averages (+ endpoints), optional | Weighted QP of §3.4, then uniform shift per anti-diagonal |
   | Symmetric averages, optional | Shift anti-diagonal $l$ by $m_l - a_l(H)$ |

3. Stop when the change between iterations is below `tol`, or after `maxiter` iterations.
4. Return $G^{\text{clean}}_l = a_l(H_{\text{final}})$.

## 5. Implementation notes

The functions in `src/projections.jl` correspond to the sets above:

- `project_psd!`: §3.1. The shifted variant is a view `X[1:end-1, 2:end]`.
- `project_hankel!`: §3.2.
- `project_positive_curvature!`: §3.4. It solves the weighted QP with JuMP and HiGHS, then
  applies the uniform shift. The weights $n_l$ in the objective are what make the step an
  exact Frobenius projection. The QP iteration count is capped, because HiGHS can cycle
  on degenerate inputs.
- `project_symmetrized!`: §3.5. It shifts whole anti-diagonals, so the step is exact even
  when $H$ is not Hankel.
