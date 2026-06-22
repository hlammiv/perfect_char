# Literature Review: Perfect Actions on Discrete Subgroups of SU(N) via Multiple Characters

## 1. Executive summary

The continuum limit of lattice gauge theory is reached by sending the bare coupling to its critical value while removing lattice-spacing (cutoff) artifacts. Two complementary programs control these artifacts: **Symanzik improvement** [1,2], which cancels powers of the lattice spacing $a$ order by order with tuned higher-dimension operators, and the **perfect / fixed-point (FP) action** program [3,4,5], which targets the renormalization-group (RG) renormalized trajectory where cutoff effects vanish to all orders. For asymptotically free theories the FP action is computable nonperturbatively as a *classical* saddle-point problem [3], and for SU(3) it has been constructed and shown to be classically (and 1-loop quantum) perfect [5].

A separate, originally pragmatic, idea replaces the continuous gauge group by a **finite subgroup** so that each link is a small integer label rather than a floating-point matrix [9,10,11]. This always fails at weak coupling: because no compact non-Abelian Lie group admits an asymptotically dense chain of discrete subgroups, a fixed finite subgroup with the Wilson action undergoes a first-order **freezing (bulk) transition** before the scaling window [9,10,11,12,18]. For the largest crystal-like SU(3) subgroup, $S(1080)=\Sigma(360\times3)$ (the Valentiner group $3.A_6$), the plain Wilson action freezes at $\beta_f = 3.935(5)$ [22], far below the scaling region $\beta\gtrsim6$.

The modern revival of finite groups is driven by **quantum and tensor-network simulation**, where a finite per-link Hilbert space is mandatory [22,23,30]. The established remedy for freezing is to **modify the single-plaquette action with extra group characters**, each with its own coupling: $S_{\rm plaq}=\sum_r\beta_r\,\mathrm{Re}\,\chi_r(U_p)$ [13,21,24,25,26]. A negative coupling on $\mathrm{Re}\,\mathrm{Tr}\,U_p^2$ already pushes $S(1080)$ to $a\approx0.08$ fm and reproduces the SU(3) glueball spectrum and deconfinement temperature to percent level [22,23].

**Where this project sits.** This project pursues a *perfect (or strongly improved) discrete action*: it treats the character-space couplings $\{\beta_c\}$ on a single fixed subgroup (e.g. $S(1080)$ or $S(108)$) as a multi-parameter manifold and tunes many of them simultaneously to approximate continuum SU(3). The gap it targets is that prior work has (i) added essentially *one* extra character (adjoint, or $\mathrm{Tr}\,U_p^2$) [13,24], or (ii) fixed several couplings *perturbatively* by character-expansion / group-space decimation matching [25,26], or (iii) used subduction as a *diagnostic* of which terms to add [21]. What is not yet done systematically is a full **variational/perfect-action optimization over many character couplings on one fixed subgroup, tuned against multiple SU(3) observables at once** (Casimir scaling across irreps, glueballs, $T_c$). That is the natural "perfect-character action" framing for this project, building directly on the perfect-action lineage [3,4,5,6,7,8] and the digitization line [21,22,23,24,25,26].

## 2. Perfect and fixed-point actions

### 2.1 Wilson's RG and the perfect action

The conceptual foundation is Wilson's block-spin renormalization group [1 — *to be re-numbered; see References*]. A block-spin transformation maps a fine lattice to a coarser one while preserving long-distance physics; the bare action flows in an infinite-dimensional space of couplings. A critical system stays critical under blocking, so its action flows along the **critical surface** toward an RG **fixed point (FP)**; the unstable manifold leaving the FP is the **renormalized trajectory (RT)**. Actions sitting exactly on the RT are *perfect*: they reproduce continuum physics with **no cutoff artifacts at all**, even at correlation lengths of $O(a)$. Wilson and Kogut already noted that perfect actions exist *off* the critical surface (along the massive RT), not only at criticality [W1]. Because the RT is generically nonlocal, the practical question — open for roughly two decades — was whether a *local, parametrizable* approximation could be simulated.

### 2.2 The Hasenfratz–Niedermayer construction

The decisive advance [H1] exploits the fact that for an **asymptotically free** theory the critical surface lies at bare coupling $\beta=\infty$ ($g=0$). In that limit the blocking recursion — a functional integral over fine-lattice fluctuations — collapses to a **saddle-point / minimization equation**, a purely *classical* field-theory problem:
$$ S^{\rm FP}(V)=\min_{U}\big[\,S^{\rm FP}(U)+T(U,V)\,\big], $$
where $T$ is the blocking kernel, $V$ the coarse and $U$ the fine field. Because this is classical minimization, $S^{\rm FP}$ can be computed **nonperturbatively without Monte Carlo**, by analytic work plus deterministic minimization on small lattices. The resulting FP action is **classically perfect**: rotational symmetry is restored, the free spectrum is exact, classical solutions are artifact-free — including **scale-invariant lattice instantons**, which the Wilson action lacks. A free parameter in the blocking kernel can be tuned so the interaction is **short-ranged (local)**, with couplings decaying exponentially. In the 2-d O(3) $\sigma$-model, cutoff effects were nearly eliminated down to correlation lengths of only $\sim3$ lattice spacings [H1].

### 2.3 Locality, parametrization, and numerical tests

A perfect action is in principle an infinite tower of couplings; practical success hinges on **truncation** to a few short-range terms, optimized via the blocking-kernel smoothing parameter. For the O(3) model, numerical work verified excellent scaling of the mass gap and clean topological behavior with the parametrized FP action [DP]; a closely related study of FP instantons and topological charge is [BI]. (Note: the often-cited "Blatter et al." reference to *Nucl. Phys. B 530, 185* conflates two distinct papers — see References and the Unverified-leads remarks.)

### 2.4 Fixed-point actions for SU(N) gauge theory

The program was carried to non-Abelian gauge theory: the **classically perfect FP action for SU(3) pure gauge theory** was constructed in [H2], with scale-invariant instantons, an exact free spectrum, and no topological artifacts. It is argued to be **1-loop quantum perfect** — no $g^2 a^n$ cutoff effects — with the static $q\bar q$ potential free of inverse-power cutoff effects at leading order. A few-parameter parametrization and nonperturbative scaling tests (string tension from torelon masses at fixed physical volume) showed the approximate FP action scales within statistical errors where the Wilson action shows $\sim$10% violations [H3].

### 2.5 Blocking from the continuum; quantum perfection and Ginsparg–Wilson

A complementary, fully solvable route is **blocking directly from the continuum** [BW]. For free (Gaussian) fields the RT is computed analytically: averaging continuum fields over hypercubes yields the **perfect free fermion and free gluon actions**. Tuning the RG-transformation parameter, the perfect free quark action reduces to nearest-neighbor Wilson fermions (1-d configs) and the free gluon FP action reduces to the plaquette action (2-d configs); classically perfect composite fields, currents and vertices were also constructed [BW].

Classically perfect actions remove tree-level artifacts; a fully **quantum perfect** action sits exactly on the RT. The most influential offshoot came when Hasenfratz recognized that **FP fermion actions automatically satisfy the Ginsparg–Wilson relation**, realizing exact lattice chiral symmetry, the correct anomaly, and a lattice index theorem with no fine-tuning [H4] — reviving the GW relation and connecting perfect actions to overlap/domain-wall fermions. A comprehensive review of the program is [H5].

## 3. Symanzik improvement and improved gauge actions

### 3.1 The Symanzik program

Symanzik's two-part 1983 work [S1] organizes lattice cutoff effects via an **effective continuum action**, $S_{\rm eff}=S_{\rm cont}+a\,S_1+a^2 S_2+\dots$, where each $S_k$ is an integral of local operators of dimension $4+k$ respecting the lattice symmetries. Adding a complete basis of higher-dimension irrelevant operators with tuned coefficients cancels the $S_k$ order by order. The program is systematic but open-ended — each order needs more, higher-dimension operators.

### 3.2 On-shell improvement: Lüscher–Weisz

Implementing improvement in gauge theory is subtle because of operator proliferation and the role of the equations of motion. Lüscher and Weisz [LW] introduced **on-shell improvement**: improve only gauge-invariant spectral quantities (stable-particle masses, the static potential, scattering amplitudes), so operators differing only by the equations of motion may be dropped. The leading $O(a^2)$ improvement of SU(N) needs, besides the $1\times1$ plaquette, three six-link loops — the planar $2\times1$ **rectangle**, the bent **chair**, and the 3-d **parallelogram** — with tree-level coefficients $c_0=5/3$, $c_1=-1/12$ (rectangle), $c_2=c_3=0$; the one-loop quantum coefficients were computed in weak coupling [LW]. The widely used tree-level **Lüscher–Weisz** action keeps plaquette + rectangle.

### 3.3 Tadpole (Lepage–Mackenzie) improvement

A largely orthogonal idea [LM] traces poor convergence of *bare* lattice perturbation theory to large unphysical **tadpole** contributions and to using $g_0^2$ as the expansion parameter. The remedy is (i) **mean-field/tadpole improvement** of links by dividing each by the mean link $u_0=\langle\frac1N\mathrm{Re}\,\mathrm{Tr}\,U_p\rangle^{1/4}$ (or the Landau-gauge link), and (ii) re-expanding in a **renormalized coupling** defined from a physical quantity (e.g. $\alpha_V$ from the static potential). Tadpole improvement is cheap and routinely combined with Symanzik improvement. A complementary analytic route via the Cayley parametrization is [PE].

### 3.4 RG-improved actions: Iwasaki and DBW2

A different philosophy builds improved actions from the RG rather than from a fixed order in $a$. Both standard RG actions use plaquette + $2\times1$ rectangle (normalized $c_0+8c_1=1$):

- **Iwasaki** [IW]: from a block-spin RG analysis aiming to sit near the RT after blocking; $c_1\approx-0.331$. It improves rotational invariance of the static potential and suppresses Wilson-fermion artifacts.
- **DBW2** ("Doubly Blocked Wilson in 2-coupling space") [TARO]: from the nonperturbative RG flow of SU(3) in the (plaquette, rectangle) plane (Swendsen blocking + Schwinger–Dyson); the flow collapses onto a narrow attractive stream, and DBW2 takes a strongly negative $c_1\approx-1.4069$. It strongly suppresses chiral-symmetry-breaking artifacts at the cost of larger fluctuations.

These RG actions are heuristic approximations to a genuine perfect action on the RT.

### 3.5 The bridge: perfect actions are all-orders on-shell Symanzik improved

Because an FP action removes cutoff effects on the equations of motion to all orders, it **coincides with an on-shell tree-level Symanzik-improved action to all orders in $a^2$** [PA; also 5,BW]. A perfect action is thus the all-orders endpoint of the (on-shell, classical) Symanzik program, reached nonperturbatively rather than order by order. In practice FP/perfect actions must be **truncated** to short-range parametrizations ("quasi-perfect" actions) [QP], reintroducing residual artifacts — placing truncated-perfect and finite-order Symanzik-improved actions in the same practical family. **For this project, this is the key conceptual hinge:** a multi-character single-plaquette action and a finite-order improved action are two parametrizations of the same approximately-perfect target.

## 4. Discrete subgroups in lattice Monte Carlo and the freezing problem

### 4.1 The idea and its fatal limitation

Replacing a continuous gauge group by a finite subgroup makes the Wilson action a lookup table on a small set of group labels — appealing for cheap Monte Carlo, and now essential for quantum/tensor-network registers. The fatal limitation, recognized immediately, is the **freezing (bulk) transition**: for any compact Lie group other than U(1) there is **no asymptotically dense chain of discrete subgroups**, so a fixed finite subgroup reproduces the continuum only down to a finite minimum lattice spacing / up to a finite $\beta$ [9,18].

### 4.2 Early SU(2) and SU(3) studies

- **Rebbi** [RE] performed Monte Carlo on 4-d theories with finite non-Abelian groups, finding a generic **two-phase structure with a first-order transition**; for SU(2) subgroups the transition moves to weaker coupling (larger $\beta$) as the group order grows.
- **Petcher–Weingarten** [PW] gave the definitive SU(2) study, simulating the SU(2) images of the tetrahedral, octahedral and icosahedral symmetry groups — the **binary tetrahedral (24)**, **binary octahedral (48)** and **binary icosahedral (120)** groups. Each shows a first-order freezing transition; the **binary icosahedral group (order 120) tracks SU(2) furthest into weak coupling**, making it the best fixed SU(2) subgroup. They gave an analytic model reproducing the critical couplings.
- For SU(3), **Bhanot–Rebbi** [BR] studied five crystal-like subgroups $S(60),S(108),S(216),S(648),S(1080)$ with the Wilson action; **all** freeze (first-order order–disorder) before the scaling window, so none can reach the continuum. **Grosse–Kühnelt** [GK] independently confirmed bulk transitions for non-Abelian SU(3) subgroups.

### 4.3 Physics of freezing and the critical coupling

The mechanism is geometric: a Metropolis proposal is drawn from a finite set, so there is a lower bound on the achievable $\Delta S$. As $\beta\to\infty$ the Boltzmann factor suppresses any update above that gap, acceptance collapses, and links freeze onto the action-minimizing element. The closed form for the SU(2) subgroup critical coupling is [HA, Eq. (30)]:
$$ \beta_c(N)\approx\frac{\ln(1+\sqrt2)}{1-\cos(2\pi/N)}, $$
with $N$ the cyclic order of the nearest neighbour of the identity. This gives $\beta_c\approx1.8$ ($N=6$), $3.0$ ($N=8$), $4.6$ ($N=10$), $6.6$ ($N=12$), consistent with the binary T/O/I ordering ($T<O<I$, order-120 latest) [HA,JA].

### 4.4 Strategies to mitigate freezing

1. **Modified / extended (Symanzik-like) actions.** Add extra Wilson-loop or higher-representation terms sharing the continuum limit but lowering the action gap between competing discrete configurations. For $S(1080)$, Alexandru–Bedaque–Brett–Lamm used $S=-\sum_p(\frac{\beta_0}{3}\mathrm{Re}\,\mathrm{Tr}\,U_p+\beta_1\mathrm{Re}\,\mathrm{Tr}\,U_p^2)$ with **negative $\beta_1$**; along $\beta_1=-0.1267\,\beta_0+0.253$ they reached $a\approx0.08$ fm and reproduced the SU(3) glueball spectrum and deconfinement temperature to percent level [23].
2. **Mixed fundamental–adjoint actions.** The classic route to control the bulk transition (Section 7).
3. **Group-space decimation (effective actions).** Integrate out the off-subgroup part of SU(3) to derive an effective action on the 1080-element subgroup that mimics SU(3) deeper into weak coupling [FL] — revived for quantum computing as group-space decimation [25].
4. **Abandon group structure: asymptotically dense subsets.** Since no dense subgroup exists, Hartung et al. discretize SU(2)$\cong S^3$ with non-group point sets (linear/octahedral, "volleyball"/hypercube, and **Fibonacci-spiral** lattices, with Voronoi-volume weights for the non-Fibonacci sets); an (unweighted) generalized Fibonacci lattice distributes points most isotropically and gives the largest unfrozen $\beta$-range, allowing essentially **unfrozen Monte Carlo at arbitrary $\beta$** — at the cost of losing closure under multiplication (a problem for Hamiltonian/gate implementations) [HA,JA]. (Voronoi weighting is applied to the linear/volleyball partitionings, not the Fibonacci ones, where it actually lowers $\beta_c$ [JA].)

**Bottom line:** fixed finite subgroups always freeze at finite coupling; freezing is *pushed back* — not removed — by modified/extended, mixed, and decimation actions, and *removed* only by dropping the group structure.

## 5. Finite subgroups of SU(3)

### 5.1 Classification

The classification originates with Blichfeldt and is reviewed/corrected in the modern flavour-physics literature [GL,LU]. Finite subgroups of SU(3) fall into: (A) abelian (diagonal) groups; (B) groups embedding SU(2)/U(2); the infinite series (C) $\cong(\mathbb Z_m\times\mathbb Z_{m'})\rtimes\mathbb Z_3$ and (D) $\cong(\mathbb Z_n\times\mathbb Z_{n'})\rtimes S_3$ (the $\Delta(3n^2)$, $\Delta(6n^2)$ series are special cases); and a finite list of **exceptional** groups $\Sigma(60)\cong A_5$, $\Sigma(168)\cong\mathrm{PSL}(2,7)$, $\Sigma(36\times3)$, $\Sigma(72\times3)$, $\Sigma(216\times3)$, $\Sigma(360\times3)$ [GL,LU]. The notation $\Sigma(n\varphi)$ means order $n\cdot\varphi$ with $\varphi$ the order of the center; groups with a nontrivial $\mathbb Z_3$ center are the $\Sigma(N\times3)$ groups, of order $3N$.

### 5.2 The lattice S(N) dictionary

The lattice-gauge community labels each group by its order $N$. The dictionary is: $S(216)=\Sigma(72\times3)$ (order 216, verified directly in a 2025 quantum-gates paper [OP]), $S(648)=\Sigma(216\times3)$ (order 648), $S(1080)=\Sigma(360\times3)$ (order 1080). The order-108 lattice group is identified with $\Sigma(36\times3)$ (order 108) by some authors — see open questions; confirm against the project's `groups/S108` file.

### 5.3 The Valentiner group $S(1080)$

$\Sigma(360\times3)$ of order 1080 is the **Valentiner group**, the perfect (quasisimple) triple cover $3.A_6$ of the alternating group $A_6$, with center $\mathbb Z_3$ and inner automorphism group $A_6$ (SmallGroups id [1080,260]) [VA]. It has four faithful 3-dimensional irreps (two complex-conjugate pairs), plus faithful 6- and 15-dimensional complex-conjugate pairs [GL]. It is the **largest crystal-like discrete subgroup of SU(3)** [BR,9].

### 5.4 Character tables

The most directly useful modern reference for lattice work is Assi–Lamm [21], giving complete character tables for the four non-Abelian crystal-like $\mathbb Z_3$-center subgroups $\Sigma(108),\Sigma(216),\Sigma(648),\Sigma(1080)$ (Tables IV–VII). Salient structure:
- **$\Sigma(108)$:** 14 conjugacy classes/irreps; four singlets, several 3-dim irreps (and conjugates), two 4-dim irreps. Only the fundamental **3** subduces one-to-one to an SU(3) irrep.
- **$\Sigma(216)$:** singlets, a **2**, eight **3**'s, a complex pair **6/6***, and an **8** ($4\cdot1+2^2+8\cdot3^2+2\cdot6^2+8^2=216$).
- **$\Sigma(648)$:** 1's, 2's, many 3's and 6's, three **8**'s, and **9/9*** — relevant because the SU(3) adjoint is 8-dimensional.
- **$\Sigma(1080)$:** irreps of dimension 1, 3 (four faithful 3's), 5, 6, 8, 9, 10, 15, with characters involving $\mu_{1,2}=(1\mp\sqrt5)/2$ (golden-ratio entries reflecting the $A_6$/icosahedral content).

Character tables for exceptional groups also appear in the flavour literature [GL,LU] and machine-readable GAP/SmallGroups listings [LG].

### 5.5 As continuum approximants

With the plain Wilson action, every crystal-like SU(3) subgroup freezes before scaling [BR], located at $\beta_f=3.935(5)$ for $S(1080)$ [9]. A **modified/improved action is therefore mandatory** for continuum-relevant use [9,23]. Assi–Lamm systematize *why* via **subduction**: observables, being class functions, deviate most for SU(3) irreps that subduce to *direct sums* of subgroup irreps; Casimir scaling of static potentials in $\Sigma(360\times3)$ agrees at percent level for one-to-one-subducing irreps (e.g. the fundamental) but deviates significantly for others — a diagnostic of digitization quality that also guides which extra action terms to add [21].

## 6. Character expansion and the single-plaquette action manifold

### 6.1 General framework

A single-plaquette action depends on $U_p=U_1U_2U_3^{-1}U_4^{-1}$; gauge invariance forces the Boltzmann weight to be a **class function** of $U_p$. By the Peter–Weyl theorem, any square-integrable class function on a compact group expands in **irreducible characters** $\chi_r(U)=\mathrm{Tr}\,D_r(U)$:
$$ e^{-S(U)}=\sum_r d_r c_r\,\frac{\chi_r(U)}{\chi_r(1)}\quad\Longleftrightarrow\quad S(U)=\sum_r\beta_r\,\chi_r(U), $$
with $d_r=\chi_r(1)$ the dimension. Character **orthogonality**, $\int_G\chi_r(U)\chi_s(U)^*\,dU=\delta_{rs}$, projects out the coefficients, $\beta_r=\int_G e^{-S(U)}\chi_r(U)^*\,dU$. This is the cornerstone of the strong-coupling expansion, reviewed by Drouffe–Zuber [DZ] and rooted in the lattice gauge series of Balian–Drouffe–Itzykson [BDI1,BDI3]. For a finite group, the number of independent character couplings equals the number of conjugacy classes.

### 6.2 Wilson action and Bessel coefficients

For SU(2) Wilson, $\exp[(\beta/2)\chi_{1/2}(U)]=\sum_j(2/\beta)(2j+1)I_{2j+1}(\beta)\,\chi_j(U)$, with $I_n$ modified Bessel functions [TO,DZ]. The normalized coefficients $D_j=I_{2j+1}/I_1$ decay rapidly at small $\beta$, suppressing higher representations — the expansion parameter is effectively $D_{1/2}\sim\beta/4$. For SU(N) the coefficients are ratios of Bessel-function determinants.

### 6.3 Heat-kernel (Villain) action

The **heat-kernel action** is the diffusion kernel on the group manifold, with exactly factorized coefficients $c_r=\exp[-\tfrac{t}{2}C_2(r)]$ ($C_2$ the quadratic Casimir) [MO,DZ]. It is the **fixed point of the Migdal recursion** — integrating out the shared link of two adjacent plaquettes reproduces the same functional form because $e^{-tC_2}$ factors are multiplicative under convolution [MI]. This makes 2-d Yang–Mills exactly solvable: $Z=\sum_r(d_r)^\chi e^{-(g^2A/2)C_2(r)}$ ($\chi$ the Euler character). The Villain action is the U(1) special case. The heat-kernel action coincides with Wilson's only as $\beta\to\infty$.

### 6.4 The manifold of single-plaquette actions

$S(U)=\sum_r\beta_r\chi_r(U)$ shows the space of gauge-invariant single-plaquette actions is an (infinite-dimensional) **manifold parametrized by the representation couplings $\{\beta_r\}$**. Wilson keeps only the fundamental; the heat-kernel action fixes all $\beta_r$ by a single coupling via the Casimir; **mixed actions** retain a finite set. The canonical example is the Bhanot–Creutz mixed fundamental–adjoint action (Section 7). Moving on this manifold is exactly how improved and perfect actions add couplings to higher representations to remove cutoff artifacts while staying on or near the RT [H1,S1]. **This is precisely the manifold this project explores for a fixed discrete group.**

### 6.5 Convergence and higher representations

The character/strong-coupling series **converges with finite radius**, bounded by the nearest singularity — typically the bulk/deconfinement coupling $\beta_c$ [DZ]. At leading order the Wilson loop obeys an area law $W\sim e^{-\sigma A}$ with $\sigma=-\ln D_f(\beta)$ set by the fundamental coefficient [BDI3,DZ]; higher representations decorate the minimal tiling at higher order. The rapid decay of $\beta_r$ in the rep label controls the expansion; its breakdown near $\beta_c$ signals the scaling/continuum regime. For finite groups there is no continuous Laplacian, so the analogue of a "heat-kernel" coefficient (a Casimir-weighted form) is not immediate — an open question for defining a "perfect character" target on a discrete group.

## 7. Mixed and multi-representation actions

### 7.1 The fundamental–adjoint plane and the Bhanot–Creutz diagram

The canonical mixed action adds the adjoint character to the fundamental. The defining study is **Bhanot–Creutz** [BC], which (for SU(2)) follows the first-order transitions of the pure SO(3) and pure $\mathbb Z(2)$ theories inward in the $(\beta_f,\beta_a)$ plane: the two lines **merge at a triple point and continue as a single first-order line that terminates at a critical endpoint** near $(\beta_f,\beta_a)\approx(1.5,0.9)$. The "specific-heat peak" / crossover on the pure-Wilson axis is the **shadow of this nearby bulk singularity** — a lattice artifact, not a continuum transition; one can deform the action *around* the endpoint and reach the continuum without crossing any singularity. (The paper is specific to SU(2) with SU(2)/SO(3)/$\mathbb Z(2)$ limits; the generic SU(N)/$\mathbb Z(N)$ form is a standard later generalization.)

### 7.2 SU(3) in the fundamental–adjoint plane

The same structure recurs for SU(3): a first-order bulk line runs from the pure-adjoint axis toward larger $\beta_f$ and ends at a critical endpoint at $(\beta_f,\beta_a)\approx(4.00(7),2.06(8))$ [HE]. Measuring the string tension and glueball masses along the $N_t=4$ thermal line, Heller found $\sqrt\sigma/T_c$ roughly constant while only the $0^{++}$ glueball mass softens near the endpoint — consistent with an ordinary critical point whose continuum limit is a free $\phi^4$ scalar with a diverging correlation length *only* in the $0^{++}$ channel. A key diagnostic: for thin lattices ($N_t=2$) the thermal deconfinement line joins the bulk line, but for $N_t\geq4$ the deconfinement line **splits cleanly away** and moves to weaker coupling as $N_t$ grows — distinguishing a genuine finite-$T$ transition from the bulk artifact [HE].

### 7.3 Order of the endpoint, and center defects

Whether the bulk endpoint is *ordinary* critical or *tricritical* has been debated. Grady [GR] argues from scaling along the first-order line that the SU(2) endpoint is **tricritical**, implying a second-order bulk line continues to the Wilson axis (which would make the conventional finite-$T$ SU(2) transition partly bulk in origin). The competing picture attributes the non-universal behaviour to lattice $\mathbb Z(2)$ **monopoles and vortices**: Gavai–Mathur showed that **suppressing $\mathbb Z(2)$ monopoles restores universality** [GM]. In the Villain form with monopoles suppressed, the SU(2) mixed action is **self-dual under $\beta_f\leftrightarrow\beta_a$**, and the deconfinement transition stays in the 3-d Ising class throughout the plane in the continuum; the argument generalizes to SU(N) with $\mathbb Z(N)$ defects [GM]. Improved mean-field and strong-coupling analyses confirm the *qualitative topology* of the bulk diagram is robust, though the endpoint location is sensitive to higher-order terms [DG].

### 7.4 Earliest SU(3) multi-character action

The earliest concrete SU(3) multi-character action is **Bhanot 1982** [BH], adding the 8-dimensional adjoint character to the fundamental Wilson term and mapping the resulting phase structure/artifacts.

### 7.5 Lesson for this project

The extended coupling space is both a **diagnostic** (the bulk first-order surface and its endpoint bound where continuum extrapolations are valid — improved actions must stay on the weak-coupling side and route around the endpoint) and a **constructive tool** (tuning adjoint/higher-rep couplings suppresses center-defect artifacts and smooths the approach to asymptotic scaling). The SU(3) endpoint coordinates and the $N_t$-dependent splitting of thermal vs. bulk lines are practical diagnostics — directly relevant to distinguishing physical signals from the noted noisy $S(1080)$-improved phase diagram at large $|\beta_2|$.

## 8. Synthesis: multi-character actions as a route to a perfect discrete action

The threads above converge on a single picture relevant to this project.

1. **The target.** The renormalized trajectory / perfect action is the cutoff-free ideal [W1,H1]; for SU(3) the FP action is explicitly constructed and is on-shell tree-level Symanzik improved to all orders [H2,PA]. In practice it is *truncated* to a short-range parametrization [QP].
2. **The basis.** Any gauge-invariant single-plaquette action is $\sum_r\beta_r\chi_r(U_p)$ [DZ,Peter–Weyl]. The heat-kernel action — the Migdal fixed point with Casimir-weighted coefficients — is the canonical "perfect" single-plaquette action and a natural benchmark [MO,MI].
3. **The constraint.** On a *finite* subgroup, the Wilson action freezes before scaling [BR,9,18]; the number of available characters is finite (= number of conjugacy classes). So the discrete single-plaquette manifold is *finite-dimensional*, and one must spend those couplings to (a) delay freezing and (b) approach the perfect/continuum action.
4. **The hazard.** Turning on extra-representation couplings generically opens first-order bulk lines (Bhanot–Creutz/Heller) [BC,HE]; continuum extrapolation must avoid them, and center-defect ($\mathbb Z(N)$) physics governs the endpoint structure [GM]. This plausibly underlies the noted phase-diagram noise at large $|\beta_2|$ for $S(1080)$-improved.

**Prior art that already does (or nearly does) this.** The ingredients are anticipated:
- **One extra character.** Bhanot 1982 (adjoint) [BH]; the $S(1080)$ studies add $\mathrm{Re}\,\mathrm{Tr}\,U_p^2$ with negative $\beta_1$ [9,23].
- **Principled multi-coupling recipes (perturbative).** Group-space decimation derives the discrete effective action by integrating out transverse fluctuations [FL,25]; character-expansion matching fixes couplings by matching continuum to discrete order-by-order, validated to fifth order in strong-coupling PT [26].
- **Diagnostics for which terms to add.** Subduction / cumulant methods identify the digitization-error-dominating irreps and measure quality via Casimir scaling across all $\Sigma(360\times3)$ irreps [21].
- **Modern perfect-action practice on continuous SU(3).** Machine-learned RG-improved gauge actions and classically perfect gradient flows show the perfect-action philosophy is current, but applied to continuous SU(3), not finite groups [HO].

**What remains novel.** No prior work appears to perform a full **variational / perfect-action optimization over many character couplings on one fixed subgroup** (e.g. $S(1080)$ or $\Sigma(360\times3)$), tuned *simultaneously against multiple SU(3) observables* — Casimir scaling across irreps, the glueball spectrum, and $T_c$ — i.e. treating the discrete character-space couplings as a multi-parameter "perfect discrete action," beyond a single adjoint/$\mathrm{Tr}\,U_p^2$ term and beyond perturbative matching. That is the gap this project can claim, building directly on the digitization line [9,23,25,26,21] and the FP/Symanzik framework [H1,H2,PA].

## 9. Modern motivation: digitization for quantum simulation

Quantum (and tensor-network) simulation requires a **finite Hilbert space per link**. Replacing a continuous group $G$ by a finite subgroup $H\subset G$ gives a $|H|$-dimensional link space with *exact* group multiplication, inversion, and gauge invariance, mapping naturally onto qubits/qudits — the original NuQS program using $S(1080)$ (11 qubits/link) [9] and the SU(2) analog using binary octahedral/icosahedral groups [HA]. The freezing obstruction (Sections 4–5) then *forces* improvement: $S(1080)$ with the Wilson action freezes at $\beta_f=3.935(5)$ before scaling [9], so the **modified multi-character action is not optional but mandatory** [9,23,25,26,21].

Alternative digitizations share the finite-link-space goal but differ in route: **non-subgroup asymptotically dense partitionings** (Fibonacci-spiral point sets) evade freezing at the cost of losing group closure (problematic for Hamiltonian/gate use) [HA]; **magnetic-basis / group-manifold partitionings** preserve the electric spectrum and approximate the fundamental commutation relations for tensor-network/Hamiltonian work [GB]; and truncated-irrep, $q$-deformed, and monomer–dimer schemes round out the landscape. Finite-subgroup digitization *plus multi-character improvement* is distinguished by exact gauge invariance and a clean classical-Monte-Carlo proof of principle: the modified $S(1080)$ action reproduces continuum SU(3) glueballs and $T_c$ at percent level [23].

## 10. Open questions and how this project can contribute

1. **Multi-coupling optimization.** Has anyone optimized *more than one* extra character coupling (beyond fundamental+adjoint) on a single fixed subgroup, against *multiple* observables simultaneously? This project can deliver exactly that — a variational perfect-character action — and should verify the gap against the full reference lists of [26,21].
2. **A discrete "heat-kernel"/perfect target.** For finite groups there is no continuous Laplacian; what plays the role of the quadratic Casimir in defining a "perfect" single-plaquette coefficient $\beta_r$? Can an RG-blocking transformation be defined *directly in character space* on a finite subgroup, or only matched perturbatively (decimation [25], character expansion [26])?
3. **How far can freezing be pushed?** Quantify how far multi-character improvement defers the $S(1080)$ freezing transition beyond $a\approx0.08$ fm [23], and whether the large-$|\beta_2|$ phase-diagram noise reflects proximity to a bulk transition (cf. Bhanot–Creutz/Heller endpoint structure [BC,HE]).
4. **Subduction-targeted couplings.** Use the subduced direct sums [21] to choose couplings that independently correct specific SU(3) irreps' Casimir scaling; test whether independent control of several irreps is achievable for $S(1080)$ vs $\Sigma(360\times3)$.
5. **Convergence and bulk lines.** Map how the convergence radius and induced first-order bulk lines move when several higher-representation couplings are turned on together, relative to the scaling region.
6. **Removal vs. deferral; reflection positivity.** Does multi-coupling improvement remove the bulk freezing transition entirely or only push it to weaker coupling, and is the improved action still reflection-positive / free of spurious phases?
7. **Resource trade-off.** Extra plaquette terms cost more gates/qudit operations; is there a sweet spot for the number of characters versus non-subgroup partitionings (Fibonacci [HA], magnetic basis [GB])?
8. **Group identity check.** Confirm whether the project's order-108 group is $\Sigma(36\times3)$ against `groups/S108`, and pin the precise freezing $\beta$ and transition order for $S(108),S(216),S(648)$ under the improved actions (only $S(1080)$, $\beta_f=3.935(5)$, is quoted in the sources).
9. **Self-duality as a check.** Does the $\beta_f\leftrightarrow\beta_a$ self-duality (Villain, monopoles suppressed) and its SU(N)/$\mathbb Z(N)$ generalization [GM] survive in multi-character improved actions, providing a consistency check on tuned couplings?
10. **Modern ML parametrizations.** Could ML-tuned multi-operator parametrizations [HO] yield better higher-order coefficient choices for the discrete setup than analytic Iwasaki/DBW2 values?

## 11. References

*(Bibliographic corrections from the verification stage have been applied. Confirmed and corrected references only; unverifiable leads are listed separately.)*

1. K. G. Wilson and J. Kogut (1974). *The renormalization group and the $\epsilon$ expansion.* Physics Reports **12**, 75–199. https://doi.org/10.1016/0370-1573(74)90023-4 — cited above as [W1].
2. P. Hasenfratz and F. Niedermayer (1994). *Perfect lattice action for asymptotically free theories.* Nucl. Phys. B **414**, 785–814. arXiv:hep-lat/9308004. https://arxiv.org/abs/hep-lat/9308004 — cited as [H1].
3. T. DeGrand, A. Hasenfratz, P. Hasenfratz, F. Niedermayer (1995). *The classically perfect fixed point action for SU(3) gauge theory.* Nucl. Phys. B **454**, 587–614. arXiv:hep-lat/9506030. https://arxiv.org/abs/hep-lat/9506030 — cited as [H2].
4. T. DeGrand, A. Hasenfratz, P. Hasenfratz, F. Niedermayer (1996). *Fixed point actions for SU(3) gauge theory.* Phys. Lett. B **365**, 233–238. arXiv:hep-lat/9508024. https://arxiv.org/abs/hep-lat/9508024 — cited as [H3]. *(Venue corrected: published in Phys. Lett. B 365, 1996; the "Nucl. Phys. B454, 615" attribution belongs to a sibling paper, hep-lat/9506030 = ref. 3.)*
5. W. Bietenholz and U.-J. Wiese (1996). *Perfect lattice actions for quarks and gluons.* Nucl. Phys. B **464**, 319–352. arXiv:hep-lat/9510026. https://arxiv.org/abs/hep-lat/9510026 — cited as [BW]. *(Page range corrected to 319–352.)*
6. P. Hasenfratz (1998). *Lattice QCD without tuning, mixing and current renormalization.* Nucl. Phys. B **525**, 401–409. arXiv:hep-lat/9802007. https://arxiv.org/abs/hep-lat/9802007 — cited as [H4]. *(Title corrected; the spurious phrase "mirror fermions and the Ginsparg–Wilson relation" removed.)*
7. P. Hasenfratz (1998). *Prospects for perfect actions.* Nucl. Phys. B Proc. Suppl. **63**, 53–58. arXiv:hep-lat/9709110. https://arxiv.org/abs/hep-lat/9709110 — companion proceedings to [H4].
8. P. Hasenfratz (1998). *Perfect actions: from the theoretical background to recent developments.* Prog. Theor. Phys. Suppl. **131**, 189–231. https://academic.oup.com/ptps/article/doi/10.1143/PTPS.131.189/1926843 — cited as [H5].
9. A. Alexandru, P. F. Bedaque, S. Harmalkar, H. Lamm, S. Lawrence, N. C. Warrington (NuQS) (2019). *Gluon Field Digitization for Quantum Computers.* Phys. Rev. D **100**, 114501. arXiv:1906.11213. https://arxiv.org/abs/1906.11213
10. C. Rebbi (1980). *Phase structure of non-Abelian lattice gauge theories.* Phys. Rev. D **21**, 3350. https://journals.aps.org/prd/abstract/10.1103/PhysRevD.21.3350 — cited as [RE].
11. D. Petcher and D. H. Weingarten (1980). *Monte Carlo calculations and a model of the phase structure for gauge theories on discrete subgroups of SU(2).* Phys. Rev. D **22**, 2465. https://journals.aps.org/prd/abstract/10.1103/PhysRevD.22.2465 — cited as [PW].
12. G. Bhanot and C. Rebbi (1981). *Monte Carlo simulations of lattice models with finite subgroups of SU(3) as gauge groups.* Phys. Rev. D **24**, 3319–3322. https://journals.aps.org/prd/abstract/10.1103/PhysRevD.24.3319 — cited as [BR].
13. G. Bhanot (1982). *SU(3) lattice gauge theory in 4 dimensions with a modified Wilson action.* Phys. Lett. B **108**, 337–340. https://doi.org/10.1016/0370-2693(82)91207-2 — cited as [BH].
14. H. Grosse and H. Kühnelt (1981). *Phase structure of lattice gauge theories for non-abelian subgroups of SU(3).* Phys. Lett. B **101**, 77–81. https://doi.org/10.1016/0370-2693(81)90500-4 — cited as [GK].
15. H. Flyvbjerg (1984). *Group space decimation: a way to simulate QCD by the 1080 element subgroup of SU(3)?* Nucl. Phys. B **243**, 350–364. https://www.sciencedirect.com/science/article/pii/0550321384900336 — cited as [FL].
16. T. Hartung, T. Jakobs, K. Jansen, J. Ostmeyer, C. Urbach (2022). *Digitising SU(2) gauge fields and the freezing transition.* Eur. Phys. J. C **82**, 237. arXiv:2201.09625. https://arxiv.org/abs/2201.09625 — cited as [HA]. *(EPJC article number is 237, not 192.)*
17. T. Jakobs, T. Hartung, K. Jansen, J. Ostmeyer, C. Urbach (2022). *Digitizing SU(2) gauge fields and what to look out for when doing so.* PoS LATTICE2022, 015. arXiv:2212.09496. https://arxiv.org/abs/2212.09496 — cited as [JA]. *(Binary T/O/I freeze at increasing $\beta_c$, order-120 latest; unweighted Fibonacci-lattice partitionings give the largest unfrozen $\beta$-range; Voronoi weighting applies to the linear/volleyball — not Fibonacci — partitionings.)*
18. (Freezing-obstruction statement) — see ref. 16 [HA].
19. K. Symanzik (1983). *Continuum limit and improved action in lattice theories: (I) Principles and $\phi^4$ theory; (II) O(N) non-linear $\sigma$ model in perturbation theory.* Nucl. Phys. B **226**, 187–204 and 205–227. https://doi.org/10.1016/0550-3213(83)90468-6 — cited as [S1].
20. M. Lüscher and P. Weisz (1985). *On-shell improved lattice gauge theories.* Commun. Math. Phys. **97**, 59–77; Erratum **98** (1985) 433. https://link.springer.com/article/10.1007/BF01206178 — cited as [LW].
21. B. Assi and H. Lamm (2024). *Digitization and subduction of SU(N) gauge theories.* Phys. Rev. D **110**, 074511. arXiv:2405.12204; FERMILAB-PUB-24-0177-T. https://arxiv.org/abs/2405.12204
22. (See ref. 9 for the $\beta_f=3.935(5)$ $S(1080)$ freezing result.) [9]
23. A. Alexandru, P. F. Bedaque, R. Brett, H. Lamm (2022). *The spectrum of qubitized QCD: glueballs in a S(1080) gauge theory.* Phys. Rev. D **105**, 114508. arXiv:2112.08482. https://arxiv.org/abs/2112.08482
24. (Mixed fundamental–adjoint single-extra-character action) — see ref. 13 [BH] and ref. 27 [BC].
25. Y. Ji, H. Lamm, S. Zhu (NuQS) (2020). *Gluon Field Digitization via Group Space Decimation for Quantum Computers.* Phys. Rev. D **102**, 114513. arXiv:2005.14221. https://arxiv.org/abs/2005.14221
26. Y. Ji, H. Lamm, S. Zhu (2023). *Gluon Digitization via Character Expansion for Quantum Computers.* Phys. Rev. D **107**, 114503. arXiv:2203.02330. https://arxiv.org/abs/2203.02330
27. G. Bhanot and M. Creutz (1981). *Variant actions and phase structure in lattice gauge theory.* Phys. Rev. D **24**, 3212–3217. https://journals.aps.org/prd/abstract/10.1103/PhysRevD.24.3212 — cited as [BC]. *(SU(2)-specific; limits SU(2)/SO(3)/$\mathbb Z(2)$. The generic SU(N) mixed fundamental–adjoint character form is a standard later generalization, not this paper's explicit content.)*
28. U. M. Heller (1995). *SU(3) lattice gauge theory in the fundamental–adjoint plane and scaling along the Wilson axis.* Phys. Lett. B **362**, 123–127. arXiv:hep-lat/9508009. https://arxiv.org/abs/hep-lat/9508009 — cited as [HE]. *(Companion: Nucl. Phys. B Proc. Suppl. 47, 262, arXiv:hep-lat/9509010.)*
29. M. Grady (2005). *Critical or tricritical point in mixed-action SU(2) lattice gauge theory?* Nucl. Phys. B **713**, 204–218. arXiv:hep-lat/0404015. https://arxiv.org/abs/hep-lat/0404015 — cited as [GR].
30. R. V. Gavai and M. Mathur (1999). *$Z_2$ monopoles, vortices and the universality of the SU(2) deconfinement transition.* Phys. Lett. B **458**, 331–337. arXiv:hep-lat/9905030. https://arxiv.org/abs/hep-lat/9905030 — cited as [GM].
31. S. Datta and R. V. Gavai (1997). *Stability of the bulk phase diagram of the SU(2) lattice gauge theory with fundamental–adjoint action.* Phys. Lett. B **392**, 172–176. arXiv:hep-lat/9610022. https://arxiv.org/abs/hep-lat/9610022 — cited as [DG].
32. G. P. Lepage and P. B. Mackenzie (1993). *On the viability of lattice perturbation theory.* Phys. Rev. D **48**, 2250–2264. arXiv:hep-lat/9209022. https://arxiv.org/abs/hep-lat/9209022 — cited as [LM].
33. V. Periwal (1996). *Improving lattice perturbation theory.* Phys. Rev. D **53**, 2605. arXiv:hep-lat/9507007. https://arxiv.org/abs/hep-lat/9507007 — cited as [PE].
34. Y. Iwasaki (1983, publ. 2011). *Renormalization group analysis of lattice theories and improved lattice action: II. Four-dimensional non-abelian SU(N) gauge model.* arXiv:1111.7054 (UTHEP-118, 1983). https://arxiv.org/abs/1111.7054 — cited as [IW].
35. QCD-TARO Collaboration (Ph. de Forcrand, T. Takaishi, et al.) (2000). *Renormalization group flow of SU(3) lattice gauge theory — numerical studies in a two coupling space.* Nucl. Phys. B **577**, 263–278. arXiv:hep-lat/9911033. https://arxiv.org/abs/hep-lat/9911033 — cited as [TARO].
36. A. Papa (1998). *Fixed point actions and on-shell tree-level Symanzik improvement.* Phys. Lett. B **437**, 123–130. arXiv:hep-lat/9803014. https://arxiv.org/abs/hep-lat/9803014 — cited as [PA].
37. W. Bietenholz (1998). *Perfect and quasi-perfect lattice actions.* arXiv:hep-lat/9802014. https://arxiv.org/abs/hep-lat/9802014 — cited as [QP].
38. J.-M. Drouffe and J.-B. Zuber (1983). *Strong coupling and mean field methods in lattice gauge theories.* Physics Reports **102**, 1–119. https://doi.org/10.1016/0370-1573(83)90034-0 — cited as [DZ].
39. R. Balian, J. M. Drouffe, C. Itzykson (1974). *Gauge fields on a lattice. I. General outlook.* Phys. Rev. D **10**, 3376–3395. https://link.aps.org/doi/10.1103/PhysRevD.10.3376 — cited as [BDI1].
40. R. Balian, J. M. Drouffe, C. Itzykson (1975). *Gauge fields on a lattice. III. Strong-coupling expansions and transition points.* Phys. Rev. D **11**, 2104–2119. https://journals.aps.org/prd/abstract/10.1103/PhysRevD.11.2104 — cited as [BDI3].
41. D. Tong (2018). *Lecture Notes on Gauge Theory, Section 4: Lattice Gauge Theory.* University of Cambridge DAMTP. https://www.damtp.cam.ac.uk/user/tong/gaugetheory/4lattice.pdf — cited as [TO].
42. P. Menotti and E. Onofri (1981). *The action of SU(N) lattice gauge theory in terms of the heat kernel on the group manifold.* Nucl. Phys. B **190** [FS3], 288–300. https://doi.org/10.1016/0550-3213(81)90560-5 (also https://cds.cern.ch/record/134274) — cited as [MO].
43. A. A. Migdal (1975). *Recursion equations in gauge field theories.* Sov. Phys. JETP **42**, 413 (Zh. Eksp. Teor. Fiz. **69**, 810). — cited as [MI].
44. W. Grimus and P. O. Ludl (2012). *Finite flavour groups of fermions.* J. Phys. A: Math. Theor. **45**, 233001. arXiv:1110.6376. https://arxiv.org/abs/1110.6376 — cited as [GL]. *(Corrected: vol. 45, art. 233001, 2012 — not "47, 203001, 2014".)*
45. P. O. Ludl (2011). *Comments on the classification of the finite subgroups of SU(3).* J. Phys. A **44**, 255204. arXiv:1101.2308. https://arxiv.org/abs/1101.2308 — cited as [LU].
46. P. O. Ludl (2017). *GAP listing of the finite subgroups of U(3) of order smaller than 2000.* arXiv:1702.00005. https://arxiv.org/abs/1702.00005 — cited as [LG].
47. S. Osorio Perez, E. M. Murairi, E. J. Gustafson, H. Lamm (2025). *Primitive Quantum Gates for an SU(3) Discrete Subgroup: $\Sigma(72\times3)$.* arXiv:2511.17437. https://arxiv.org/abs/2511.17437 — cited as [OP].
48. M. D'Elia, F. Farchioni, A. Papa (1995). *Scaling and topology in the 2-d O(3) $\sigma$-model on the lattice with the fixed point action.* Nucl. Phys. B **456**, 313–338. arXiv:hep-lat/9505004. https://arxiv.org/abs/hep-lat/9505004 — cited as [DP]. *(Corrected attribution: this title/arXiv id is by D'Elia–Farchioni–Papa, not Blatter et al.)*
49. M. Blatter, R. Burkhalter, P. Hasenfratz, F. Niedermayer (1996). *Instantons and the fixed point topological charge in the two-dimensional O(3) $\sigma$ model.* Phys. Rev. D **53**, 923. arXiv:hep-lat/9508028. https://arxiv.org/abs/hep-lat/9508028 — cited as [BI]. *(The Blatter et al. paper, distinct from ref. 48.)*
50. P. Hasenfratz, V. Laliena, F. Niedermayer (1998). *The index theorem in QCD with a finite cut-off.* arXiv:hep-lat/9801021. https://arxiv.org/abs/hep-lat/9801021 — supporting ref. for the FP/Ginsparg–Wilson connection.

### Modern parametrization (confirmed lead)

51. K. Holland, A. Ipp, D. I. Müller, U. Wenger (2025; rev. 2026). *Machine-learned renormalization-group-improved gauge actions and classically perfect gradient flows.* Phys. Rev. D (2026), DOI 10.1103/k41k-2pnc. arXiv:2504.15870. https://arxiv.org/abs/2504.15870 — cited as [HO]. *(Confirmed. Gauge-equivariant-CNN parametrization of the SU(3) fixed-point action together with its **classically perfect gradient flow** — flow free of tree-level cutoff effects to all orders in a; gradient-flow observables show <1% discretization error up to a≈0.14 fm. Continuous SU(3); no discrete subgroups or character expansions.)*

52. C. J. Garofalo, P. Bauer et al. (2024). *Digitizing lattice gauge theories in the magnetic basis: reducing the breaking of the fundamental commutation relations.* Eur. Phys. J. C (2024). arXiv:2311.11928. https://arxiv.org/abs/2311.11928 — cited as [GB]. *(Author list/venue as supplied by the survey; exact details to-be-confirmed.)*

### Unverified leads

The following were referenced in the gathered material but not independently citation-verified; treat as starting points, not established citations:
- "DeGrand, Hasenfratz, Hasenfratz, Niedermayer, *Towards a perfect fixed point action for SU(3) gauge theory*, arXiv:hep-lat/9412058" — quoted in the digitization survey as a perfect-action lead; not separately verified here. The verified SU(3) FP-action references are refs. 3 and 4 above; use those unless 9412058 is confirmed.
- The "Wikipedia / groupprops / Valentiner (1889)" *Valentiner group* identification (order-1080 $=3.A_6$, center $\mathbb Z_3$, SmallGroups [1080,260]) is standard group theory but cited only to encyclopedic sources [VA]; corroborate against ref. 44 [GL] (which gives the $\Sigma(360\times3)$ irreps) for a peer-reviewed source.
- Precise EPJC volume/page for ref. 52 and the exact PRL coordinates for ref. 51 were supplied by the survey but not re-verified.

## 12. How to extend the search

Precise queries and databases to go deeper:

- **INSPIRE-HEP** (`inspirehep.net`):
  - `title "fixed point action" and title SU(3)` ; `title "perfect action" and gauge`
  - `title digitization and (S(1080) or Sigma(360))` ; `title "discrete subgroup" and SU(3) and lattice`
  - `title "fundamental-adjoint" and (SU(3) or SU(2))` ; `title "bulk transition" and lattice gauge`
  - `refersto recid <id of Alexandru 2019 / Assi-Lamm 2024>` to crawl the forward-citation tree of the digitization line.
- **arXiv** (`arxiv.org`, hep-lat / quant-ph):
  - `cat:hep-lat AND (abs:"discrete subgroup" AND abs:"character") `
  - `cat:hep-lat AND abs:"Sigma(360" ` ; `abs:"heat kernel" AND abs:"finite group" AND lattice`
  - `cat:quant-ph AND abs:"gauge field digitization"` ; `abs:"Valentiner" AND lattice`
- **Google Scholar** "Cited by" on: Alexandru et al. 2019 [9]; Assi–Lamm 2024 [21]; Ji–Lamm–Zhu 2020/2023 [25,26]; Hartung et al. 2022 [16]; Bhanot–Creutz 1981 [27]; Hasenfratz–Niedermayer 1994 [H1].
- **Specific gaps to chase:**
  - "multi-representation" or "multiple character" lattice action + finite group + variational/optimization — to confirm the project's novelty claim.
  - "$q$-deformed" / "magnetic basis" / "Fibonacci" digitization for **SU(3)** specifically (most current work is SU(2)).
  - Character tables and conjugacy classes for $S(108)/\Sigma(36\times3)$ to settle the order-108 identity (cross-check ref. 46 [LG] GAP listing against `groups/S108`).
  - Confirm ref. 51 [HO] (ML RG-improved actions) and ref. 52 [GB] (magnetic basis) bibliographic details via INSPIRE/arXiv before citing as established.
