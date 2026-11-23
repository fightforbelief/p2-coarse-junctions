# How many states does a coarse junction need?

Yunjie Guo and Hanqing Yang. Independent Researcher, Suzhou, China.

A pulse enters a junction through one of several ports. A model with too few internal states can predict more transfer to other ports than the measured dynamics allow. This project uses endpoint measurements at two observation lags to lower-bound the number of states required *inside the junction*, even when the surrounding cells have arbitrary hidden dynamics.

[Read the manuscript (PDF)](paper/main.pdf) · [LaTeX source](paper/main.tex) · [Diffusion protocol](protocol.md) · [Measured-trajectory protocol](empirical/protocol.md) · [Reviewer assessment](REVIEW.md)

## Results

For a time-homogeneous reversible Markov kernel with a deterministic observation map, the paper derives the exact minimum weighted cross-port transfer compatible with the one-step masses, junction mass, and state budget, when pair costs form a semimetric. The proof gives an attaining discrete-time model. It also gives the smallest observable fitting error when a surrogate is allowed to change its one-step predictions. Continuous-time models obey the finite-lag bound. A separate result gives a sharp infimum for their initial curvature.

For four equally loaded ports with one-step mass 0.2 per port and junction mass one, three internal states force a two-step cross-port mass of at least 1/15. An observation of 0.05 forces any three-state model to miss at least one of the five mass measurements by approximately 0.00902.

A reflecting diffusion benchmark covers 84 combinations of lag, resolution, port geometry and anisotropy. With separated ports the observations certify at least four junction states in 19 of 42 cases. With full-width ports meeting at corners the same all-pair score produces only the trivial lower bound in all 42 cases. A continuum heat-kernel calculation explains the geometric difference. In the designated independent-sampling experiment, 100,000 paths per initial ensemble recovered the four-state lower bound in 200 of 200 repeats. These counts describe a lower bound on model size. The underlying continuum has infinitely many states.

Two public recordings supply a separate measurement check. The passive C-Trap recording gave a four-state plug-in comparison at 12.8 microseconds and a trivial count of one at five longer lags. The static harmonic feedback trace gave two at the four shortest lags. Using complete excursions, the homogeneous one-state null was rejected for both records (Holm-adjusted p = 0.0002 each). The confidence statement for full-path tests uses observations beyond two-lag pair counts. [Inputs and replication](empirical/README.md)

## Using the bound

The function takes stationary **joint masses**. To turn a conditional probability into a joint mass, multiply by the mass of the source region. Count each unordered port pair once at the double lag, with the same mass scale and cost matrix.

```python
from junction import envelope

floor, groups, hub_masses = envelope([0.2] * 4, total=1, k=3)
print(floor)  # 0.06666666666666668, up to roundoff
```

Pass a symmetric nonnegative distance matrix as `weights` for a custom finite-lag score. `curvature=True` accepts general symmetric nonnegative pair weights with zero diagonal. Supply interface fluxes for the sharp continuous-time curvature calculation, or one-step joint masses for the relaxed finite-lag bound. At curvature level some reported masses can be zero, indicating a limiting allocation. Exact enumeration targets junctions with a small number of ports. Denser graphs may require other optimization methods.

## Reproducing the study

Two runs on the designated Windows research host produced byte-identical diffusion and sampling tables. The review run passed eight original and five trajectory tests. The original diffusion experiments use Python 3.12.10 and its standard library. From the project root on Windows, create the project environment if needed and run the audit, tests and benchmark.

```powershell
py -3.12 -m venv .venv
.\execute.ps1
.\.venv\Scripts\python.exe verify.py
```

On other platforms, with Python 3.12 or later, the corresponding entry points are `python -m unittest discover -v -s tests`, `python benchmark.py`, and `python verify.py`. The Monte Carlo experiment uses the Python 3.12 binomial sampler and seed 20261121. For a byte-level comparison with the archived output, run in a separate checkout and compare the two CSV hashes against the archived manifest. Floating-point libraries and interpreter versions can affect the last digits.

| File | Contents |
| --- | --- |
| [diffusion.csv](results/diffusion.csv) | All population results. Columns include observed transfer, full-coarse Chapman–Kolmogorov prediction, sharp and relaxed floors for each state budget, and the resulting lower bounds. |
| [sampling.csv](results/sampling.csv) | Counts, confidence limits and state-count lower bounds for each of the 1,200 simulated repetitions. Each repetition uses two independent ensembles. |
| [run.json](results/run.json) | Environment and SHA-256 digests of the experiment inputs and outputs. |
| [junction.py](junction.py) | Port partitions, optimal mass allocation and balanced-port formula. |
| [tests](tests/) | Closed-form, fractional-allocation, attaining-kernel, fitting-distance, graph-coloring, input-domain, confidence-limit and independent heat-solver checks. |

The paper figure reads the archived diffusion table. Rebuild it with Matplotlib 3.10.8 via `python paper/figure.py`, then run `tectonic paper/main.tex`. The figure is already included for builds that only need TeX.

## Research scope

The state-count bound targets time-homogeneous reversible dynamics observed through disjoint labels. A connected aperiodic kernel attains the sharp value when there is spare port capacity and unmonitored exterior mass. Continuous-time embeddability remains a stricter constraint. The two-lag confidence guarantee uses independent preparations and known region masses. The separate measured-trajectory tests use the chronology of visits and explicitly allow serial dependence. Fredkin and Rice derive aggregate state-count bounds from dwell times. Noé and colleagues model observed kinetics as a projection of full dynamics. Duan and colleagues estimate global aggregation structure through nonnegative rank. The result here fixes one hub, uses two-lag endpoint marginals, optimizes the scalar transfer over all admissible reversible kernels, and constructs a kernel attaining the optimum. This specific extremal result and its fitting-error consequence support submission as an applied-probability or model-reduction methods paper. The diffusion examples illustrate the theorem and delimit the measurement geometries in which it has power. [Fredkin and Rice](https://doi.org/10.6028/jres.090.053) · [Noé et al.](https://arxiv.org/abs/1309.3220) · [Duan et al.](https://doi.org/10.1137/18M1220790)

Correspondence may be directed to Hanqing Yang at hanqing.yang189@gmail.com.
