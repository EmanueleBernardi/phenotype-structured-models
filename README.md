# Phenotype-structured integro-differential models

Numerical simulations for the paper:

> E. Bernardi, T. Lorenzi, A. Tosin  
> *Derivation and quasi-invariant asymptotics of phenotype-structured integro-differential models*  
> Preprint: [arXiv:2510.15646](https://arxiv.org/abs/2510.15646)

---

## Overview

This repository contains the MATLAB code used to produce the numerical simulations presented in the paper. Starting from a stochastic agent-based model for phenotype-structured biological populations, we derive a mesoscopic integro-differential equation (IDE) and, via a quasi-invariant limit, a non-local Fokker-Planck-type equation (PDE). The simulations compare the three levels of description across different values of the scaling parameter $\varepsilon$.

## Repository structure
```plaintext
../
│
├── src/
│   ├── parameters.m     % Set and save all simulation parameters
│   ├── sol_MC.m         % Monte Carlo simulation of the agent-based model
│   ├── sol_IDE.m        % Numerical solution of the IDE
│   ├── sol_PDE.m        % Numerical solution of the PDE
│   ├── binary.m         % One time step of the agent-based dynamics (called by sol_MC.m)
│   └── Plotting.m       % Generate and save figures
│
└── data/                % Created automatically by parameters.m
```

## Usage

All scripts are in the `src/` folder. Run them in the following order from MATLAB:

1. **`parameters.m`** — set all parameters and save them to `data/parameters.mat`. The `data/` folder is created automatically if it does not exist.
2. **`sol_MC.m`** — run the Monte Carlo simulation for three values of $\varepsilon$.
3. **`sol_IDE.m`** — solve the IDE for three values of $\varepsilon$.
4. **`sol_PDE.m`** — solve the PDE.
5. **`Plotting.m`** — generate and save the figure.

To reproduce the simulations for a different value of `alpha`, set the desired value in `parameters.m` and rerun all scripts. Output files are named with the corresponding `alpha` value to avoid overwriting.

---

## Parameters

The key parameters are set in `parameters.m`:

| Parameter | Description | Default value |
|-----------|-------------|---------------|
| `vm` | Fittest phenotypic trait | `1.5` |
| `R` | Radius of the mollifier $\psi$ | `5` |
| `delta` | Transition layer amplitude of $\psi$ | `0.5` |
| `alpha` | Drift coefficient | `-0.3` |
| `beta` | Diffusion coefficient | `0.4` |
| `Tf` | Final time | `10` |
| `N` | Number of agents (Monte Carlo) | `1e5` |
| `rho0` | Initial mass | `0.3` |
| `v_0` | Initial mean phenotype | `0` |
| `s0` | Inverse of initial variance | `2` |

The three values of $\varepsilon$ used in the simulations are $\varepsilon^2 \in \{1, 10^{-1}, 10^{-2}\}$.

---

## Reproducibility notes

**Choice of `alpha`.** In the paper, simulations are carried out for three values of the drift coefficient: `alpha = -0.3` (negative drift), `alpha = 0` (zero drift), and `alpha = 0.3` (positive drift). To reproduce the figures in the paper, run the full pipeline three times with these values. The output files are saved with the corresponding `alpha` value in the filename to avoid overwriting.

**Computational cost.** The Monte Carlo simulation (`sol_MC.m`) with `N = 1e5` or more agents is computationally intensive and may take several hours on a standard laptop. We recommend running it on a performant workstation.

---

## Requirements

- MATLAB R2020a or later
- No additional toolboxes required

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Citation

If you use this code, please cite:

```bibtex
@misc{bernardi2025phenotype,
  author    = {Bernardi, Emanuele and Lorenzi, Tommaso and Tosin, Andrea},
  title     = {Derivation and quasi-invariant asymptotics of phenotype-structured integro-differential models},
  year      = {2025},
  eprint    = {2510.15646},
  archivePrefix = {arXiv},
  primaryClass  = {math.AP}
}
```

---

## Contact

Emanuele Bernardi — [emanuele.bernardi@polito.it](mailto:emanuele.bernardi@polito.it)  
Personal website: [emanuelebernardi.github.io/personal-site](https://emanuelebernardi.github.io/personal-site/)