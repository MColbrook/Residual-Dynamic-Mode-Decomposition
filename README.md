# Residual-Dynamic-Mode-Decomposition

Computation of spectral properties of Koopman operators associated with discrete-time autonomous dynamical systems. Highlights include: verified computation of spectra with error control (and avoidance of spurious modes) and computation of spectral measures with explicit high-order convergence.

This repository will grow as further papers are written and if you are interested in collaborating, please get in touch at: m.colbrook@damtp.cam.ac.uk

The code includes **"main_routines"** that are used across the papers (see that subfolder for the additional README file). Each paper has a gallery of examples.

To get started on a simple example, try the file **"Duffing_example.m"**

For a pedagogical review of DMD methods, including code, see the repository: https://github.com/MColbrook/DMD-Multiverse

For measure-preserving EDMD, see the repository: https://github.com/MColbrook/Measure-preserving-Extended-Dynamic-Mode-Decomposition

Code for the papers:

1. M.J. Colbrook, A. Townsend, *"Rigorous data-driven computation of spectral properties of Koopman operators for dynamical systems"* in **"Examples_gallery_1"**. Paper can be found here: http://www.damtp.cam.ac.uk/user/mjc249/pdfs/RigorousKoopman.pdf<br>
Please cite using the following bibtex: @article{colbrook2024rigorous,
  title={Rigorous data-driven computation of spectral properties of Koopman operators for dynamical systems},
  author={Colbrook, Matthew J and Townsend, Alex},
  journal={Communications on Pure and Applied Mathematics},
  volume={77},
  number={1},
  pages={221--283},
  year={2024},
  publisher={Wiley Online Library}
}

2. M.J. Colbrook, L. Ayton, M. Szőke, *"Residual Dynamic Mode Decomposition: Robust and verified Koopmanism"* in **"Examples_gallery_2"**. Paper can be found here: https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/residual-dynamic-mode-decomposition-robust-and-verified-koopmanism<br>
Please cite using the following bibtex: @article{colbrook2023residual,
  title={Residual dynamic mode decomposition: robust and verified {K}oopmanism},
  author={Colbrook, Matthew J and Ayton, Lorna J and Sz{\H{o}}ke, M{\'a}t{\'e}},
  journal={Journal of Fluid Mechanics},
  volume={955},
  pages={A21},
  year={2023},
  publisher={Cambridge University Press}
}

3. M.J. Colbrook, Q. Li, R.V. Raut, A. Townsend, *"Beyond expectations: Residual Dynamic Mode Decomposition and
Variance for Stochastic Dynamical Systems"* in **"Examples_gallery_3"**. Paper can be found here: [https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/residual-dynamic-mode-decomposition-robust-and-verified-koopmanism](https://link.springer.com/article/10.1007/s11071-023-09135-w)<br>
Please cite using the following bibtex: @article{colbrook2024beyond,
  title={Beyond expectations: residual dynamic mode decomposition and variance for stochastic dynamical systems},
  author={Colbrook, Matthew J and Li, Qin and Raut, Ryan V and Townsend, Alex},
  journal={Nonlinear Dynamics},
  volume={112},
  number={3},
  pages={2037--2061},
  year={2024},
  publisher={Springer}
}

4. M.J. Colbrook, *"Another look at Residual Dynamic Mode Decomposition with fewer Snapshots than Dictionary Size"* in **"Examples_gallery_4"**.

**Datasets** (needed for some of the examples) can be found here: https://www.dropbox.com/sh/xj59e5in7dfsobi/AAAfkxqa1x9WFSTgrvqoqqRqa?dl=0

Some of the code for setting up the examples makes use of Chebfun, which can be found at https://www.chebfun.org/.
