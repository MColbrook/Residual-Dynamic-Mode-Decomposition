# Residual-Dynamic-Mode-Decomposition

Computation of spectral properties of Koopman operators associated with discrete-time autonomous dynamical systems. Highlights include: verified computation of spectra with error control (and avoidance of spurious modes) and computation of spectral measures with explicit high-order convergence.

This repository will grow as further papers are written and if you are interested in collaborating, please get in touch at: m.colbrook@damtp.cam.ac.uk

Code for the papers:

1. M.J. Colbrook, A. Townsend, *"Rigorous data-driven computation of spectral properties of Koopman operators for dynamical systems"* in **"Examples_gallery_1"**. Paper can be found here: http://www.damtp.cam.ac.uk/user/mjc249/pdfs/RigorousKoopman.pdf<br>
Please cite using the following bibtex: @article{colbrook2021rigorous,
  title={Rigorous data-driven computation of spectral properties of {K}oopman operators for dynamical systems},
  author={Colbrook, Matthew J. and Townsend, Alex},
  journal={arXiv preprint arXiv:2111.14889},
  year={2021}
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
Variance for Stochastic Dynamical Systems"* in **"Examples_gallery_3"**.

The code includes **"main_routines"** that are used across the papers (see that subfolder for the additional README file). Each paper has a gallery of examples.

To get started on a simple example, try the file Duffing_example.m

**Datasets** (needed for some of the examples) can be found here: https://www.dropbox.com/sh/xj59e5in7dfsobi/AAAfkxqa1x9WFSTgrvqoqqRqa?dl=0

Some of the code for setting up the examples makes use of Chebfun, which can be found at https://www.chebfun.org/.
