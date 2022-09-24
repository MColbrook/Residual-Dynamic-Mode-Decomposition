Here is a high-level summary of each main routine and what they do. Full details are given in the above .m files.

**ErgodicMoments:** Given a vector X of a signal and a non-negative integer N, this function computes the autocorrelations of X, or moments of the Koopman operator, from -N to N. This code uses the FFT for rapid computation.

**IsomMeas:** This function computes smoothed spectral measures of an isometry using the ResDMD matrices. These ResDMD matrices are "G" corresponding to <\psi_k,\psi_j>, "A" corresponding to <K\psi_k,\psi_j>, "L" corresponding to <K\psi_k,K\psi_j>. Details of how these three matrices are constructed from the data are given in the ResDMD paper, and can be seen in the various examples given (i.e., G=PSI_X'WPSI_X, A=PSI_X'WPSI_Y, L=PSI_Y'WPSI_Y). IsomMeas also takes as input the vector f, smoothing parameter epsilon, and angles theta\in[-pi,pi].

**KoopPseudoSpec:** Given the ResDMD matrices G, A, and L, and a grid of points z_pts, this function computes the minimal residuals for spectral parameters over the grid. This is used to compute pseudospectra and spectra of Koopman operators.

**MomentMeas:** Given the autocorrelations (which could be computed by ErgodicMoments, for example), this function computes a smoothed approximation of the associated spectra measure in the form of a Fourier series (represented using chebfun). The user can also specify the choice of filter function used.

**kernel_ResDMD:** This implements the choice of basis that uses kernel EDMD. The inputs are X and Y snapshots to form the dictionary (Xa and Ya) and X and Y snapshots to form the ResDMD matrices (Xb and Yb).
