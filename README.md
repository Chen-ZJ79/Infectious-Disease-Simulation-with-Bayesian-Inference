# Bayesian SIR Parameter Estimation (MH + Gibbs + Adaptive MCMC)

This repository provides an R implementation for Bayesian parameter estimation in a **Susceptible–Infected–Removed (SIR)**-type epidemic model using **Metropolis–Hastings (MH)** and **Gibbs sampling**.

## Results

### 1️⃣ Standard Metropolis–Hastings Sampling
<img src="Rplot.pdf" alt="MH Trace Plots" width="600"/>

Each panel shows the parameter trajectory across iterations.
The red line marks the true simulated value.
This version uses a fixed proposal covariance, resulting in slower convergence.

### 2️⃣ Adaptive MH + Gibbs Sampling
<img src="Rplot01.pdf" alt="Adaptive MH Trace Plots" width="600"/>

Adaptive MH automatically updates the proposal covariance using the running chain covariance,
leading to smoother mixing and faster convergence.
Parameters stabilize around their true values more efficiently.