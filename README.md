# Bayesian-Inference-and-Data-Assimilation-Using-Kalman-Filtering
This project implements Bayesian data assimilation in MATLAB to analyse 90 days of agricultural commodity log-price and rainfall data. A linear Gaussian state-space model links observed log-prices to an unobserved environmental condition, which evolves according to daily rainfall and the previous day’s condition.

The environmental state is combined with two unknown parameters—the baseline log-price and the rainfall effect—into an augmented state vector. A Kalman filter recursively computes forecast and analysis means and covariance matrices, updating estimates and their uncertainty as new observations are assimilated.

Results include trajectories of the filtered state and parameter means with pointwise 95% credible bands, posterior means and standard deviations, and 95% credible intervals for both parameters. A joint posterior credible ellipse visualises their dependence, while the posterior distribution of a derived quantity, θ = μ − 10β, demonstrates how uncertainty is propagated through a linear combination of the parameters.

## Project Report
[View full report](DataAssimilation_2.pdf)

## License
This project is provided for viewing and evaluation purposes only.
Reuse, modification, or distribution is not permitted without explicit permission.
