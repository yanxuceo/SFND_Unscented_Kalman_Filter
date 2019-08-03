# Unscented Kalman Filter Introduction

## Background
###
Recall that the extended Kalman filter uses the Jacobian matrix to linearize non-linear functions. The UKF, on the other hand, does not need to linearize non-linear functions; instead, the UKF takes representative points from a Gaussian distribution. These points wil be plugged into the non-linear equations as you'll see in the lectures.

<img src="screenshots/basic principle.png" width="460" height="280" />

