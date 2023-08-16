# Preferential Sampling in Environmental Networks

This package present a general framework for detecting the preferential
sampling in an environmental process across space and/or time.

This is based on the work [A general theory for preferential sampling in
environmental
networks](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-13/issue-4/A-general-theory-for-preferential-sampling-in-environmental-networks/10.1214/19-AOAS1288.full)
by Joe Watson, James V. Zidek and Gavin Shaddick. This is achieved by
considering the joint distribution of an environmental process with a
site-selection process that considers where and when sites are placed to
measure the process. The environmental process may be spatial, temporal
or spatio-temporal in nature. By sharing random effects between the two
processes, the joint model is able to establish whether site placement
was stochastically dependent of the environmental process under study.

## The Spatio-temporal models

The data are defined by a process indexed by space and time

$$
Y(s, t) = \{y(s, t), (s, t) \in \mathcal{D} \subset \mathbb{R}^2\times \mathbb{R}\}
$$

and are observed at $n$ spatial locations or areas and at $T$ time
points.

When considering the spatio-temporal data, we usually need to define a
nonseparable spatio-temporal covariance function
$cov(y(s_i, t), y(s_j, u))$. In practice, to overcome the computational
complexity of the nonseparable models, we can simply assume separability
so that the space-time covariance function is decomposed into the
product of a purely spatial and a purely temporal term.
