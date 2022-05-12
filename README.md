# Generative Adversarial Networks for the standardization of peripheral blood cell images

## Abstract

The wireless communications technologies have experienced a exponential development during the last decades. 5G is a prominent exponent whose one of its crucial
components is the Massive MIMO technology. By supporting multiple streams of signals it allows a revamped signal reconstruction in terms of mobile traffic size, data rate,
latency and reliability. In this thesis work, we isolated this technology into a SIMO
approach (Single-Input Multiple-Output) to explore a Machine Learning modeling to address the so-called Channel Prediction problem. Generally, the algorithms available to
perform Channel Estimation in FDD and TDD deployments incur computational complexity downsides and require explicit feedback from client devices, which is typically
prohibitive.
This thesis work focuses on Channel Prediction by aims of employing Machine and
Deep Learning models in order to reduce the computational complexity by further relying
in statistical modeling/learning. We explored the cross-Frequency Subband prediction
intra-TTI (Transmission Time Interval) by means of proposing 3 three models. These
intended to leverage frequency Multipath Components correlations along TTIs. The
first two ones are Probabilistic Principal Components Analysis (PPCA) and its Bayesian
approach, Bayesian Principal Components Analysis (BPCA). Then, we implemented a
Deep Learning Variational Encoder-Decoder (VED) architecture. These three models
intended to deal with the hugely high-dimensional space of the 4 datasets used by its
intrinsic dimensionality reduction.
By evaluating the performance of all models, we highlight that the PPCA model
outperformed the others in all the datasets, yielding a MSE error of 0.0029, 0.0112,
0.0584 and 0.0247

´´´bash
├───Baseline_90TTI
├───Baseline_joint
├───Baseline_UE1
├───Baseline_UE2
└───VED

´´´
