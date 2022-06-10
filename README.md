# Cross Band Prediction

Cross band prediction Master Thesis

## Abstract 

The wireless communications technologies have experienced a exponential development
during the last decades. 5G is a prominent exponent whose one of its crucial
components is the Massive MIMO technology. By supporting multiple streams of signals
it allows a revamped signal reconstruction in terms of mobile traffic size, data rate,
latency and reliability. In this thesis work, we isolated this technology into a SIMO
approach (Single-Input Multiple-Output) to explore a Machine Learning modeling to address
the so-called Channel Prediction problem. Generally, the algorithms available to
perform Channel Estimation in FDD and TDD deployments incur computational complexity
downsides and require explicit feedback from client devices, which is typically
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
0.0584 and 0.0247.





``` bash
├───AR Modeling
│   └───Approach_1_SB
├───Baseline
│   └───PPCA
│       ├───Baseline_90TTI
│       │   ├───Cartesian
│       │   ├───Polar
│       │   └───Statistical Modeling
│       ├───Baseline_joint
│       ├───Baseline_UE1
│       │   └───PPCA-Holdout_UE1_files
│       ├───Baseline_UE2
│       └───Dummy_Baseline
├───Data
│   ├───data2UEjoint
│   ├───data_90TTI
│   │   └───2ndpolarisation
│   ├───data_UE1_600
│   │   └───2ndpolarisation
│   └───data_UE2_600
├───Sketches VAE
├───Stats_Figures
└───VAE
    ├───Approach_1_SB
    ├───Dummyex
    └───VED-Polar
   ```
