# Bayesian Low-Tubal-Rank Robust Tensor Factorization with Multi-Rank Determination
Matlab source codes for Bayesian Tubal-Rank Tensor Factorization (BTRTF)

## Introduction

This is the PyTorch implementation of the [RotatE](https://openreview.net/forum?id=HkgEQnRqYQ) model for knowledge graph embedding (KGE). We provide a toolkit that gives state-of-the-art performance of several popular KGE models. The toolkit is quite efficient, which is able to train a large KGE model within a few hours on a single GPU.

%[Algorithms]%
The matlab codes provided here implement the BTRTF algorithm presented in the 
following paper:

## Usage
Image denoising with BTRTF on example images from the dataset:
```
Demo_BTRTF.m
```

**Using the library**
## Descriptions of the files in this repository 

 - Demo_BTRTF.m provides example usage of BTRTF for image denoising.
 - BTRTF.m implements BTRTF with multi-rank determination described as Alg.1 in [1].
 - myPSNR.m computes the PSNR value from the ground truth and the estimated low-rank component.

The python libarary is organized around 3 objects:

 - TrainDataset (dataloader.py): prepare data stream for training
 - TestDataSet (dataloader.py): prepare data stream for evluation
 - KGEModel (model.py): calculate triple score and provide train/test API

The run.py file contains the main function, which parses arguments, reads data, initilize the model and provides the training loop.

## Citation

If you find our codes helpful, please consider cite the following [paper](https://ieeexplore.ieee.org/abstract/document/8740980):
```
@article{
    zhou2019bayes,
    title={Bayesian Low-Tubal-Rank Robust Tensor Factorization with Multi-Rank Determination},
    author={Yang Zhou and Yiu-ming Cheung},
    journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
    year={2019},
    doi={10.1109/TPAMI.2019.2923240},
}
```
