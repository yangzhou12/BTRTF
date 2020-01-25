# BTRTF
Matlab source codes of the Bayesian Tubal-Rank Tensor Factorization (BTRTF) algorithm presented in the paper [Bayesian Low-Tubal-Rank Robust Tensor Factorization with Multi-Rank Determination](https://ieeexplore.ieee.org/abstract/document/8740980).

## Usage
Image denoising with BTRTF on example images from the BSD500 dataset:
```
Demo_BTRTF.m
```

## Descriptions of the files in this repository 
 - *Demo_BTRTF.m* provides example usage of BTRTF for image denoising.
 - *BTRTF.m* implements BTRTF with multi-rank determination described as Alg.1 in [paper](https://ieeexplore.ieee.org/abstract/document/8740980).
 - *myPSNR.m* computes the PSNR value from the ground truth and the estimated low-rank component.

## Citation
If you find our codes helpful, please consider cite the following [paper](https://ieeexplore.ieee.org/abstract/document/8740980):
```
@article{
    zhou2019BTRTF,
    title={Bayesian Low-Tubal-Rank Robust Tensor Factorization with Multi-Rank Determination},
    author={Yang Zhou and Yiu-ming Cheung},
    journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
    year={2019},
    doi={10.1109/TPAMI.2019.2923240},
}
```
