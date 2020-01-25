# BTRTF
Matlab Code for Bayesian Tubal-Rank Tensor Factorization (BTRTF)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Matlab source codes for                              %
%        Bayesian Low-Tubal-Rank Robust Tensor Factorization (BTRTF)          %
%                                                                             %
% Author: 		Yang ZHOU                                                     %
% Email: 		youngzhou12@gmail.com                                         %
% Release date: June 22, 2019                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[Algorithms]%
The matlab codes provided here implement the BTRTF algorithm presented in the 
following paper:

    Yang Zhou, Yiu-ming Cheung. 
    Bayesian Low-Tubal-Rank Robust Tensor Factorization with Multi-Rank Determination. 
    IEEE Transactions on Pattern Analysis and Machine Intelligence, 
    DOI:10.1109/TPAMI.2019.2923240, 2019.
---------------------------

%[Usages]%
Please refer to the file "Demo_BTRTF.m", which provides example usage of 
BTRTF for image denoising
---------------------------

%[Descriptions of the files in this package]%
1. Demo_BTRTF.m provides example usage of BTRTF for image denoising.
2. BTRTF.m implements BTRTF with multi-rank determination described as Alg.1 in [1].
3. myPSNR.m computes the PSNR value from the ground truth and the estimated low-rank component.
---------------------------

%[Restriction]%
In all documents and papers reporting research work that uses the matlab codes 
provided here, the respective author(s) must reference the following paper: 

[1] Yang Zhou, Yiu-ming Cheung. 
    Bayesian Low-Tubal-Rank Robust Tensor Factorization with Multi-Rank Determination. 
    IEEE Transactions on Pattern Analysis and Machine Intelligence, 
    DOI:10.1109/TPAMI.2019.2923240, 2019.
---------------------------
