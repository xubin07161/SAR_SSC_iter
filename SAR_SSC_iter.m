%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  SSCiter is an algorithm for SAR image despeckling. 
%  This algorithm reproduces the results from the article:
%  [1] B. Xu et al. 'An Iterative SAR Image Filtering Method Using
%       Nonlocal Sparse Model'
%  Please refer to this paper for a more detailed description of the algorithm.
%
%  BASIC USAGE EXAMPLES:
%
%     1) Using the default parameters
% 
%      img_filtered = SAR_SSC_iter(img,'ENL',ENL)
% 
%  INPUT ARGUMENTS (OPTIONAL):
%
%     1) img : The input SAR image should be intensity image.
%
%     2) ENL : The equivalent number of looks. The ENL can be obtained by 
%              by supervised or unsupervised estimation. For a homogeneous 
%              region, the ENL can be calculated by ENL=(mean)^2/var
%
%  OUTPUTS:
%     1) img_filtered  : The filtered intensity image                                             
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2015 Bin Xu.
% All rights reserved.
% This work should only be used for nonprofit purposes.
%
% AUTHORS:
%     Bin Xu, email: xubin07161@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%