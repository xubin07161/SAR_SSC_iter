clear all,close all,clc;

load Monarch(ENL=1).mat;
ENL=1;
tic;
E_Img = SAR_SSC_iter(N_Img,'ENL',ENL);
toc;
figure,imshow(sqrt([N_Img,E_Img]),[0 1]);