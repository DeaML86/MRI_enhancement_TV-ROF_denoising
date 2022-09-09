%% This file reproduceses the results on the Split Bregman method
% 
% The Split Bregman Method for L1-Regularized Problems
% T.Goldstein and S.Osher https://doi.org/10.1137/080725891
%
% The script is inspired also by the work of Benjamin Trémoulhéac
% Vincenzo Schiano
% Nov 2017

clc; clear all; close all;
addpath(genpath('utils_functions'))



%% Set up test image and noise
%f = double(phantom); 
f_true = double(imread('Lena512.png'));
[n,m] = size(f_true);
sigma=15;
g = f_true + randn(m,n)*sigma;


%% Compute and print

mu=0.05;
lambda=0.1;
Tol=5*10^(-7);

[g_denoise_sbi,k1,rel1,l1] = SplitB_IsotropicTV(g,lambda,mu,Tol, 50, f_true);
Str1='SBI-TV';
[outStr1,nEr1,psnr1]=printErPnsImages(Str1, g_denoise_sbi, f_true);
fprintf(outStr1);

[g_denoise_si2,l2,k2] = SplitIsotropic2(g,f_true,Tol,lambda,mu);
Str2='SI2';
[outStr2,nEr2,psnr2]=printErPnsImages(Str2, g_denoise_si2, f_true);
fprintf(outStr2);


Theta=0.2; % 0.2- 10 - 20
g_denoise_rof = ROFdenoise(g,Theta);
StrROF='ROFdenoise';
[outStrROF,nErROF,psnrROF]=printErPnsImages(StrROF, g_denoise_rof, f_true);
fprintf(outStrROF);


%mu = 20;
fprintf('Split Isotropic Algorithm \n')
g_denoise_atv = SB_ATV(g,mu);
ATV='ATV';
[outStrATV,nErATV,psnrATV]=printErPnsImages(ATV, reshape(g_denoise_atv,n,m), f_true);
fprintf(outStrATV);

g_denoise_itv = SB_ITV(g,mu);
ITV='ITV';
[outStrITV,nErITV,psnrITV]=printErPnsImages(ITV, reshape(g_denoise_itv,n,m), f_true);
fprintf(outStrITV);

Tol=5*10^(-3); mu = 20;
beta=2.5;%beta=0.005;
[g_denoise_si2lse,l3,k3] = SplitIsotropic2Lse(g,f_true,Tol,lambda,mu,beta);
Str3='SI2Lse';
[outStr3,nEr3,psnr3]=printErPnsImages(Str3, g_denoise_si2lse, f_true);
fprintf(outStr3);




niter=100; lambda= 1.3; % 1 - 1.7 
g_denoise_tvl1 = TVL1denoise(g,lambda,niter);
StrTV='TVL1denoise';
[outStrTV,nErTV,psnrTV]=printErPnsImages(StrTV, g_denoise_tvl1, f_true);
fprintf(outStrTV);








%% Display

%Fig 1
figure('Name','Comparison between denoisers'); colormap gray;
subplot(331); imagesc(f_true); axis off image; title('Ground truth');
subplot(332); imagesc(g); axis off image; title('Noisy');
subplot(333); imagesc(reshape(g_denoise_atv,n,m)); axis image off; title(ATV);
subplot(334); imagesc(reshape(g_denoise_itv,n,m)); axis image; title(ITV);
subplot(335); imagesc(g_denoise_sbi); axis image; title(Str1);
subplot(336); imagesc(g_denoise_si2); axis image; title(Str2);
subplot(337); imagesc(g_denoise_si2lse); axis off image; title(Str3);
subplot(338); imagesc(g_denoise_tvl1); axis image; title(StrTV);
subplot(339); imagesc(g_denoise_rof); axis image; title(StrROF);


%Fig 2
figure('Name','Error vs. iteration')
i=(0:k1); semilogy(i, l1(i+1)); hold on;
i=(0:k2); semilogy(i, l2(i+1))
i=(0:k3); semilogy(i, l3(i+1))

xlabel('Number of Iterations')
ylabel('Error')
legend(Str1, Str2, Str3, 'Location', 'southeast')
legend('boxoff')



figure('Name','Comparison of Rel. Err. and PDNR')
c = categorical({ATV, ITV, Str1, Str2, Str3, StrTV, StrROF});
subplot(211);
errors = [nErATV nErITV nEr1 nEr2 nEr3 nErTV nErROF];
bar(c,errors,'BaseValue',min(errors))
title('Err');
subplot(212);
Psnrs = [psnrATV psnrITV psnr1 psnr2 psnr3 psnrTV psnrROF];
bar(c,Psnrs,'BaseValue',min(Psnrs))
title('PSNR Err');







