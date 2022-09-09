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
%f_true = double(phantom); 
%f_true = double(imread('pout.tif'));
f_true = double(imread('Lena512.png'));
[n,m] = size(f_true);
sigma=15;
g = f_true + randn(n,m)*sigma;


%% Compute and print

mu=0.05;

lambda=0.1;
Tol=5*10^(-9);

[g_denoise_sbiDef,k0,rel0,l0] = SplitB_IsotropicTV_Definition(g,lambda,mu,Tol, 30, f_true);
Str0='SBI-TV-Def';
[outStr0,nEr0,psnr0]=printErPnsImages(Str0, g_denoise_sbiDef, f_true);
fprintf(outStr0);

[g_denoise_sbiRev1,k1,rel1,l1] = SplitB_IsotropicTV_rev02(g,lambda,mu,Tol, 50, f_true);
Str1='SBI-TV-rev2';
[outStr1,nEr1,psnr1]=printErPnsImages(Str1, g_denoise_sbiRev1, f_true);
fprintf(outStr1);


[g_denoise_sbi,k3,rel3,l3] = SplitB_IsotropicTV(g,lambda,mu,Tol, 50, f_true);
Str3='SBI-TV';
[outStr3,nEr3,psnr3]=printErPnsImages(Str3, g_denoise_sbi, f_true);
fprintf(outStr3);


Tol=5*10^(-7);
[g_denoise_si2,l2,k2] = SplitIsotropic2(g,f_true,Tol,lambda,mu);
Str2='SI2';
[outStr2,nEr2,psnr2]=printErPnsImages(Str2, g_denoise_si2, f_true);
fprintf(outStr2);








%% Display

%Fig 1
figure('Name','Comparison between denoisers'); colormap gray;
subplot(231); imagesc(f_true); axis off image; title('Ground truth');
subplot(232); imagesc(g); axis off image; title('Noisy');
subplot(233); imagesc(g_denoise_sbiDef); axis image; title(Str0);
subplot(234); imagesc(g_denoise_sbiRev1); axis image; title(Str1);
subplot(235); imagesc(g_denoise_si2); axis image; title(Str2);
subplot(236); imagesc(g_denoise_sbi); axis image; title(Str3);



%Fig 2
figure('Name','Error vs. iteration')
i=(0:k0); semilogy(i, l0(i+1)); hold on;
i=(0:k1); semilogy(i, l1(i+1))
i=(0:k2); semilogy(i, l2(i+1))
i=(0:k3); semilogy(i, l3(i+1))
xlabel('Number of Iterations')
ylabel('Error')
legend(Str0, Str1, Str2, Str3, 'Location', 'southeast')
legend('boxoff')



%Fig 3
figure('Name','Rel. Error vs. iteration')
i=(0:k0); semilogy(i, rel0(i+1)); hold on;
i=(0:k1); semilogy(i, rel1(i+1))
i=(0:k3); semilogy(i, rel3(i+1))
xlabel('Number of Iterations')
ylabel('Relative Error')
legend(Str0, Str1, Str3, 'Location', 'southeast')
legend('boxoff')


%Fig 4
figure('Name','Comparison of Rel. Err. and PDNR')
c = categorical({Str0, Str1, Str2, Str3});
subplot(211);
errors = [nEr0 nEr1 nEr2 nEr3];
bar(c,errors,'BaseValue',min(errors))
title('Err');
subplot(212);
Psnrs = [psnr0 psnr1 psnr2 psnr3];
bar(c,Psnrs,'BaseValue',min(Psnrs))
title('PSNR Err');







