%% This file reproduceses the results of the Split Bregman on Shepp Logan
% 
% based on: The Split Bregman Method for L1-Regularized Problems
% T.Goldstein and S.Osher https://doi.org/10.1137/080725891
%
% The script is inspired also by the work of Benjamin Trémoulhéac
% Vincenzo Schiano
% Feb 2017

clc; clear all; close all;
addpath(genpath('utils_functions'))


%% Set up test image
X_true = phantom;

%Noise
[n,m] = size(X_true);
sigma=0.05;
X_noise = X_true + sigma.*randn(n,m);

%s=.0095; % noise level (NB actual Rician stdev depends on signal, see ricestat)
%X_noise = ricernd(X_true, s); %add rician noise




%% Compute and print

mu=1.9; %0.05;

lambda= 0.05; %mu*2;
Tol=5*10^(-5);



fprintf('Split Isotropic Algorithm by Definition \n')
[u0,k0,rel0,l0] = SplitB_IsotropicTV_Definition(X_noise,lambda,mu,Tol, 30, X_true);
Str0='SBI-TV-Def';
[outStr0,nEr0,psnr0]=printErPnsImages(Str0, u0, X_true);
fprintf(outStr0);

[u1,k1,rel1,l1] = SplitB_IsotropicTV_rev01(X_noise,lambda,mu,Tol, 50, X_true);
Str1='SBI-TV-rev1';
[outStr1,nEr1,psnr1]=printErPnsImages(Str1, u1, X_true);
fprintf(outStr1);

[u2,k2,rel2,l2] = SplitB_IsotropicTV_rev02(X_noise,lambda,mu,Tol, 50, X_true);
Str2='SBI-TV-rev2';
[outStr2,nEr2,psnr2]=printErPnsImages(Str2, u2, X_true);
fprintf(outStr2);



[u3,k3,rel3,l3] = SplitB_IsotropicTV(X_noise,lambda,mu,Tol, 50, X_true);
Str3='SBI-TV';
[outStr3,nEr3,psnr3]=printErPnsImages(Str3, u3, X_true);
fprintf(outStr3);



fprintf('Split Isotropic Algorithm by  Jacqueline Bush\n')
Tol=5*10^(-3);
[u,l,k] = SplitIsotropic2(X_noise,X_true,Tol,lambda,mu);
Str='SI2';
[outStr,nEr,psnr]=printErPnsImages(Str, u, X_true);
fprintf(outStr);



%% Display

figure('Name','Comparison among denoiser'); colormap gray;
subplot(241); imagesc(X_true); axis off image; title('Original');
subplot(242); imagesc(X_noise); axis off image; title('Noisy');
subplot(243); imagesc(u0); axis off image off; title(Str0);
subplot(244); imagesc(u1); axis off image; title(Str1);
subplot(245); imagesc(u2); axis off image; title(Str2);
subplot(246); imagesc(u3); axis off image; title(Str3);
subplot(247); imagesc(u); axis off image; title(Str);




%Fig 2
figure('Name','Error vs. iteration')
i=(0:k0); semilogy(i, l0(i+1)); hold on;
i=(0:k1); semilogy(i, l1(i+1))
i=(0:k2); semilogy(i, l2(i+1),'--')
i=(0:k3); semilogy(i, l3(i+1), '--')
i=(0:k); semilogy(i, l(i+1))
xlabel('Number of Iterations')
ylabel('Error')
legend(Str0, Str1, Str2, Str3, Str, 'Location', 'southeast')
legend('boxoff')


%Fig 3
figure('Name','Rel. Error vs. iteration')
i=(0:k0); semilogy(i, rel0(i+1)); hold on;
i=(0:k1); semilogy(i, rel1(i+1))
i=(0:k2); semilogy(i, rel2(i+1),'--')
i=(0:k3); semilogy(i, rel3(i+1), '--')
xlabel('Number of Iterations')
ylabel('Relative Error')
legend(Str0, Str1, Str2, Str3, 'Location', 'southeast')
legend('boxoff')


%Fig 4
figure('Name','Comparison of Rel. Err. and PDNR')
c = categorical({Str0, Str1, Str2, Str3, Str});
subplot(211);
errors = [nEr0 nEr1 nEr2 nEr3 nEr];
bar(c,errors,'BaseValue',min(errors)),set(gca,'yscale','log')
title('Err');
subplot(212);
Psnrs = [psnr0 psnr1 psnr2 psnr3 psnr];
bar(c,Psnrs,'BaseValue',min(Psnrs))
title('PSNR Err');










