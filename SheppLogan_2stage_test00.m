% Simulation of shepp logan phantom
% Demo 5
clear all; close all; clc;
addpath(genpath('utils_functions'))

%% Dimensions to use in simulation
s=.0095; % noise level (NB actual Rician stdev depends on signal, see ricestat)
%Serena C. : s3=0.0057 about 3%noise
%            s5=0.0095 about 5% noise
%            s9=0.017 about 9% noise
%Vin:        max useful test: s=0.05

%% Set up test image
X_true = phantom;
hist_true=imhist(X_true);


%Other Noise
N = size(X_true,1); n = N^2;
g = X_true(:) + 0.05*max(X_true(:))*randn(n,1); %0.05  TODO: 0.09!!!!
X_noise=reshape(g,N,N);


X_noise = ricernd(X_true, s); %add rician noise





%% Compute and print

%mu=1;lambda=.09;Tol=10^(-4);beta=2;
%mu=1.9;lambda=.05;Tol=10^(-3);beta=.4;
%mu=1.9;lambda=0.06;Tol=10^(-3);
%mu=0.5;lambda=.04;
mu=1;lambda=0.000003; Tol=10^(-7);

Tol=10^(-3);

fprintf('Tol = %.4f \n', Tol)


fprintf('Split Isotropic Algorithm \n')
[u0,l0,k0] = SplitIsotropic2(X_noise,X_true,Tol,lambda,mu);
Str0='SBI-TV-old';
[outStr0,nEr0,psnr0]=printErPnsImages(Str0, u0, X_true);
fprintf(outStr0);
fprintf('Number of iterations %.0f \n', k0)



[u2,k2,rel2,l2] = SplitB_IsotropicTV(X_noise,lambda,mu,Tol, 50, X_true);
Str2='SBI-TV';
[outStr2,nEr2,psnr2]=printErPnsImages(Str2, u2, X_true);
fprintf(outStr2);


%% Match with Noise
kMax=50;
ImageToMatch=X_noise;
HistMatchType='ExactColtuc';
[u3,k3,rel3,l3, J3] = SB_IsotropicTV_2stage(X_noise,lambda,mu,ImageToMatch, HistMatchType,Tol, kMax, X_true);
Str3='SIA Exact Noise';
[outStr3,nEr3,psnr3]=printErPnsImages(Str3, u3, X_true);
fprintf(outStr3);

HistMatchType='UniformHistMatch';
[u4,k4,rel4,l4, J4] = SB_IsotropicTV_2stage(X_noise,lambda,mu,ImageToMatch, HistMatchType,Tol, kMax, X_true);
Str4='SIA Approx Noise';
[outStr4,nEr4,psnr4]=printErPnsImages(Str4, u4, X_true);
fprintf(outStr4);


%% Match with True
ImageToMatch=X_true;
HistMatchType='HistEq';
[u5,k5,rel5,l5, J5] = SB_IsotropicTV_2stage(X_noise,lambda,mu,ImageToMatch, HistMatchType,Tol, kMax, X_true);
Str5='SIA Approx True';
[outStr5,nEr5,psnr5]=printErPnsImages(Str5, u5, X_true);
fprintf(outStr5);

HistMatchType='ExactColtuc';
[u6,k6,rel6,l6, J6] = SB_IsotropicTV_2stage(X_noise,lambda,mu,ImageToMatch, HistMatchType,Tol, kMax, X_true);
Str6='SIA Exact True';
[outStr6,nEr6,psnr6]=printErPnsImages(Str6, u6, X_true);
fprintf(outStr6);




%mu=1.9;
Tol=10^(-3);
lambda=.005;
beta=.4;
%beta=0.005;
[u1,l1,k1] = SplitIsotropic2Lse(X_noise,X_true,Tol,lambda,mu,beta);
Str1='SI2Lse';
[outStr1,nEr1,psnr1]=printErPnsImages(Str1, u1, X_true);
fprintf(outStr1);




%% Display

figure('Name','Comparison among denoiser'); colormap gray;
zoom on;
subplot(3,3,1); imagesc(X_true); axis off image; title('Original'); zoom(2);
subplot(3,3,2); imagesc(X_noise); axis off image; title('Noisy'); zoom(2);
subplot(3,3,3); imagesc(u0); axis off image off; title(Str0); zoom(2);
subplot(3,3,4); imagesc(u1); axis image off; title(Str1); zoom(2);
subplot(3,3,5); imagesc(u2); axis image off; title(Str2); zoom(2);
subplot(3,3,6); imagesc(u3); axis image off ; title(Str3); zoom(2);
subplot(3,3,7); imagesc(u4); axis image off; title(Str4); zoom(2);
subplot(3,3,8); imagesc(u5); axis off image; title(Str5); zoom(2);
subplot(3,3,9); imagesc(u6); axis off image; title(Str6); zoom(2);
tightfig;


%Fig 2
figure('Name','Error vs. iteration')
i=(0:k0); semilogy(i, l0(i+1)); hold on;
i=(0:k1); semilogy(i, l1(i+1),'--')
i=(0:k2); semilogy(i, l2(i+1),'--')
xlabel('Number of Iterations')
ylabel('Error')
legend(Str0, Str1, Str2, 'Location', 'southeast')
legend('boxoff')

%Fig 2bis
figure('Name','Error vs. iteration')
i=(0:k0); semilogy(i, l0(i+1)); hold on;
i=(0:k1); semilogy(i, l1(i+1))
i=(0:k2); semilogy(i, l2(i+1))
i=(0:k3); semilogy(i, l3(i+1),'--')
i=(0:k4); semilogy(i, l4(i+1),'-.')
i=(0:k5); semilogy(i, l5(i+1),':')
i=(0:k6); semilogy(i, l6(i+1),'-.')
xlabel('Number of Iterations')
ylabel('Error')
legend(Str0, Str1, Str2, Str3, Str4, Str5, Str6, 'Location', 'southeast')
legend('boxoff')



%Fig 3
figure('Name','Jr vs. iteration')
i=(0:k3); semilogy(i, J3(i+1)); hold on;
i=(0:k4); semilogy(i, J4(i+1),'--')
i=(0:k5); semilogy(i, J5(i+1),':')
i=(0:k6); semilogy(i, J6(i+1),'-.')
xlabel('Number of Iterations')
ylabel('J')
legend(Str3, Str4, Str5, Str6, 'Location', 'southeast')
legend('boxoff')


%Fig 3
figure('Name','Rel. Error vs. iteration')
i=(0:k2); semilogy(i, rel2(i+1)); hold on;
i=(0:k3); semilogy(i, rel3(i+1))
i=(0:k4); semilogy(i, rel4(i+1),'--')
i=(0:k5); semilogy(i, rel5(i+1),':')
i=(0:k6); semilogy(i, rel6(i+1),'-.')
xlabel('Number of Iterations')
ylabel('Relative Error')
legend(Str2, Str3, Str4, Str5, Str6, 'Location', 'southeast')
legend('boxoff')


%Fig 4
figure('Name','Comparison of Rel. Err. and PDNR')
c = categorical({Str0, Str1, Str2, Str3, Str4, Str5, Str6});
subplot(211);
errors = [nEr0 nEr1 nEr2 nEr3 nEr4 nEr5 nEr6];
bar(c,errors,'BaseValue',min(errors)),set(gca,'yscale','log')
title('Err');
subplot(212);
Psnrs = [psnr0 psnr1 psnr2 psnr3 psnr4 psnr5 psnr6];
bar(c,Psnrs,'BaseValue',min(Psnrs))
title('PSNR Err');




%figure(3)
%plot(X(100,:))
%hold on;
%plot(u3(100,:))
%xlabel('Distance along profile')
%ylabel('Pixel intensity (Gray Value)')
%title('Intensity Profile (Isotropic TV-LSE)')
%legend('Original Image','Restored Image', 'Location', 'south')
%legend('boxoff')
