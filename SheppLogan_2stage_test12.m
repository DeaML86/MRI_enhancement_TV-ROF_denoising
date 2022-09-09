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
mu=1;
kMax=25;
MatchType={'None', 'UniformHistMatch',  'UniformHistMatch'};
Str= {'SplitB TV Std', 'SbTV Approx Noise', 'SbTV Approx True'};
nOfTests=size(MatchType,2);


for s_nois=0.0045:0.0045:0.02

X_noise = ricernd(X_true, s); %add rician noise

for lamb_i=1:6

lambda=3*10^(-1-lamb_i);
Tol=10^(-2-lamb_i);
fprintf('Tol = %.4f \n', Tol)

ImageToMatch=X_noise;
for i=1:nOfTests-2
    HistMatchType=string(MatchType(i));
    [u{i},k{i},rel{i},l{i}, J{i}] = SB_IsotropicTV_2stage(X_noise,lambda,mu,ImageToMatch, HistMatchType,Tol, kMax, X_true);
    [outStr{i},nEr(i),psnr(i)]=printErPnsImages(string(Str{i}), u{i}, X_true);
    fprintf(string(outStr{i}));
end

ImageToMatch=X_true;

for i = (nOfTests-2) : nOfTests
    fprintf(string(Str{i}));
    HistMatchType = string(MatchType{i});
    [u{i},k{i},rel{i},l{i}, J{i}] = SB_IsotropicTV_2stage(X_noise,lambda,mu,ImageToMatch, HistMatchType,Tol, kMax, X_true);
    [outStr{i},nEr(i),psnr(i)]=printErPnsImages(string(Str{i}), u{i}, X_true);
    fprintf(string(outStr{i}));
end



%% Display

figure('Name', ['Comparison among denoiser' 'Noise'  num2str(s) 'Lambda' num2str(lambda,'%.1e') ] ); colormap gray;
subplot(3,3,1); imagesc(X_true); axis off image; title('Original'); zoom(2);
subplot(3,3,2); imagesc(X_noise); axis off image; title('Noisy'); zoom(2);
for i=1:nOfTests
    subplot(3,3,2+i); imagesc(u{i}); axis image off; title(string(Str{i})); zoom(2);
end
tightfig;


%Fig 2
figure('Name',['Error vs. iteration' 'Noise'  num2str(s) 'Lambda' num2str(lambda,'%.1e')])
i=(0:k{1}); semilogy(i, l{1}(i+1)); hold on;
for j=2:nOfTests
    i=(0:k{j}); semilogy(i, l{j}(i+1),'--')
end
xlabel('Number of Iterations')
ylabel('Error')
legend(Str, 'Location', 'southeast')
legend('boxoff')



%Fig 3
figure('Name',['Jr vs. iteration' 'Noise'  num2str(s) 'Lambda' num2str(lambda,'%.1e')])
i=(0:k{1}); semilogy(i, J{1}(i+1)); hold on;
for j=2:nOfTests
    i=(0:k{j}); semilogy(i, J{j}(i+1),'-.')
end
xlabel('Number of Iterations')
ylabel('J')
legend(Str, 'Location', 'southeast')
legend('boxoff')


%Fig 3
figure('Name',['Rel. Error vs. iteration' 'Noise'  num2str(s) 'Lambda' num2str(lambda,'%.1e')])
i=(0:k{1}); semilogy(i, rel{1}(i+1)); hold on;
for j=2:nOfTests
    i=(0:k{j}); semilogy(i, rel{j}(i+1),'-.')
end
xlabel('Number of Iterations')
ylabel('Relative Error')
legend(Str, 'Location', 'southeast')
legend('boxoff')


%Fig 4
figure('Name',['Comparison of Rel. Err. and PDNR' 'Noise'  num2str(s) 'Lambda' num2str(lambda,'%.1e')])
c = categorical(Str);
subplot(2,1,1);
errors = nEr;
bar(c,errors,'BaseValue',min(errors)),set(gca,'yscale','log')
title('Err');
subplot(2,1,2);
Psnrs = psnr;
bar(c,Psnrs,'BaseValue',min(Psnrs))
title('PSNR Err');

end
end

%figure(3)
%plot(X(100,:))
%hold on;
%plot(u3(100,:))
%xlabel('Distance along profile')
%ylabel('Pixel intensity (Gray Value)')
%title('Intensity Profile (Isotropic TV-LSE)')
%legend('Original Image','Restored Image', 'Location', 'south')
%legend('boxoff')
