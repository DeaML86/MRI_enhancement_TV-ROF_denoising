function [u, k, relErr, err, J]=SB_ROF_VARIABLE_2stage(X_noise,lambdaSP,lambdaL1 ,mu,Theta, ImageToMatch, DenoiserType, HistMatchType, Tol, kMax, X_true)
%SplitB_IsotropicTV - Split Bregman Isotropic TV Denoising.
% This implement the 2stage method described 
%
% SplitB_IsotropicTV(f,lambda,mu,[Tol, kMax, X]) takes an image f,
% and applies the Split Bregman Isotropic TV Denoising Histogram matching
% with tolerance Tol, lambda and mu.
% If setting the ground truth X, also calculates the difference form the
% ground truth.
% The code is based on:
% - the article: THE SPLIT BREGMAN METHOD FOR L1 REGULARIZED PROBLEMS 
% of TOM GOLDSTEIN, STANLEY OSHER,
% - the Senior Thesis - Bregman Algorithms of Jacqueline Bush
% - the work on a similar code of Serena Crisci

% Vincenzo Schiano
% Feb 2018

% This code is released under the Gnu Public License (GPL). 
% For more information, see
% http://www.gnu.org/copyleft/gpl.html


%% stage1
k=0;
switch DenoiserType
    case 'StdSB'
        [u, k, relErr, J, err] = SB_IsotropicTV(X_noise,lambdaSP,mu,Tol, kMax, X_true);
    case 'ROF'
        [u, k, relErr, J, err] = ROFdenoiseNew(X_noise, Theta, Tol, kMax, X_true);
    case 'TVl1'
        [u, k, relErr, J, err] = TVL1denoiseNew(X_noise, lambdaL1, Tol, kMax, X_true);
    otherwise
        %No matching
        fprintf('******\n No Denoiser!! \n **********\n');
        u=X_noise;
        relErr=norm(u,2); 
        err=norm(u-X_true,2)/norm(X_true,2);
        J=0;
      
end



if (k(end)-1) > kMax
    fprintf('Stopped because n.iter > %d ; norm(up-u)/norm(u) = %.1e\n',kMax, nn(k+1));
end
    

%% stage2


%%Apply histogram Equalization 
%%Idea: can we apply it before the calculation of d's and b's?
uOld = u ;
switch HistMatchType
    case 'ExactColtuc'
        u=im2double(exact_histogram(im2uint8(u),imhist(im2uint8(ImageToMatch))));
    case 'UniformHistMatch' %Approx case
        u=imhistmatch(u,ImageToMatch,2^8); 
    case 'PolynomialHistMatch'
        u=imhistmatch(u,imhist(ImageToMatch),'method','polynomial'); 
    case 'HistEq'
        u=histeq(u,imhist(ImageToMatch));        
    otherwise
        %No matching
        fprintf('******\n No Matching!! \n **********\n');
end


if strcmp(HistMatchType,'ExactColtuc') || strcmp(HistMatchType,'UniformHistMatch') || strcmp(HistMatchType,'PolynomialHistMatch') || strcmp(HistMatchType,'HistEq')
    k=k+1;
    [ux,uy]= gradient(u);
    J(k+1)=mu/2*norm(u-X_true)^2 + sum(reshape(sqrt(ux.^2+uy.^2),[],1));
    fprintf('J.=%g \t',J(k+1));
    relErr(k+1) = norm(u-uOld,2)/norm(u,2);
    err(k+1)= norm(u-X_true,2)/norm(X_true,2);
    fprintf('err.=%g \t',err(k+1));
    fprintf('rel.err.=%g \n',relErr(k+1));
end



end
