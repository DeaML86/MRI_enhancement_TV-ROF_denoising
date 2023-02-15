function [u, k, relErr, err, J]=SB_IsotropicTV_MultiStage(X_noise,lambda,mu, ImageToMatch, HistMatchType, Tol, kMax, X_true)
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

%% Initialization
tic
f=X_noise;  %f: original input image
[n,m]=size(f); %size of input image

u=f;
dx=zeros(n,m);
dy=zeros(n,m);
bx=zeros(n,m);
by=zeros(n,m);




%% Pre-first step
k=0;
err(1) = norm(f-X_true,2)/norm(X_true,2);
relErr(1) = norm(u,2);
[ux,uy]= gradient(u);% u at step 1
J(1)=mu/2*norm(u-f)^2 + sum(reshape(sqrt(ux.^2+uy.^2),[],1));



%% While -> SplitBregmanIsotropicIteration
while relErr(k+1) > Tol
    fprintf('it. %g ',k);
    k=k+1;
    uOld = u ;

    %compute u's %using dx, dx, i.e. grad d a
    i=2:(n-1);
    j=2:(m-1);
    u(i,j)=(lambda/(mu+4*lambda))*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)+dx(i,j-1)-dx(i,j)+dy(i-1,j)-dy(i,j)-bx(i,j-1)+bx(i,j)-by(i-1,j)+by(i,j))+(mu/(mu+4*lambda))*f(i,j);

    [ux,uy]= gradient(u);% u at step k+1
    J(k+1)=mu/2*norm(u-f)^2 + sum(reshape(sqrt(ux.^2+uy.^2),[],1));
    fprintf('J.=%g \t',J(k+1));
    
    %Compute the sx's  %s=zeros(size(u)) - preallocation non required.
    s = sqrt( abs(ux+bx).^2 + abs(uy+by).^2);
    
    %Compute the d's with u at step k!!!
    %dx = max(s-1/lambda,0).*((ux+bx)./s);
    %dy = max(s-1/lambda,0).*((uy+by)./s);
    dx = (s*lambda.*(ux+bx))./(s*lambda + 1);
    dy = (s*lambda.*(uy+by))./(s*lambda+1);

    %Compute the b's
    bx = bx + (ux - dx);
    by = by + (uy - dy);
    
    
    switch HistMatchType
        case 'ExactColtuc'
            u=im2double(exact_histogram(im2uint8(u),imhist(im2uint8(ImageToMatch))));
        case 'UniformHistMatch'
            u=imhistmatch(u,ImageToMatch); 
        case 'PolynomialHistMatch'
            %Works only for Matlab R2018?
            u=imhistmatch(u,imhist(ImageToMatch),'method','polynomial'); 
        case 'HistEq'
            u=histeq(u,imhist(ImageToMatch));        
        otherwise
            %No matching
            fprintf('******\n No Matching!! \n **********\n');
    end

    

    relErr(k+1) = norm(u-uOld,2)/norm(u,2);
    err(k+1)= norm(u-X_true,2)/norm(X_true,2);
    fprintf('err.=%g \t',err(k+1));
    fprintf('rel.err.=%g \n',relErr(k+1));
    
    if k>kMax
        fprintf('Stopped because n.iter > %d ; norm(up-u)/norm(u) = %.1e\n',kMax, relErr(k+1));
        return; 
    end
    
end
toc
fprintf('Stopped because norm(up-u)/norm(u) = %.1e\n', relErr(k+1));
end
