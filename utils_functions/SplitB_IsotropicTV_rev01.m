function [u, k, nn, varargout]=SplitB_IsotropicTV_Definition(f,lambda,mu, varargin)
%SplitB_IsotropicTV - Split Bregman Isotropic TV Denoising.
% This tries to implement the exaxt definiton of the methond described by
% T.Goldstein and S.Osher https://doi.org/10.1137/080725891
%
% SplitB_IsotropicTV(f,lambda,mu,[Tol, kMax, X]) takes an image f,
% and applies the Split Bregman Isotropic TV Denoising with
% tolerance Tol, lambda and mu.
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
[n,m]=size(f); %size of input image

u=f;
dx=zeros(n,m);
dy=zeros(n,m);
bx=zeros(n,m);
by=zeros(n,m);


kMax=20;
Tol=10^(-3);

if size(varargin,2) >= 1; Tol=varargin{1}; end
if size(varargin,2) >= 2; kMax=varargin{2}; end
if size(varargin,2) >= 3; X=varargin{3}; end %X true!

%% Pre-first step
k=0;
if size(varargin,2)==3; p(1) = norm(f-X,2)/norm(X,2); varargout{1}=p; end
nn(1) = norm(u,2);
[ux,uy]= gradient(u); %returns the x and y components of the two-dimensional numerical gradient of matrix u.

%% While -> SplitBregmanIsotropicIteration
while nn(k+1) > Tol
    fprintf('it. %g ',k);
    k=k+1;
    uOld = u ;

    %Compute the sx's  %s=zeros(size(u)) - preallocation non required.
    s = sqrt( abs(ux+bx).^2 + abs(uy+by).^2);
    
    %Compute the d's with u at step k!!!
    %dx_new = max(s-1/lambda,0).*((ux+bx)./s);
    %dy_new = max(s-1/lambda,0).*((uy+by)./s);
    dx_new = (s*lambda.*(ux+bx))./(s*lambda + 1);
    dy_new = (s*lambda.*(uy+by))./(s*lambda+1);
    
    
    %compute u's %using dx, dx, i.e. grad d a
    i=2:(n-1);
    j=2:(m-1);
    u(i,j)=(lambda/(mu+4*lambda))*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)+dx(i,j-1)-dx(i,j)+dy(i-1,j)-dy(i,j)-bx(i,j-1)+bx(i,j)-by(i-1,j)+by(i,j))+(mu/(mu+4*lambda))*f(i,j);
    
    
    [ux,uy]= gradient(u);% u at step k+1
    
    dx=dx_new;
    dy=dy_new;
    
    %Compute the b's
    bx = bx + (ux - dx);
    by = by + (uy - dy);
    
    nn(k+1) = norm(u-uOld,2)/norm(u,2);
    if size(varargin,2)==3
        p(k+1)= norm(u-X,2)/norm(X,2);
        fprintf('err.=%g \t',p(k+1));
        varargout{1}=p;  
    end
    
    fprintf('rel.err.=%g \n',nn(k+1));
    
    if k>kMax
        fprintf('Stopped because n.iter > %d ; norm(up-u)/norm(u) = %.1e\n',kMax, nn(k+1));
        return; 
    end
    
end
toc
fprintf('Stopped because norm(up-u)/norm(u) = %.1e\n', nn(k+1));

%psnr_iso= 10*log10(max(max(X))^2/mse(X,u));

end
