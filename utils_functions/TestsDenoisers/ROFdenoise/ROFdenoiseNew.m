%% ROFdenoise
%
%  This denoising method is based on total-variation, originally proposed by
%  Rudin, Osher and Fatemi. In this particular case fixed point iteration
%  is utilized.
%
%  For the included image, a fairly good result is obtained by using a
%  theta value around 12-16. A possible addition would be to analyze the
%  residual with an entropy function and add back areas that have a lower
%  entropy, i.e. there are some correlation between the surrounding pixels.
% 
%  Philippe Magiera & Carl L�ndahl, 2008
%

function [u, k, nn, J, err] = ROFdenoiseNew(Image, Theta, tol, nbrOfMaxIterations, X_true,norma)
  
if(nargin<3) || isempty(tol)
	tol=10^(-5);
end
g = 1; dt = 1/4; 

if(nargin<3) || isempty(tol)
	nbrOfMaxIterations = 5;
end

if(nargin<6) || isempty(norma)
	norma=1;
end

[Image_h, Image_w] = size(Image); 
Image = im2double(Image);
p = zeros(Image_h,Image_w,2);
d = zeros(Image_h,Image_w,2);


%% Zero step
tic

u=Image;
k=0;
nn(1)= norm(Image,norma);
%norma
err(1) = norm(Image-X_true,norma)/norm(X_true,norma);

% Calculate forward derivatives
du(:,:,2) = u(:,[2:Image_w, Image_w])-u;
du(:,:,1) = u([2:Image_h, Image_h],:)-u;

% Iterate
d(:,:,1) = (1+(dt/Theta/g).*abs(sqrt(du(:,:,1).^2+du(:,:,2).^2)));
d(:,:,2) = (1+(dt/Theta/g).*abs(sqrt(du(:,:,1).^2+du(:,:,2).^2)));
p = (p-(dt/Theta).*du)./d;

while nn(k+1) > tol
    k=k+1;
    fprintf('it. %g ',k);
    uOld = u ;
    
    p_slice_1=p((1:Image_h-1), :, 1);
    %p_slice_1=p((2:Image_h-1), :, 1);zeros(1, Image_w)
    div_p= [ p_slice_1; zeros(1, Image_w)]  -  [zeros(1, Image_w);p_slice_1]; 
    
    p_slice_2=p(:, (1:Image_w-1), 2);
    %p_slice_2=p(:, (2:Image_w-1), 2);
    div_p=[p_slice_2 zeros(Image_h, 1)] - [zeros(Image_h, 1) p_slice_2] + div_p;
    
%     % Handle boundaries
%     div_p(:,1) = p(:,1,2);
%     div_p(:,Image_w) = -p(:,Image_w-1,2);
%     div_p(1,:) = p(1,:,1);
%     div_p(Image_h,:) = -p(Image_h-1,:,1);
    
    % Update u
    u = Image-Theta*div_p;

    % Calculate forward derivatives
    du(:,:,2) = u(:,[2:Image_w, Image_w])-u;
    du(:,:,1) = u([2:Image_h, Image_h],:)-u;

    % Iterate
    d(:,:,1) = (1+(dt/Theta/g).*abs(sqrt(du(:,:,1).^2+du(:,:,2).^2)));
    d(:,:,2) = (1+(dt/Theta/g).*abs(sqrt(du(:,:,1).^2+du(:,:,2).^2)));
    p = (p-(dt/Theta).*du)./d;
    
    
    nn(k+1) = norm(u-uOld,norma)/norm(u,norma);
    err(k+1) = norm(u-X_true,norma)/norm(X_true,norma);
    J(k+1)=0;
    fprintf('rel.err.=%g \n',nn(k));
    
    if k>=nbrOfMaxIterations
        fprintf('Stopped because n.iter >= %d ; norm(up-u)/norm(u) = %.1e\n',nbrOfMaxIterations, nn(k));
        break ; 
    end
      
end

toc
fprintf('Stopped because norm(up-u)/norm(u) = %.1e\n', nn(k));

%psnr_iso= 10*log10(max(max(X))^2/mse(X,u));

end
