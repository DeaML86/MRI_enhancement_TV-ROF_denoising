function [u,p,k]=SplitIsotropic2Lse(f,X,Tol,lambda,mu,beta)
tic
n=size(f,1); 
m=size(f,2);
n1=n-1; 
m1=m-1;
u=f; 
dx=zeros(n,m); 
dy=zeros(n,m); 
bx=zeros(n,m); 
by=zeros(n,m);
dxx=zeros(n,m);
dyy=zeros(n,m);
dxy=zeros(n,m);
bxx=zeros(n,m);
byy=zeros(n,m);
bxy=zeros(n,m);
k=0;
p(1) = norm(f-X,2)/norm(X,2);
nn(1) = norm(u,2);
while nn(k+1) >Tol
k=k+1;

[k nn(k) Tol];

u1 = u ;
for i=3:(n1-1)
for j=3:(m1-1)
  
u(i,j)=(lambda/(mu+22*lambda))*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-u(i+2,j)+(4*u(i+1,j))+(4*u(i-1,j))-u(i-2,j)-u(i,j+2)+(4*u(i,j+1))+(4*u(i,j-1))-u(i,j-2)-u(i+2,j+2)+(4*u(i+1,j+1))+(4*u(i-1,j-1))-u(i-2,j-2)+dx(i-1,j)-dx(i,j)+dy(i,j-1)-dy(i,j)-bx(i-1,j)+bx(i,j)-by(i,j-1)+by(i,j)-dxx(i-2,j)+(2*dxx(i-1,j))-dxx(i,j)-dyy(i,j-2)+(2*dyy(i,j-1))-dyy(i,j)-dxy(i-2,j-2)+(2*dxy(i-1,j-1))-dxy(i,j)+bxx(i-2,j)-(2*bxx(i-1,j))+bxx(i,j)+byy(i,j-2)-(2*byy(i,j-1))+byy(i,j)+bxy(i-2,j-2)-(2*bxy(i-1,j-1))+bxy(i,j))+(mu/(mu+22*lambda))*f(i,j);

end
end
s=zeros(size(u));
%Compute the sx’s
for i=2:n1
for j=2:m1
s(i,j) = sqrt( abs( (u(i+1,j)-u(i-1,j))/2 + bx(i,j))^2 + abs((u(i,j+1)-u(i,j-1))/2 + by(i,j))^2);
end
end
for i=1:n1
j=1;
s(i,j) = sqrt( abs( (u(i+1,j)-u(i,j)) + bx(i,j))^2 + abs((u(i,j+1)-u(i,j)) + by(i,j))^2);
end
for j=2:m1
i=1;
s(i,j) = sqrt( abs( (u(i+1,j)-u(i,j)) + bx(i,j))^2 + abs((u(i,j+1)-u(i,j)) + by(i,j))^2);
end
for i=2:n
j=m;
s(i,j) = sqrt( abs( (u(i,j)-u(i-1,j)) + bx(i,j))^2 + abs((u(i,j)-u(i,j-1)) + by(i,j))^2);
end
for j=2:m1
i=n;
s(i,j) = sqrt( abs( (u(i,j)-u(i-1,j)) + bx(i,j))^2 + abs((u(i,j)-u(i,j-1)) + by(i,j))^2);
end
for i=n
    j=1;
s(i,j) = sqrt( abs( (u(i,j)-u(i-1,j)) + bx(i,j))^2 + abs((u(i,j)-u(i,j+1)) + by(i,j))^2);    
end
for i=1
j=m;
  s(i,j) = sqrt( abs( (u(i,j)-u(i+1,j)) + bx(i,j))^2 + abs((u(i,j)-u(i,j-1)) + by(i,j))^2);  
end


%Compute the d’s
for i=2:n1
for j=1:m
dx(i,j)= (s(i,j)*lambda*((u(i+1,j)-u(i-1,j))/2+bx(i,j)))/(s(i,j)*lambda + 1);
end
end
for j=1:m
i=1;
dx(i,j)= (s(i,j)*lambda*((u(i+1,j)-u(i,j))+bx(i,j)))/(s(i,j)*lambda + 1);
end
for j=1:m
i=n;
dx(i,j)= (s(i,j)*lambda*((u(i,j)-u(i-1,j))+bx(i,j)))/(s(i,j)*lambda + 1);
end
for i=1:n
for j=2:m1
dy(i,j) = (s(i,j)*lambda*(((u(i,j+1)-u(i,j-1))/2 + by(i,j))))/(s(i,j)*lambda +1 );
end
end
for i=1:n
j=1;
dy(i,j) = (s(i,j)*lambda*(((u(i,j+1)-u(i,j)) + by(i,j))))/(s(i,j)*lambda +1 );
end
for i=1:n
j=m;
dy(i,j) = (s(i,j)*lambda*(((u(i,j)-u(i,j-1))/2 + by(i,j))))/(s(i,j)*lambda +1 );
end
%Compute the b’s
for i=2:n1
for j=1:m
bx(i,j)= bx(i,j) + ((u(i+1,j)-u(i-1,j))/2 - dx(i,j));
end
end
for j=1:m
i=1;
bx(i,j)= bx(i,j) + ((u(i+1,j)-u(i,j)) - dx(i,j));
end
for j=1:m
i=n;
bx(i,j)= bx(i,j) + ((u(i,j)-u(i-1,j)) - dx(i,j));
end
for i=1:n
for j=2:m1
by(i,j) = by(i,j) + ((u(i,j+1)-u(i,j-1))/2 - dy(i,j));
end
end
for i=1:n
j=1;
by(i,j) = by(i,j) + ((u(i,j+1)-u(i,j)) - dy(i,j));
end
for i=1:n
j=m;
by(i,j) = by(i,j) + ((u(i,j)-u(i,j-1)) - dy(i,j));
end
t=zeros(size(u));
%aggiunto da Serena 
%Compute the tx's
for i=2:n1
for j=2:m1
    t(i,j)=2+u(i+1,j)-2*u(i,j)+u(i-1,j)+u(i,j+1)-2*u(i,j)+u(i,j-1)+bxx(i,j)+byy(i,j); 
end
end
for i=1:(n1-1)
j=1;
t(i,j)=2+u(i+2,j)-2*u(i+1,j)+u(i,j)+u(i,j+2)-2*u(i,j+1)+u(i,j)+bxx(i,j)+byy(i,j);
end
for j=2:(m1-1)
i=1;
t(i,j)=2+u(i+2,j)-2*u(i+1,j)+u(i,j)+u(i,j+2)-2*u(i,j+1)+u(i,j)+bxx(i,j)+byy(i,j);
end
for i=3:n
j=m;
t(i,j)=2+u(i-2,j)-2*u(i-1,j)+u(i,j)+u(i,j-2)-2*u(i,j-1)+u(i,j)+bxx(i,j)+byy(i,j);
end
for j=3:m1
i=n;
t(i,j)=2+u(i-2,j)-2*u(i-1,j)+u(i,j)+u(i,j-2)-2*u(i,j-1)+u(i,j)+bxx(i,j)+byy(i,j);
end
for i=n
j=2;
t(i,j)=2+u(i-2,j)-2*u(i-1,j)+u(i,j)+u(i,j+2)-2*u(i,j+1)+u(i,j)+bxx(i,j)+byy(i,j);  
end
for i=n1:n
j=1;
t(i,j)=2+u(i-2,j)-2*u(i-1,j)+u(i,j)+u(i,j+2)-2*u(i,j+1)+u(i,j)+bxx(i,j)+byy(i,j);  
end
for j=m
i=2;
t(i,j)=2+u(i+2,j)-2*u(i+1,j)+u(i,j)+u(i,j-2)-2*u(i,j-1)+u(i,j)+bxx(i,j)+byy(i,j);  
end
for j=m1:m
i=1;
t(i,j)=2+u(i+2,j)-2*u(i+1,j)+u(i,j)+u(i,j-2)-2*u(i,j-1)+u(i,j)+bxx(i,j)+byy(i,j);  
end

%Compute the dxx's
for i=2:n1
for j=1:m
    dxx(i,j)=((lambda*t(i,j))*(u(i+1,j)-2*u(i,j)+u(i-1,j) +bxx(i,j))-beta)/(lambda*t(i,j)+beta);
end
end
for j=1:m
i=1;
dxx(i,j)=((lambda*t(i,j))*(u(i+2,j)-2*u(i+1,j)+u(i,j) +bxx(i,j))-beta)/(lambda*t(i,j)+beta);
end
for j=1:m
i=n;
dxx(i,j)=((lambda*t(i,j))*(u(i-2,j)-2*u(i-1,j)+u(i,j) +bxx(i,j))-beta)/(lambda*t(i,j)+beta);
end
%Compute the dyy's
for i=1:n
for j=2:m1
    dyy(i,j)=((lambda*t(i,j))*(u(i,j+1)-2*u(i,j)+u(i,j-1) +byy(i,j))-beta)/(lambda*t(i,j)+beta);
end
end
for i=1:n
j=1;
dyy(i,j)=((lambda*t(i,j))*(u(i,j+2)-2*u(i,j+1)+u(i,j) +byy(i,j))-beta)/(lambda*t(i,j)+beta);
end
for i=1:n
j=m;
dyy(i,j)=((lambda*t(i,j))*(u(i,j-2)-2*u(i,j-1)+u(i,j) +byy(i,j))-beta)/(lambda*t(i,j)+beta);
end
%Compute the dxy's
for i=2:n1
for j=2:m1
    dxy(i,j)= u(i+1,j+1) -2*u(i,j)+u(i-1,j-1)+bxy(i,j)- (beta/lambda);
end
end
for i=1:(n1-1)
j=1;
 dxy(i,j)= u(i+2,j+2) -2*u(i+1,j+1)+u(i,j)+bxy(i,j)- (beta/lambda);
end
for j=2:(m1-1)
i=1;
dxy(i,j)= u(i+2,j+2) -2*u(i+1,j+1)+u(i,j)+bxy(i,j)- (beta/lambda);
end
for i=3:n
j=m;
dxy(i,j)= u(i-2,j-2) -2*u(i-1,j-1)+u(i,j)+bxy(i,j)- (beta/lambda);
end
for j=3:m1
i=n;
dxy(i,j)= u(i-2,j-2) -2*u(i-1,j-1)+u(i,j)+bxy(i,j)-(beta/lambda);
end
for i=n
j=2;
dxy(i,j)= u(i-2,j+2) -2*u(i-1,j+1)+u(i,j)+bxy(i,j)-(beta/lambda);
end
for i=n1:n
j=1;
dxy(i,j)= u(i-2,j+2) -2*u(i-1,j+1)+u(i,j)+bxy(i,j)-(beta/lambda);
end
for j=m
i=2;
dxy(i,j)= u(i+2,j-2) -2*u(i+1,j-1)+u(i,j)+bxy(i,j)-(beta/lambda);    
end
for j=m1:m
i=1;
dxy(i,j)= u(i+2,j-2) -2*u(i+1,j-1)+u(i,j)+bxy(i,j)-(beta/lambda);
end

%Compute the bxx's
for i=2:n1
for j=1:m
bxx(i,j)= bxx(i,j) + ((u(i+1,j)-2*u(i,j)+u(i-1,j)) - dxx(i,j));
end
end
for j=1:m
i=1;
bxx(i,j)= bxx(i,j) + ((u(i+2,j)-2*u(i+1,j)+u(i,j)) - dxx(i,j));
end
for j=1:m
i=n;
bxx(i,j)= bxx(i,j) + ((u(i-2,j)-2*u(i-1,j)+u(i,j)) - dxx(i,j));
end
%Compute the byy's
for i=1:n
for j=2:m1
byy(i,j) = byy(i,j) + ((u(i,j+1)-2*u(i,j)+ u(i,j-1)) - dyy(i,j));
end
end
for i=1:n
j=1;
byy(i,j) = byy(i,j) + ((u(i,j+2)-2*u(i,j+1)+u(i,j)) - dyy(i,j));
end
for i=1:n
j=m;
byy(i,j) = byy(i,j) + ((u(i,j)-2*u(i,j-1)+u(i,j-2)) - dyy(i,j));
end
%Compute the bxy's
for i=2:n1
for j=2:m1
    bxy(i,j)=bxy(i,j)+((u(i+1,j+1)-2*u(i,j)+u(i-1,j-1))- dxy(i,j)); 
end
end
for i=1:(n1-1)
j=1;
 bxy(i,j)=  bxy(i,j)+ ((u(i+2,j+2)-2*u(i+1,j+1)+u(i,j))-dxy(i,j));
end
for j=2:(m1-1)
i=1;
bxy(i,j)=  bxy(i,j)+ ((u(i+2,j+2)-2*u(i+1,j+1)+u(i,j))-dxy(i,j));
end
for i=3:n
j=m;
bxy(i,j)= bxy(i,j)+((u(i-2,j-2)-2*u(i-1,j-1)+u(i,j))-dxy(i,j));
end
for j=3:m
i=n;
bxy(i,j)=  bxy(i,j)+((u(i-2,j-2)-2*u(i-1,j-1)+u(i,j))-dxy(i,j));
end
for i=n
j=2;
bxy(i,j)=  bxy(i,j)+((u(i-2,j+2)-2*u(i-1,j+1)+u(i,j))-dxy(i,j));
end
for i=n1:n
j=2;
bxy(i,j)=  bxy(i,j)+((u(i-2,j+2)-2*u(i-1,j+1)+u(i,j))-dxy(i,j));
end
for j=m
i=2;
bxy(i,j)=  bxy(i,j)+((u(i+2,j-2)-2*u(i+1,j-1)+u(i,j))-dxy(i,j));
end
for j=m1:m
i=2;
bxy(i,j)=  bxy(i,j)+((u(i+2,j-2)-2*u(i+1,j-1)+u(i,j))-dxy(i,j));
end

nn(k+1) = norm(u-u1,2)/norm(u,2);
p(k+1)= norm(u-X,2)/norm(X,2);
end
toc;
%psnr_iso= 10*log10(max(max(X))^2/mse(X,u))

end