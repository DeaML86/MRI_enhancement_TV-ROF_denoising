function [u,p,k]=SplitIsotropic2(f,X,Tol,lambda,mu)
tic
n=size(f,1); 
m=size(f,2);
n1=n-1; 
m1=m-1;
u=f; 
dx =zeros(n,m); 
dy=zeros(n,m); 
bx=zeros(n,m); 
by=zeros(n,m);
k=0;
p(1) = norm(f-X,2)/norm(X,2);
nn(1) = norm(u,2);
while nn(k+1) >Tol
k=k+1;

%[k nn(k) Tol]

u1 = u ;
for i=2:n1
for j=2:m1
    
u(i,j)=(lambda/(mu+4*lambda))*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)+dx(i-1,j)-dx(i,j)+dy(i,j-1)-dy(i,j)-bx(i-1,j)+bx(i,j)-by(i,j-1)+by(i,j))+(mu/(mu+4*lambda))*f(i,j);%+ max((lambda/(mu+4*lambda))*(u(i+1,j)-2*u(i,j)+u(i,j-1)),(lambda/(mu+4*lambda))*(u(i,j+1)-2*u(i,j)+u(i,j-1)));

end
end
s=zeros(size(u));
%Compute the sx�s
for i=2:n1
for j=2:m1
s(i,j) = sqrt( abs( (u(i+1,j)-u(i-1,j))/2 + bx(i,j))^2 + abs((u(i,j+1)-u(i,j-1))/2 + by(i,j))^2);
end
end
for i=1:n1
j=1;
s(i,j) = sqrt( abs( (u(i+1,j)-u(i,j)) + bx(i,j))^2 + abs((u(i,j+1)-u(i,j)) + by(i,j))^2);
end
for j=1:m1
i=1;
s(i,j) = sqrt( abs( (u(i+1,j)-u(i,j)) + bx(i,j))^2 + abs((u(i,j+1)-u(i,j)) + by(i,j))^2);
end
for i=2:n
j=m;
s(i,j) = sqrt( abs( (u(i,j)-u(i-1,j)) + bx(i,j))^2 + abs((u(i,j)-u(i,j-1)) + by(i,j))^2);
end
for j=2:m
i=n;
s(i,j) = sqrt( abs( (u(i,j)-u(i-1,j)) + bx(i,j))^2 + abs((u(i,j)-u(i,j-1)) + by(i,j))^2);
end
for j=1:m1
i=1;
s(i,j) = sqrt( abs( (u(i+1,j)-u(i,j)) + bx(i,j))^2 + abs((u(i,j+1)-u(i,j)) + by(i,j))^2);
end
%Compute the d�s
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
%Compute the b�s
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
nn(k+1) = norm(u-u1,2)/norm(u,2);
p(k+1)= norm(u-X,2)/norm(X,2);
    if k>300; return; end
end
toc
%psnr_iso= 10*log10(max(max(X))^2/mse(X,u));

end