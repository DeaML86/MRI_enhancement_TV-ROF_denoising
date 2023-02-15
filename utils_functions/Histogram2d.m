function [f]=Histogram2d(D1,D2,title1,title2)

X = D1;
Y = D2;
data = [X,Y];
nfont=15;

f=figure();
hold on
histogram2(X, Y,2^8,DisplayStyle="tile") ;
axis([0 1 0 1]);
ylabel(title1,'FontSize',nfont+3);
xlabel(title2,'FontSize',nfont+3);


% cx= colorbar();

% yellow	[1 1 0]	[255 255 0]
% magenta	[1 0 1]	[255 0 255]
% cyan	[0 1 1]	[0 255 255]
% red	[1 0 0]	[255 0 0]
% green	[0 1 0]	[0 255 0]
% blue	[0 0 1]	[0 0 255]
% white	[1 1 1]	[255 255 255]
% black	[0 0 0]

n2 = 50;               %// number of colors
n1=60;

R1 = linspace(0,0.4,n1);  %// Red from 1 to 0
R2 = linspace(0.4,1,n2);  %// Red from 0 to 1
G1 = linspace(0,1,n1);  %// Green from 0 to 1
G2 = linspace(1,0.9,n2);  %// Green from 0 to 1
B1 = linspace(0.85,0.1,n1);  %// Blue from 0 to 1
B2 = linspace(0.1,0,n2);  %// Blue from 0 to 1

R=[R1,R2];
B=[B1,B2];
G=[G1,G2];
map=[R(:), G(:), B(:)];

map='parula';
colormap( map );

clim([1 20])
% cx.Ruler.Scale = 'log';
% cx.Ruler.MinorTick = 'on';

ax = gca; 
ax.YAxis.FontSize = nfont;
ax.XAxis.FontSize =nfont;
t=linspace(0,1);
plot(t,t,'red','LineWidth',1)
plot(t,t+0.1,'--r','LineWidth',1)
plot(t,t-0.1,'--r','LineWidth',1)
hold off