% Graph paper

clc; % Clear command window
close all; % Close figures
% Clear loaded functions: % (avoid "clear all")
clear all;
clear('functions');
clear('classes');
clear('java');
clear('global');
%clear('import');  % Not from inside a function!
clear('variables');
SetParamScript;
slice=75;
end_nbins=22;
outImage= strcat('outImages/Nbins_graph') ;

lambdaSP=4.95*10^(-2);
lambdaL1= 2.5;
Theta= (2/7)*10^(-1);
Tol=10^(-3);
lambda=lambdaSP;
mu=1; kMax=700;

for i=1:4
    switch i
        case 1
            ref=0; s=9; shift=0; 
        case 2
            ref=0; s=9; shift=5; 
        case 3
            ref=3; s=9; shift=5;
        case 4
            ref=3; s=9; shift=-5; 
    end
    [y_Psnr_SB{i},y_Psnr_FP{i},y_RelErr_SB{i},y_RelErr_FP{i}]= ComparisonGraph_Nbins_Yreturn(end_nbins,Theta, Tol, kMax,lambdaSP,mu,ref,s,shift,slice);
end

x_graph=1:end_nbins;
Choice1=8;
Choice2=16;
figure_Comparison_Nbins=figure('Name','Psnr and RelErr vs. Nbins');

subplot(1,2,1);
hold on
for i=1:4
    y_graph_SB=cell2mat( y_Psnr_SB(i) );
    y_graph_FP=cell2mat( y_Psnr_FP(i) );
    plot(x_graph,y_graph_SB);
    plot(x_graph,y_graph_FP);
end
xline(Choice1,'--r');
xline(Choice2,'--r');
ylabel('Psnr (higher)');
xlabel('Nbins (exponent of 2)');
xlim([6 21]); xticks([6:21]);
legend('SB case1','FP case1','SB case2','FP case2','SB case3','FP case3','SB case4','FP case4','Location','northeast')
hold off,

subplot(1,2,2);
hold on
for i=1:4
    y_graph_SB=cell2mat( y_RelErr_SB(i) );
    y_graph_FP=cell2mat( y_RelErr_FP(i) );
    plot(x_graph,y_graph_SB);
    plot(x_graph,y_graph_FP);
end
xline(Choice1,'--r');
xline(Choice2,'--r');
ylabel('RelErr (Lower) ');
xlabel('Nbins (exponent of 2)');
xlim([6 21]); xticks([6:21]);
legend('SB case1','FP case1','SB case2','FP case2','SB case3','FP case3','SB case4','FP case4','Location','northeast')
hold off

figure_Comparison_Nbins.Position = [200 200 1100 400];
exportgraphics(figure_Comparison_Nbins, fullfile(outImage, strcat('Graph_Comparison_Nbins','.jpg')),"resolution",300)
exportgraphics(figure_Comparison_Nbins, fullfile(outImage, strcat('Graph_Comparison_Nbins','.png')),"resolution",300)