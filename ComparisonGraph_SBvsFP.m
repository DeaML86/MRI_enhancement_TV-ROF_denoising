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

slice=75;
outImage= strcat('outImages/SBvsFP_graph') ;
SetParamScript;

lambdaSP=4.95*10^(-2);
lambdaL1= 2.5;
Theta= (2/7)*10^(-1);

Kstart=2; KBins=5; kMax=30;
norma="inf";
lambda=lambdaSP;
mu=1; 
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
    [y_Psnr_SB{i},y_Psnr_FP{i},y_RelErr_SB{i},y_RelErr_FP{i}] = ComparisonGraph_SBvsFP_Yreturn( Theta,Kstart, KBins, kMax,lambdaSP,mu,ref,s,shift,slice,norma);
end

x_graph=Kstart:KBins:kMax;
figure_Comparison_SBvsFP=figure('Name','Psnr and RelErr vs. Nbins');

subplot(1,2,1);
hold on
for i=1:4
    y_graph_SB=cell2mat( y_Psnr_SB(i) );
    y_graph_FP=cell2mat( y_Psnr_FP(i) );
    plot(x_graph,y_graph_SB);
    plot(x_graph,y_graph_FP);
end
ylabel('Psnr (higher)');
xlabel('Iteration');
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
ylabel('RelErr (Lower) ');
xlabel('Iteration');
legend('SB case1','FP case1','SB case2','FP case2','SB case3','FP case3','SB case4','FP case4','Location','northeast')
hold off

figure_Comparison_SBvsFP.Position = [200 200 1100 400];
exportgraphics(figure_Comparison_SBvsFP, fullfile(outImage, strcat('Graph_Comparison_SBvsFP_',num2str(Kstart),'-',num2str(kMax),'_norm',string(norma),'.jpg')),"resolution",300)
exportgraphics(figure_Comparison_SBvsFP, fullfile(outImage, strcat('Graph_Comparison_SBvsFP_',num2str(Kstart),'-',num2str(kMax),'_norm',string(norma),'.png')),"resolution",300)