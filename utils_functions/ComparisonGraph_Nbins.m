function [figure_Comparison_Nbins]= ComparisonGraph_Nbins(end_nbins,X_noise, Theta, Tol, kMax, X_true,lambdaSP,mu,ImageToMatch)
for i=1:end_nbins
    nbins=2^i;
    [u_FP, kgraph, relErr, J, err] = ROFdenoiseNew(X_noise, Theta, Tol, kMax, X_true);
    [u_SB, k, relErr, J, err] = SB_IsotropicTV(X_noise,lambdaSP,mu,Tol, kMax, X_true);

    u_FP=imhistmatch(u_FP,ImageToMatch,nbins);
    u_SB=imhistmatch(u_SB,ImageToMatch,nbins);

    graph_relErr_FP{i} = norm(u_FP-X_true,2)/norm(X_true,2);
    
    graph_Psnr_FP{i}=psnr(u_FP,X_true);
    graph_relErr_SB{i} = norm(u_SB-X_true,2)/norm(X_true,2);
    graph_Psnr_SB{i}=psnr(u_SB,X_true);
end
x_graph=1:end_nbins;

figure_Comparison_Nbins=figure('Name','Psnr and RelErr vs. Nbins');

y_graph_SB=cell2mat(graph_Psnr_SB);
y_graph_FP=cell2mat(graph_Psnr_FP);
subplot(1,2,1);
hold on
plot(x_graph,y_graph_SB);
plot(x_graph,y_graph_FP);
ymax=max(y_graph_FP(6:end_nbins));
ymin=min(y_graph_SB(6:end_nbins));
axis([6 end_nbins ymin ymax]);
ylabel('Psnr');
xlabel('Nbins (power of 2)');
legend('SB Psnr','FP-ROF Psnr','Location','best')
hold off

y_graph_SB=cell2mat(graph_relErr_SB);
y_graph_FP=cell2mat(graph_relErr_FP);
subplot(1,2,2);
hold on
plot(x_graph,y_graph_SB);
plot(x_graph,y_graph_FP);
hold off
ymax=max(y_graph_SB(6:end_nbins));
ymin=min(y_graph_FP(6:end_nbins));
axis([6 end_nbins ymin ymax]);
ylabel('RelErr');
xlabel('Nbins (power of 2)');
legend('SB RelErr','FP-ROF RelErr','Location','best')
hold off

figure_Comparison_Nbins.Position = [200 200 1000 400];

