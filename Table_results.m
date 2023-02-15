% Table paper

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
addpath(genpath('TestsDenoisers'))

caso=1;
switch caso
        case 1
            ref=0; s=9; shift=0; 
        case 2
            ref=0; s=9; shift=5; 
        case 3
            ref=3; s=9; shift=5;
        case 4
            ref=3; s=9; shift=-5;
end

outImage= strcat('outImages/Ref_',num2str(ref),'-Noise_', num2str(s),'-Shift_',num2str(shift)) ;

slice=75;

% Set up test image
pulseSequence='t1';
file=fullfile(imagePath, pulseSequence, 't1_icbm_normal_1mm_pn0_rf20.mnc');
mriData = loadminc(file);
fprintf('Dimension MRI %dx%dx%d\n', size(mriData))
X_true = mriData(:,:,slice); % (particular slice through the brain)
%hist_true=imhist(X_true);


% Noise section
fileNoise=fullfile(imagePath, pulseSequence,strcat('t1_icbm_normal_1mm_pn', num2str(s),'_rf20.mnc') );
mriDataNoise= loadminc(fileNoise);
fprintf('Dimensione MRI %dx%dx%d\n', size(mriDataNoise))
X_noise=mriDataNoise(:,:,slice);


% ImageToMatch section
fileImgToMatch=fullfile(imagePath, pulseSequence,strcat('t1_icbm_normal_1mm_pn', num2str(ref),'_rf20.mnc') );
mriDataImgToMatch= loadminc(fileImgToMatch);
fprintf('Dimensione MRI %dx%dx%d\n', size(mriDataImgToMatch))
ImageToMatch=mriDataImgToMatch(:,:,slice-shift);



%Fix Image with data from 0 to 1.

minMRI = min(min(min(mriDataNoise(:),mriData(:))));
maxMRI = max(max(mriDataNoise(:),mriData(:)));

X_noise= double((X_noise - minMRI)/(maxMRI - minMRI)); %quantizzOnlyNormaliz  quantizzazione
X_true = double((X_true - minMRI)/(maxMRI - minMRI)); % this helps to get double from 0 to 1.
ImageToMatch= double((ImageToMatch - minMRI)/(maxMRI - minMRI));

histo_True=imhist(X_true,2^8);
Histo_True_figure=figure('Name','Histogram X True');
bar(histo_True);
title('Histogram X True','FontSize',20);
exportgraphics(Histo_True_figure, fullfile(outImage, strcat('Histogram X True','.jpg')),"resolution",300)
exportgraphics(Histo_True_figure, fullfile(outImage, strcat('Histogram X True','.png')),"resolution",300)

histo_Noise=imhist(X_noise,2^8);
Histo_Noise_figure=figure('Name','Histogram X Noise');
bar(histo_Noise)
title(strcat('Histogram X Noise ', num2str(s)),'FontSize',20);

exportgraphics(Histo_Noise_figure, fullfile(outImage, strcat('Histogram X Noise','.jpg')),"resolution",300)
exportgraphics(Histo_Noise_figure, fullfile(outImage, strcat('Histogram X Noise','.png')),"resolution",300)

Histo=Histogram2d(X_true,X_noise,  'X True','X Noise');
exportgraphics(Histo, fullfile(outImage, strcat('H2d_true_noise', num2str(s),'.jpg')),"resolution",300)
exportgraphics(Histo, fullfile(outImage, strcat('H2d_true_noise', num2str(s),'.png')),"resolution",300)

%%

mu=1; 
kMax=700;
DenoiserType={'None','None','StdSB', 'StdSB', 'StdSB','ROF', 'ROF', 'ROF','None','None'};
MatchType={'ExactColtuc','UniformHistMatch','None', 'ExactColtuc', 'UniformHistMatch', 'None', 'ExactColtuc', 'UniformHistMatch','None','None'};
Str= {'ExactMatch','ApproxMatch','SBStd','SBexact','SBapprox','FPStd','FPexact', 'FPapprox', 'Noisy img','Reference img'};


nOfTests=size(DenoiserType,2);


%lambdaSP=9*10^(-2);
lambdaSP=4.95*10^(-2);
lambdaL1= 2.5; %buono compreso tra 2 e 2,9, con 1,9 o 3,1 comportamenti diversi... malissimo per 10^(-3) a scendere 
Theta= (2/7)*10^(-1);
Tol=10^(-3);
lambda=lambdaSP;
fprintf('lamb = %f \n', lambda)
fprintf('Tol = %f \n', Tol)


%% Compute and print


for i=1:nOfTests
    HistMatchType=string(MatchType(i));
    ImageDenoiserType=string(DenoiserType(i));
    if strcmp(Str(i),'Reference img')
        [u{i},k{i},rel{i},l{i}, J{i}] = SB_ROF_VARIABLE_2stage(ImageToMatch,lambdaSP,lambdaL1,mu,Theta, ImageToMatch,ImageDenoiserType, HistMatchType,Tol, kMax, X_true);
        [outStr{i},nEr(i),Psnr(i),Ssim(i)]=printErPnsImages(string(Str{i}), u{i}, X_true);
        fprintf(string(outStr{i}));
    else
        [u{i},k{i},rel{i},l{i}, J{i}] = SB_ROF_VARIABLE_2stage(X_noise,lambdaSP,lambdaL1,mu,Theta, ImageToMatch,ImageDenoiserType, HistMatchType,Tol, kMax, X_true);
        [outStr{i},nEr(i),Psnr(i),Ssim(i)]=printErPnsImages(string(Str{i}), u{i}, X_true);
        fprintf(string(outStr{i}));
    end
end


 end_nbins=17;
 [figure_Comparison_Nbins]= ComparisonGraph_Nbins(end_nbins,X_noise, Theta, Tol, kMax, X_true,lambdaSP,mu,ImageToMatch);
 exportgraphics(figure_Comparison_Nbins, fullfile(outImage, strcat('graph_Comparison_Nbins','.jpg')),"resolution",300)
 exportgraphics(figure_Comparison_Nbins, fullfile(outImage, strcat('graph_Comparison_Nbins','.png')),"resolution",300)



% % Display
fig1=figure('Name', ['Display original noisy' 'Noise'  num2str(s) ] ); colormap gray;
subaxis(1,2,1,'SpacingVert',0,'MR',0); imagesc(X_true); axis off image; title('Original','FontSize',20); zoom(0);
subaxis(1,2,2,'SpacingVert',0,'MR',0); imagesc(X_noise); axis off image; title('Noisy','FontSize',20); zoom(0);

exportgraphics(fig1, fullfile(outImage, strcat('TrueNoise', num2str(s),'.jpg')),"resolution",300)
exportgraphics(fig1, fullfile(outImage, strcat('TrueNoise', num2str(s),'.png')),"resolution",300 )


figT=zoom_Function(X_true,50,40,70,40, 'True');
exportgraphics(figT,fullfile(outImage, strcat('True','.jpg')))
exportgraphics(figT,fullfile(outImage, strcat('True','.png')), "resolution",300 )

figN=zoom_Function(X_noise,50,40,70,40, [ 'Noise'  num2str(s)] );
exportgraphics(figN,fullfile(outImage, strcat('Noise', num2str(s),'.jpg')),"resolution",300)
exportgraphics(figN,fullfile(outImage, strcat('Noise', num2str(s),'.png')),"resolution",300 )

figTm=zoom_Function(ImageToMatch,50,40,70,40, [ 'ImageToMatch'  num2str(s)] );
exportgraphics(figTm,fullfile(outImage, strcat('ImageToMatch','.jpg')), "resolution",300)
exportgraphics(figTm,fullfile(outImage, strcat('ImageToMatch','.png')), "resolution",300 )

figCom=figure('Name', ['Comparison among denoiser' 'Noise'  num2str(s)] ); colormap gray;

for i=1:nOfTests
    subaxis(2,3,i,'SpacingVert',0.03,'MR',0.001); imagesc(u{i}); axis image off; title(string(Str{i})); zoom(0);
end
tightfig;

exportgraphics(figCom,fullfile(outImage, strcat('ComparisonAmongDenoiserNoise', num2str(s) ,'.jpg')),"resolution",300)
exportgraphics(figCom,fullfile(outImage, strcat('ComparisonAmongDenoiserNoise', num2str(s) ,'.png')),"resolution",300 )

for i=1:nOfTests
    str=strcat('Noise', num2str(s) ,'Denoiser', char(Str{i}),'.jpg');
    str1=strcat('Noise', num2str(s) ,'Denoiser', char(Str{i}),'.png');
    fls=fullfile(outImage, str);
    fls1=fullfile(outImage, str1);
    figCom=zoom_Function(u{i},50,40,70,40, ['Denoiser' char(Str{i})  'Noise'  num2str(s) ] );
   
    exportgraphics(figCom,fls, "resolution",300)
    exportgraphics(figCom,fls1, "resolution",300 )
    Histo=Histogram2d(X_true,u{i},  'X True',strcat('Noise', num2str(s) ,' Denoiser', char(Str{i})));
    exportgraphics(Histo, fullfile(outImage, strcat('H2d_True_','Noise', num2str(s) ,'Denoiser', char(Str{i}), num2str(s),'.jpg')),"resolution",300)
    exportgraphics(Histo, fullfile(outImage, strcat('H2d_True_','Noise', num2str(s) ,'Denoiser', char(Str{i}), num2str(s),'.png')),"resolution",300)

end



%Fig 2
f2=figure('Name',['Error vs. iteration' 'Noise'  num2str(s)]);
i=(0:k{1}); semilogy(i, l{1}(i+1)); hold on;
for j=2:nOfTests
    i=(0:k{j}); semilogy(i, l{j}(i+1),'--')
end
xlabel('Number of Iterations','FontSize',20)
ylabel('Error','FontSize',20)
legend(Str, 'Location', 'best')
legend('boxoff')
exportgraphics(f2,fullfile(outImage, strcat('Noise', num2str(s), 'ErrorVsIteration','.jpg')), "resolution",300)
exportgraphics(f2,fullfile(outImage, strcat('Noise', num2str(s), 'ErrorVsIteration','.png')), "resolution",300 )

%Fig 4
f4=figure('Name',['Rel. Error vs. iteration' 'Noise'  num2str(s,'%1g')]);
i=(0:k{1}); semilogy(i, rel{1}(i+1)); hold on;
for j=2:nOfTests
    i=(0:k{j}); semilogy(i, rel{j}(i+1),'-.')
end
xlabel('Number of Iterations','FontSize',20)
ylabel('Relative Error','FontSize',20)
legend(Str, 'Location', 'best')
legend('boxoff')
exportgraphics(f4,fullfile(outImage, strcat('Noise', num2str(s),'RelErVsIteration','.jpg')), "resolution",300)
exportgraphics(f4,fullfile(outImage, strcat('Noise', num2str(s),'RelErVsIteration','.png')), "resolution",300 )


%Fig 5
f5=figure('Name',['Comparison of Rel. Err., PDNR and SSIM' 'Noise'  num2str(s)]);
c = categorical(Str);
subplot(3,1,1);
errors = nEr;
bar(c,errors,'BaseValue',min(errors)),set(gca,'yscale','log')
title('Err','FontSize',20);

subplot(3,1,2);
Psnrs = Psnr;
bar(c,Psnrs,'BaseValue',min(Psnrs))
title('PSNR','FontSize',20);

subplot(3,1,3);
bar(c,Ssim,'BaseValue',min(Ssim))
title('SSIM','FontSize',20);

fig5Nmae=fullfile(outImage, strcat('Noise', num2str(s), 'PDNRbar','.jpg'));
exportgraphics(f5, fig5Nmae,"resolution",300)
exportgraphics(f5, fullfile(outImage, strcat('Noise', num2str(s), 'PDNRbar','.png')),"resolution",300 )

close all
%% Print table on Latex

k=cell2mat(k);
T = table(nEr(:),Psnr(:),Ssim(:),k(:), 'RowNames',Str,'VariableNames',{'RelErr' 'PSNR' 'SSIM' 'numOfIter'} );
% Now use this table as input in our input struct:
input.data = T;
% Set the row format of the data values 
input.dataFormat = {'%.4f',3,'%d',1};
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'l';
% Switch table borders on/off:
input.tableBorders = 1;
% Transpose Table:
input.transposeTable = 0;
% % LaTex table caption:
% input.tableCaption = 'MyTableCaption';
% % LaTex table label:
% input.tableLabel = 'MyTableLabel';
T
% Now call the function to generate LaTex code:
    latex = latexTable(input);


% % save LaTex code as file
% OutLatex='OutLatex/VerMyLatexBasicNoise';
% fid=fopen(strcat(OutLatex, num2str(s),'.tex'),'w');
% [nrows,ncols] = size(latex);
% for row = 1:nrows
%     fprintf(fid,'%s\n',latex{row,:});
% end
% fclose(fid);
% fprintf('\n... your LaTex code has been saved as ''MyLatex_____.tex'' in your working directory\n');

%figure(3)

%plot(X(100,:))
%hold on;
%plot(u3(100,:))
%xlabel('Distance along profile')
%ylabel('Pixel intensity (Gray Value)')
%title('Intensity Profile (Isotropic TV-LSE)')
%legend('Original Image','Restored Image', 'Location', 'south')
%legend('boxoff')

%miglioramenti percentuali
[fp_exact,fp_approx,sb_exact,sb_approx]=miglioramenti(nEr,Psnr,Ssim);
row={'RelErr' 'PSNR' 'SSIM'};
Str2= {'Mig% FP_exact vs FP-Std' 'Mig% FP_approx vs FP-Std' 'Mig% SB_exact vs SB-Std' 'Mig% SB_approx vs SB-Std'};
T2=table(fp_exact(:),fp_approx(:),sb_exact(:),sb_approx(:),'RowNames',row,'VariableNames',Str2 );
T2 
