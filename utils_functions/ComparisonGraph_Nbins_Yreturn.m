function [y_Psnr_SB,y_Psnr_FP,y_RelErr_SB,y_RelErr_FP]= ComparisonGraph_Nbins_Yreturn(end_nbins, Theta, Tol, kMax,lambdaSP,mu,ref,s,shift,slice)

SetParamScript;
addpath(genpath('TestsDenoisers'))

% Set up test image
pulseSequence='t1';
file=fullfile(imagePath, pulseSequence, 't1_icbm_normal_1mm_pn0_rf20.mnc');
mriData = loadminc(file);
fprintf('Dimension MRI %dx%dx%d\n', size(mriData))
X_true = mriData(:,:,slice); % (particular slice through the brain)

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

for i=1:end_nbins
    nbins=2^i;
    i
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

y_Psnr_SB=cell2mat(graph_Psnr_SB);
y_Psnr_FP=cell2mat(graph_Psnr_FP);

y_RelErr_SB=cell2mat(graph_relErr_SB);
y_RelErr_FP=cell2mat(graph_relErr_FP);

