function [str,varargout] = printErPnsImages(strPrs,X_out,X_true)

p = norm(X_out-X_true,2)/norm(X_true,2);

% psnr1= 10*log10(   (max(X_true(:)))^2   /   immse(X_true,X_out)   )
% peakval=sqrt( immse(X_true,X_out)*10^(Psnr/10) )
% massimo=(max(X_true(:)))

Psnr=psnr(X_out,X_true);

Ssim= ssim(X_out,X_true);

str3=sprintf(strcat(strPrs, ' SSIM = %.7f \n'), Ssim );

str2=sprintf(strcat(strPrs, ' PSNR = %.7f \n'), Psnr);

str1=sprintf(strcat(strPrs, ' Relative Err = %g \n'), p);

varargout{1} = p;
varargout{2} = Psnr;
varargout{3} = Ssim;

str= strcat(str1, str2,str3);


end