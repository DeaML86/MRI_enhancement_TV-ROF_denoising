function [str,varargout] = printErPnsImages(strPrs,X_out,X_true)

p = norm(X_out-X_true,2)/norm(X_true,2);
psnr1= 10*log10((max(X_true(:)))^2/immse(X_true(:),X_out(:)));

str1=sprintf(strcat(strPrs, ' Relative Err. = %g \n'), p);
str2=sprintf(strcat(strPrs, ' PSNR = %.7f \n'), psnr1);


varargout{1} = p;
varargout{2} = psnr1;

str= strcat(str1, str2);

end