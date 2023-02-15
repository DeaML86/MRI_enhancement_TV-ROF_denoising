function [fp_exact,fp_approx,sb_exact,sb_approx]=miglioramenti(RelErr,Psnr,Ssim)

for i =1:3
    switch i

        case 1
            fpstd=RelErr(6);
            a=RelErr(7);
            b=RelErr(8);
           
            sbstd=RelErr(3);
            c=RelErr(4);
            d=RelErr(5);
            
        case 2
            fpstd=Psnr(6);
            a=Psnr(7);
            b=Psnr(8);
           
            sbstd=Psnr(3);
            c=Psnr(4);
            d=Psnr(5);
           
        case 3
            fpstd=Ssim(6);
            a=Ssim(7);
            b=Ssim(8);
           
            sbstd=Ssim(3);
            c=Ssim(4);
            d=Ssim(5);
    end


    fp_exact{i}=abs( (a-fpstd)/fpstd )*100;
    fp_approx{i}=abs( (b-fpstd)/fpstd )*100;
    sb_exact{i}=abs( (c-sbstd)/sbstd )*100;
    sb_approx{i}=abs( (d-sbstd)/sbstd )*100;

end
