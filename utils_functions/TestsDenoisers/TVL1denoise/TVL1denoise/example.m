% TVL1denoise demo
%
% Manolis Lourakis 2016

im=phantom;
imshow(im);

outim=TVL1denoise(im, 1.0, 100);
figure; imshow(outim, []);
