function [ psnr ] = Mypsnr( I,ori )
%MYPSNR 此处显示有关此函数的摘要
%   此处显示详细说明
%%%%%%计算psf
clearImage=ori;
mc = max(max(clearImage));
mcm = min(min(clearImage));
clearImage1 = (clearImage-mcm)/(mc-mcm)*255;
% 
% BlurImage=I;
% mc = max(max(BlurImage));
% mcm = min(min(BlurImage));
% BlurImage1 = (BlurImage-mcm)/(mc-mcm)*255;

DeblurImage=I;
mc = max(max(DeblurImage));
mcm = min(min(DeblurImage));
DeblurImage1 = (DeblurImage-mcm)/(mc-mcm)*255;
[row,col] = size(DeblurImage);
%Obs = - 20*log10(sqrt(norm(clearImage1- BlurImage1, 'fro')^2/(row*col*255^2)));
psnr = - 20*log10(sqrt(norm(DeblurImage1- clearImage1, 'fro')^2/(row*col*255^2)));
% ss = sprintf('PSNR: %.4f dB, Restored PSNR: u %.4f dB\n',Obs,Res);

% disp(ss);
end

