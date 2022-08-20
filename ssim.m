function [ssim_val] = ssim(I,I1)
 %¶ÁÈëÍ¼Ïñ
if size(I,3)==3
    I= im2double(rgb2gray(I));
else
    I= im2double(I);
end
Imax = max(max(I));
Imin = min(min(I));
IN = (I-Imin)/(Imax-Imin)*255;
 %¶ÁÈëÇåÎúÍ¼Ïñ
if size(I1,3)==3
    I1= im2double(rgb2gray(I1));
else
    I1= im2double(I1);
end
I1max = max(max(I1));
I1min = min(min(I1));
IN1 = (I1-I1min)/(I1max-I1min)*255;
 K = [0.05 0.05];
window = ones(8);
ssim_val=ssim_index_new(IN,IN1, K, window);

