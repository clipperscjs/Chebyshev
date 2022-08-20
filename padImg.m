function [outimg]=padImg(img0)
[M,N]=size(img0);
%扩充图像为奇数*奇数
if rem(M,2)==0 && rem(N,2)==1 %行为偶数
   img=padarray(img0,[1 0],'replicate','post');
end
if rem(M,2)==1 && rem(N,2)==0 %列为偶数
   img=padarray(img0,[0 1],'replicate','post');
end
if rem(M,2)==0 && rem(N,2)==0 %均为偶数
   img=padarray(img0,[1 1],'replicate','post');
end
if rem(M,2)==1 && rem(N,2)==1 %均为奇数
    img=img0;
end
outimg=img;