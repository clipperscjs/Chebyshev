function [outimg]=padImg(img0)
[M,N]=size(img0);
%����ͼ��Ϊ����*����
if rem(M,2)==0 && rem(N,2)==1 %��Ϊż��
   img=padarray(img0,[1 0],'replicate','post');
end
if rem(M,2)==1 && rem(N,2)==0 %��Ϊż��
   img=padarray(img0,[0 1],'replicate','post');
end
if rem(M,2)==0 && rem(N,2)==0 %��Ϊż��
   img=padarray(img0,[1 1],'replicate','post');
end
if rem(M,2)==1 && rem(N,2)==1 %��Ϊ����
    img=img0;
end
outimg=img;