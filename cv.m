function CVvalue = cv(img0)
%-------��ֵ---------
img=im2double(img0);%����ͼ��
mean_img=mean2(img);
%------��׼��---------
std_img=std2(img);

CVvalue=std_img/mean_img;
