function CVvalue = cv(img0)
%-------均值---------
img=im2double(img0);%读入图像
mean_img=mean2(img);
%------标准差---------
std_img=std2(img);

CVvalue=std_img/mean_img;
