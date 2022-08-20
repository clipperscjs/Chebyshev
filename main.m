clc;
clear all;
close all;

%% degrade image
addpath('image_degraded\');
file_path0 =  '.\image_degraded\';
img_path_list0 = dir(strcat(file_path0,'*.bmp'));
img_path_list0_name = sort_nat({img_path_list0.name});
%% clear image
addpath('image_clear\')
file_path1 =  '.\image_clear\';
img_path_list1 = dir(strcat(file_path1,'*.bmp'));
img_path_list1_name = sort_nat({img_path_list1.name});
img_num = length(img_path_list0);
%%
for j = 1
    image_name0 = img_path_list0_name{j};
    fprintf('%d %s\n',j,strcat(file_path0,image_name0));
    yg0 = im2double(imread(strcat(file_path0,image_name0)));
    yg = padImg(yg0);
    % yg=rgb2gray(yg); if RGB

    image_name1 = img_path_list1_name{j};
    fprintf('%d %s\n',j,strcat(file_path1,image_name1));
    ycimage0 = im2double(imread(strcat(file_path1,image_name1)));
    ycimage = padImg(ycimage0);
    %     ycimage=rgb2gray(ycimage); if RGB
    %% Multi-scale
    [M N]=size(yg);
    M1=floor(M*(1/2)^(1/2));N1=floor(N*(1/2)^(1/2));
    M2=floor(M*(1/2));N2=floor(N*(1/2));
    yg1_0=imresize(yg,[M1 N1]);
    yg1=padImg(yg1_0);
    yg2_0=imresize(yg,[M2 N2]);
    yg2=padImg(yg2_0);

    ycimage1_0=imresize(ycimage,[M1 N1]);
    ycimage1=padImg(ycimage1_0);
    ycimage2_0=imresize(ycimage,[M2 N2]);
    ycimage2=padImg(ycimage2_0);

    yfilter= wlsFilter(yg);
    yfilter1 = wlsFilter(yg1);
    yfilter2 = wlsFilter(yg2);

    %NLM
    %t:������뾶��f:���ƿ�뾶��h:�˲�Ƶ�ʰٷֱ�,function  [output]=NLmeansfilter(input,t,f,h)
%     yfilter = NLmeansfilter(yg,5,2,10);
%     yfilter1=NLmeansfilter(yg1,5,2,10);
%     yfilter2=NLmeansfilter(yg2,5,2,10);

    %% iter
    tic
    K=5;
    [M1_0,N1_0]=size(yg1);
    [s,b] = biasCorrection_filtered_gradient_l0l2norm(yg2,yfilter2,ycimage2,K);
    s1=imresize(s,[M1_0 N1_0]);b1=imresize(b,[M1_0 N1_0]);
    [s,b] = Muti_scal_1(yg1,s1,b1,ycimage1,K);
    s0=imresize(s,[M N]);b0=imresize(b,[M N]);
    [s,b] = Muti_scal_0(yg,s0,b0,ycimage,K);
    toc
    % [b]=remove_by_Wa(yg,yfilter);

    %% Norm
    clearImg0 = yg-b;
    maxv = max(max(clearImg0));
    minv = min(min(clearImg0));
    clearImg = (clearImg0-minv)/(maxv-minv);

    bias0 = b;
    maxb = max(max(bias0));
    minb = min(min(bias0));
    bias = (bias0-minb)/(maxb-minb);

    % figure(1);
    % imshow(clearImg);
    % figure(2);
    % imshow(yg);
    % figure(3);
    % imshow(b,[]);
    % figure(4);
    % imshow(yfilter);

    fprintf('�˻�ͼPSNR=%4.2f\n',Mypsnr(yg,ycimage));
    fprintf('�˻�ͼSSIM=%6.4f\n',ssim(yg,ycimage));
    fprintf('PSNR=%4.2f\n',Mypsnr(yg-b,ycimage));
    fprintf('SSIM=%6.4f\n',ssim(yg-b,ycimage));
    %     imwrite(yg,strcat('D:\����\ʵ��ԱȽ��\ours\92\',num2str(j),'.bmp'),'bmp');
    %     imwrite(ycimage,strcat('D:\����\ʵ��ԱȽ��\ours\92\',num2str(j),'_clear.bmp'),'bmp');
    %     imwrite(clearImg,strcat('D:\����\ʵ��ԱȽ��\ours\92\',num2str(j),'_result.bmp'),'bmp');
    %     imwrite(clearImg,strcat('D:\����\ʵ��ԱȽ��\real_image\Ours\',num2str(j),'_result.jpg'),'jpg');
    %     imwrite(bias,strcat('D:\����\ƫ�ù���ʵ��\chebyshev_bias\k=9\202\',num2str(j),'_bias.bmp'),'bmp');
    %     fprintf('DE=%6.4f\n',DE(yg-b));
    %     fprintf('CV=%6.4f\n',cv(yg-b));
    %

end


