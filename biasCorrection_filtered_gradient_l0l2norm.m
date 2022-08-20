function [s,b,a]=biasCorrection_filtered_gradient_l0l2norm(y,filtered_im,ycimage,K)


%迭代
fx = [1, -1];
fy = [1; -1];
[N,M] = size(y);
sizeI2D = [N,M];
otfFx = psf2otf(fx,sizeI2D);
otfFy = psf2otf(fy,sizeI2D);
Denormin1 = abs(otfFx).^2 + abs(otfFy ).^2;
iter = 1;  %迭代次数
% xk_iter=40;
%对b给初值
C=M;%列数对应x
B=N;%行数对应y
[W]=Get_chebeyshev(K,C,B);
b=filtered_im;
a = (W'*W)\(W'*b(:));
Wa0 =W*a;
Wa=zeros(B,C);
for i=1:C
    Wa(:,i)=Wa0(B*(i-1)+1:B*i,1);
end
% figure(2); imshow(Wa,[]);
b=Wa;

%对s给初值
s=y-Wa;
% figure(1); imshow(s,[]);

alpha = 0.1;
beta = 0.05;
gamma = 2000;
delte = 22;

u_x = real(ifft2(otfFx.*fft2(s)));  % u_x(y)  d的初始值
u_y = real(ifft2(otfFy.*fft2(s)));  % u_y(y)
t_x = zeros(N,M);   % t赋初值
t_y = zeros(N,M);






psnrS(1) = Mypsnr(y-b,ycimage);
psnrT(1) = Mypsnr(y-Wa,ycimage);
psnrK(1) = Mypsnr(s,ycimage);
while 1
    % for iter=1:xk_iter
    %更新s
    Denormin = 1+(alpha+beta*delte)*Denormin1;  %  分母
    x1 = u_x - t_x;
    x2 = u_y - t_y;
    Ft2 = conj(otfFx).*fft2(x1) + conj(otfFy).*fft2(x2);
    Fs = (fft2(y)-fft2(b)+alpha*((fft2(y)-fft2(b)).*Denormin1)+beta*delte*Ft2)./Denormin;
    s = real(ifft2(Fs));
    
    
    %更新ux uy
    Du_x = real(ifft2(otfFx.*Fs));%deltax*s
    Du_y = real(ifft2(otfFy.*Fs));
    %     u_x = Du_x + t_x; u_y = Du_y + t_y;
    
    %0范数求解
    norm_0 = u_x.^2;
    ij = norm_0 >= 1/delte;
    nij = ~ij;
    u_x(ij) = u_x(ij);
    u_x(nij) = 0;
    
    norm_0 = u_y.^2;
    ij = norm_0 >= 1/delte;
    nij = ~ij;
    u_y(ij) = u_y(ij);
    u_y(nij) = 0;
    
    
    % Soft shrink formula 软阈值函数
    %     norm_d2 = sqrt(u_x.^2 + u_y.^2); ij = norm_d2 >= 15; nij = ~ij;
    %     u_x(ij) = u_x(ij) - 15*u_x(ij)./norm_d2(ij); u_x(nij) = 0; u_y(ij) =
    %     u_y(ij) - 15*u_y(ij)./norm_d2(ij); u_y(nij) = 0;
    
    % 更新tx ty
    t_x = t_x + Du_x - u_x;
    t_y = t_y + Du_y - u_y;
    % 更新b
    Denormin2 = 1+(alpha+gamma)*Denormin1;  %  分母
    Fb = (fft2(y)-fft2(s)+(alpha*(fft2(y)-fft2(s))+gamma*fft2(Wa)).*Denormin1)./Denormin2;
    b = real(ifft2(Fb));
    
    %更新a
    a = (W'*W)\(W'*b(:));
    Wa0=W*a;
    Wa=zeros(B,C);
    for i=1:C
        Wa(:,i)=Wa0(B*(i-1)+1:B*i,1);
    end
    psnrS(iter+1) = Mypsnr(y-b,ycimage);
    psnrT(iter+1) = Mypsnr(y-Wa,ycimage);
    psnrK(iter+1) = Mypsnr(s,ycimage);
    
    d=psnrS(iter+1)-psnrS(iter);
    iter = iter+1;
    if d<0&iter>=10
        break;
    end
    
end
% figure(8);
% plot(1:1:iter,psnrS,'b');grid on;
% hold on; plot(1:1:iter,psnrT,'r'); hold on; plot(1:1:iter,psnrK,'g');
% legend('y-b','y-Wa','s');
% xlabel('迭代次数');
% ylabel('PSNR');
% title('尺度2');

