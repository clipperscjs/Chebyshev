function [b]= remove_by_Wa(y,filtered_im)
[N,M]=size(y);
b = filtered_im;
K=11;
C=M;
B=N;
[W]=Get_chebeyshev(K,C,B);
a = (W'*W)\(W'*b(:));
Wa0 =W*a;
Wa=zeros(B,C);
for i=1:C
     Wa(:,i)=Wa0(B*(i-1)+1:B*i,1);
end
b=Wa;
