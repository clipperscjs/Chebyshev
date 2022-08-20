function  [output]=NLmeansfilter(input,t,f,h)                
[m n]=size(input); 
%t:搜索框半径；f:相似框半径；h:滤波频率百分比
Output=zeros(m,n);
input2 = padarray(input,[f f],'symmetric');              
 %将边缘对称折叠上去  
 % f :加宽的宽度值
 kernel = make_kernel(f);                              %计算得到一个高斯核，用于后续的计算
 kernel = kernel / sum(sum(kernel));             %sum：对矩阵k的每一列求和，k表示矩阵k的总和
 h=h*h;
 for i=1:m
 for j=1:n
         i1 = i+ f;
         j1 = j+ f;      
         W1= input2(i1-f:i1+f , j1-f:j1+f);		
         wmax=0; 
         average=0;
         sweight=0;
         rmin = max(i1-t,f+1);      %确定相似框的边长，其实可以在代入数据的时候就设置搜索框和边界框边长的大小关系，就不用求最大最小值了
         rmax = min(i1+t,m+f);      %这段主要是控制不超出索引值
         smin = max(j1-t,f+1);
         smax = min(j1+t,n+f);
         for r=rmin:1:rmax
         for s=smin:1:smax                 
                if(r==i1 && s==j1) continue; end;      
                W2= input2(r-f:r+f , s-f:s+f);                
                d = sum(sum(kernel.*(W1-W2).*(W1-W2)));
                w=exp(-d/h);                 
                if w>wmax      
                   wmax=w;                   
                end
                sweight = sweight + w;
                average = average + w*input2(r,s);                                  
         end 
         end
        average = average + wmax*input2(i1,j1);
        sweight = sweight + wmax;       
        if sweight > 0
            output(i,j) = average / sweight;
        else
            output(i,j) = input(i,j);
        end                
 end
 end	

 
 
 function [kernel] = make_kernel(f)            %这是一个高斯核   
kernel=zeros(2*f+1,2*f+1);   
for d=1:f    
  value= 1 / (2*d+1)^2 ;    
  for i=-d:d
  for j=-d:d
    kernel(f+1-i,f+1-j)= kernel(f+1-i,f+1-j) + value ;
  end
  end
end
kernel = kernel ./ f;    

