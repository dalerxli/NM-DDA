function  outArray = DDA_ConvAccelerate( Afg , Y , AHermiFlag)

%此函数用于  矩阵Afg 乘以 向量Y
% Afg 已经被fft过，且已经储存在显存之中
% Y 为原向量
% 此算法的核心是用fft加速卷积运算 ， 详细理论参见 Goodman+Draine+Flatau 1991 年的文章，以及wiki上的卷积定理

% A的大小为 2Nx + 1, 2Ny + 1 ,2Nz + 1
% Y的大小为 Nx + 1, Ny + 1 ,Nz + 1
% 后缀g代表是gpuArray

%创建一个和Afg一样大的Yg阵 outArrayg阵
persistent Yg ;
persistent outArrayg ;
if isempty(Yg)
    Yg = gpuArray.zeros(size(Afg),'single') ;
end
if isempty(outArrayg)
    outArrayg = gpuArray.zeros(size(Afg),'single') ;
end
%将Y的值赋值给Yg 并使Yg其他的地方为0 注意Yg几乎有 7/8 的地方为0元素
n1 = size(Y,1) ;
n2 = size(Y,2) ;
n3 = size(Y,3) ;

Yg(1:n1 , 1:n2 , 1:n3 ) = Y;

%使用卷积定理 并提取结果
%AHermiFlag 是应用于这样一种情况：由于A是Hermi的，而算法中需要用到A'的乘法，故而让A直接乘以Y' 并对结果取复共轭
%详见上文提到的paper

if AHermiFlag == 0
    outArrayg = ifftn(Afg .* fftn(Yg)) ;
    outArray = gather(outArrayg(n1 : 2*n1 - 1 , n2 : 2*n2 - 1 , n3 : 2*n3 - 1)) ;
end

if AHermiFlag == 1
    outArrayg = conj( ifftn( Afg .* fftn( conj( Yg ))) );
    outArray = gather(outArrayg(n1 : 2*n1 - 1 , n2 : 2*n2 - 1 , n3 : 2*n3 - 1)) ;
end


