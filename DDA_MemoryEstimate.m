function [ memoryRequired , GPUmemoryRequired ] = DDA_MemoryEstimate( sizeArray )
%DDA_MEMORYCOMPUTE 
%此函数用于估算此Model 需要的内存以及显存 
%使用此函数时， 输入为 size( Model.struc )
Nx = sizeArray(1) - 1 ;
Ny = sizeArray(2) - 1 ;
Nz = sizeArray(3) - 1 ;

%总格点数
N = Nx * Ny * Nz ;

%matlab的一些运算需要额外内存，下面两个系数是实际运行时测得的。
realRatioCPU = 2 ;
realRatioGPU = 3 ;

% memoryPeakRequired 是因为在DDA_Compute函数内，先生成A，然后才存入GPU
GPUmemoryRequired = realRatioGPU * 64 * N * 4 / 1024^2 ;
memoryRequired = realRatioCPU * 22 * N * 4 / 1024^2 ;
memoryPeakRequired = 70 * N * 4 / 1024^2 ;
memoryPerOrienRequired = 6 * N * 4 / 1024^2 ;

fprintf('需要的基本内存为 %7.1f MB，需要的显存为 %7.1f MB。\n' ,memoryRequired ,GPUmemoryRequired ) ;
fprintf('需要的内存峰值为 %7.1f MB。 \n' , memoryPeakRequired ) ;
fprintf('每完成一个特定取向的计算，所需内存增加 %7.1f MB。 \n' ,memoryPerOrienRequired) ;
end

