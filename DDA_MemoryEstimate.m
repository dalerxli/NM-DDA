function [ memoryRequired , GPUmemoryRequired ] = DDA_MemoryEstimate( sizeArray )
%DDA_MEMORYCOMPUTE 
%�˺������ڹ����Model ��Ҫ���ڴ��Լ��Դ� 
%ʹ�ô˺���ʱ�� ����Ϊ size( Model.struc )
Nx = sizeArray(1) - 1 ;
Ny = sizeArray(2) - 1 ;
Nz = sizeArray(3) - 1 ;

%�ܸ����
N = Nx * Ny * Nz ;

%matlab��һЩ������Ҫ�����ڴ棬��������ϵ����ʵ������ʱ��õġ�
realRatioCPU = 2 ;
realRatioGPU = 3 ;

% memoryPeakRequired ����Ϊ��DDA_Compute�����ڣ�������A��Ȼ��Ŵ���GPU
GPUmemoryRequired = realRatioGPU * 64 * N * 4 / 1024^2 ;
memoryRequired = realRatioCPU * 22 * N * 4 / 1024^2 ;
memoryPeakRequired = 70 * N * 4 / 1024^2 ;
memoryPerOrienRequired = 6 * N * 4 / 1024^2 ;

fprintf('��Ҫ�Ļ����ڴ�Ϊ %7.1f MB����Ҫ���Դ�Ϊ %7.1f MB��\n' ,memoryRequired ,GPUmemoryRequired ) ;
fprintf('��Ҫ���ڴ��ֵΪ %7.1f MB�� \n' , memoryPeakRequired ) ;
fprintf('ÿ���һ���ض�ȡ��ļ��㣬�����ڴ����� %7.1f MB�� \n' ,memoryPerOrienRequired) ;
end

