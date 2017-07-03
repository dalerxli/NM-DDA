function  Model  = DDA_ModelGenerate( shape , diaArray , m , accuracy)
%函数生成一个Model
%Model结构体完全地描述了一个粒子的形状、大小、折射率、以及对应的空间点阵

%输入说明
% shape = 'cylinder' 'ellipsoid' 'brick'（若想输出球，则使得椭球三轴相等即可）
% diaArray 粒子的X Y Z方向上的线度，单位为um 比如, 'ellipsoid' [ 1 2 1] 代表一个直径分别为 1 2 1
% um的椭球 ， 'cylinder' [1 2 1] 代表一个高 2um ，直径 1 um 的圆柱
% m 粒子的复折射率
% 注意！！旋转模块的轴向为Y轴正方向。
% 轴向的确定与DDA_TargetOrienGener函数中默认的轴向有关，不选用此轴向会在取向选取不为均匀球面时出现问题
% |m|kd的值会影响计算精度 如果希望高精度计算，可以通过设置 accuracy 完成
% accuracy = 'Normal' 'High' 'Ultra' 对应的|m|kd 分别为 0.8 0.5 0.25 。缺省值为 'Normal'

%输出说明
%Model 结构体
%Model.shape 说明model的形状，字符串格式
%Model.aeff  说明model的等效粒径
%Model.m     说明model的复折射率
%Model.struc 一个三维矩阵，用于标定哪些点有偶极子
%Model.N     说明model中偶极子的个数
%Model.d     说明model中偶极子之间的距离
%Model.fileName 为了方便保存而生成的名字

%相比于 DDA_ModelGenerate_aeff ，此函数在于让用户自由输入粒子空间线度，并自适应选择d，生成空间点阵

%关于struc矩阵
%struc的长宽高为Nx , Ny , Nz
%struc的三个方向的节点数是Nx + 1 , Ny + 1 , Nz + 1
%所以 struc是一个size为( Nx + 1 , Ny + 1 , Nz + 1 ) 的零一矩阵
%而模型的长宽高为 Nx * d ，Ny * d ，Nz * d

%以下为生成Narray 此函数会从Narray中挑选可以包含模型的最小Nx Ny Nz来输出模型
twoNP1Array = [ 27 45 63 75 81 99 105 117 125 135 147 165 175 189 195 225 231 243 245 ...
    273 275 297 315 325 343 351 363 375 385 405 429 441 455 495 507 525 539 567 585 605 ...
    625 637 675 693 715 729 735 819 825 845];
Narray = ( twoNP1Array - 1 ) / 2 ;

% 默认lambda
lambda = 0.532 ;
k = 2 * pi / lambda ;

% 以下部分开始生成Model
Model.m = m ;
Model.shape = shape ;

if nargin < 4
    accuracy = 'Normal' ;
end

% 得到d
switch accuracy
    case 'Normal'
        d = 0.8 / k / abs( m ) ;
    case 'High'
        d = 0.5 / k / abs( m ) ;
    case 'Ultra'
        d = 0.25 / k / abs( m ) ;
end

%开始生成 ...
diax = diaArray(1) ;
diay = diaArray(2) ;
diaz = diaArray(3) ;

%椭球的情况
if strcmp( shape , 'ellipsoid')
    nx = ceil( diax / d ) ;
    ny = ceil( diay / d ) ;
    nz = ceil( diaz / d ) ;
    
    Nx = Narray(find(Narray > nx + 2 , 1)) ;
    Ny = Narray(find(Narray > ny + 2 , 1)) ;
    Nz = Narray(find(Narray > nz + 2 , 1)) ;
    
    Model.struc = zeros( Nx + 1 , Ny + 1 , Nz + 1 ,'logical') ;
    
    xc = round( (Nx + 2) / 2 ) ;
    yc = round( (Ny + 2) / 2 ) ;
    zc = round( (Nz + 2) / 2 ) ;
    
    for iterx = 1 : Nx + 1
        for itery = 1 : Ny + 1
            for iterz = 1 : Nz + 1
                Model.struc(iterx,itery,iterz) = ( xc - iterx)^2 / nx^2 + ( yc - itery )^2 / ny^2 ...
                    + ( zc - iterz )^2 / nz^2 <= 0.25 ;
            end
        end
    end
    
    %计算偶极子个数
    pointNum = sum( Model.struc(:) ) ;
    
    %重新生成d 以对得上体积
    d = ((4/3 * pi * diax * diay * diaz / 8) / pointNum )^(1/3) ;
    
    %Model.fileName
    dxName = num2str( diax );
    dyName = num2str( diay );
    dzName = num2str( diaz );
    
    if length( dxName ) >= 4
    dxName = dxName( 1:4 ) ;
    end 
    if length( dyName ) >= 4
    dyName = dyName( 1:4 ) ;
    end
    if length( dzName ) >= 4
    dzName = dzName( 1:4 ) ;
    end
    Model.fileName = ['ELLI_',dxName,'X',dyName,'X',dzName ,'_m',num2str(m)] ;
end

%圆柱的情况
if strcmp( shape , 'cylinder')
    nx = ceil( diax / d ) ;
    ny = ceil( diay / d ) ;
    nz = ceil( diaz / d ) ;
    
    Nx = Narray(find(Narray > nx + 2 , 1)) ;
    Ny = Narray(find(Narray > ny + 2 , 1)) ;
    Nz = Narray(find(Narray > nz + 2 , 1)) ;
    
    Model.struc = zeros( Nx + 1 , Ny + 1 , Nz + 1 ,'logical') ;
    
    xc = round( (Nx + 2) / 2 ) ;
    yc = round( (Ny + 2) / 2 ) ;
    zc = round( (Nz + 2) / 2 ) ;
    
    for iterx = 1 : Nx + 1
        for itery = 1 : Ny + 1
            for iterz = 1 : Nz + 1
                Model.struc(iterx,itery,iterz) = ( xc - iterx)^2 / nx^2 + ( zc - iterz )^2 / nz^2 <= 0.25 ...
                    && abs(yc - itery) <= abs( ny / 2 ) ;
            end
        end
    end
    
    %计算偶极子个数
    pointNum = sum( Model.struc(:) ) ;
    
    %重新生成d 以对得上体积
    d = ((pi * diax * diay * diaz / 4) / pointNum )^(1/3) ;
    
    %Model.fileName
    dxName = num2str( diax );
    dyName = num2str( diay );
    dzName = num2str( diaz );
    if length( dxName ) >= 4
    dxName = dxName( 1:4 ) ;
    end 
    if length( dyName ) >= 4
    dyName = dyName( 1:4 ) ;
    end
    if length( dzName ) >= 4
    dzName = dzName( 1:4 ) ;
    end
    Model.fileName = ['CYLIN_',dxName,'X',dyName,'X',dzName ,'_m',num2str(m)] ;
end

%方块的情况
if strcmp( shape , 'brick')
    
    nx = ceil( diax / d ) ;
    ny = ceil( diay / d ) ;
    nz = ceil( diaz / d ) ;
    
    Model.struc(1:nx , 1: ny , 1:nz) = 1;
    
    %计算偶极子个数
    pointNum = sum( Model.struc(:) ) ;
    
    %重新生成d 以对得上体积
    d = (( diax * diay * diaz ) / pointNum )^(1/3) ;
    
    %Model.fileName
    dxName = num2str( diax );
    dyName = num2str( diay );
    dzName = num2str( diaz );
    
    if length( dxName ) >= 4
    dxName = dxName( 1:4 ) ;
    end 
    if length( dyName ) >= 4
    dyName = dyName( 1:4 ) ;
    end
    if length( dzName ) >= 4
    dzName = dzName( 1:4 ) ;
    end
    Model.fileName = ['BRICK_',dxName,'X',dyName,'X',dzName ,'_m',num2str(m)] ;
end

%aeff和偶极子数N
Model.N = pointNum ;
Model.aeff = ( pointNum * d^3 / pi / 4 * 3 )^(1/3) ;
Model.d = d ;

strucFigureHandle = DDA_StrucPlot( Model ) ;
[ memoryRequired , GPUmemoryRequired ] = DDA_MemoryEstimate( size( Model.struc ) ) ;

save(['Model_', Model.fileName ,'.mat'],'Model') ;
end



