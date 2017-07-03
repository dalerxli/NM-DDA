function  Model = DDA_ModelGenerate_aeff( shape , diaArray , aeff ,m )
%此函数是DDA_ModelGenerate的aeff版
%输入粒子的等效粒径的半径aeff、复折射率m 以及粒子的形状shape 以及三个方向形状对应的格点数diaArray
%函数会自动计算对应的晶格间距d
%当不满足 |m|kd < 1时， 函数会发出警告（ 默认的波长为 lambda = 0.532 ）
%此函数支持生成 球、椭球、柱
%20170418

%输入说明
% shape = 'cylinder' 'ellipsoid' 'brick'（若想输出球，则使得椭球三轴相等即可）
% diaArray 一个 1 X 3 数组 ，分别是 x y z方向上用于表示粒子形状的格点数 [40 40 40]表示一个直径为 40*d 的球
% （ 注意！！旋转模块的轴向为Y轴正方向。比如若输入一个长短轴比为1.5的椭球 diaArray 可以写为[40 60 40]
% 轴向的确定与DDA_TargetOrienGener函数中默认的轴向有关，不选用此轴向会在取向选取不为均匀球面时出现问题 ）
% aeff 颗粒物的等效半径 ， 单位为um
% m  物体的复折射率

%输出说明
%Model 结构体
%Model.shape 说明model的形状，字符串格式
%Model.aeff  说明model的等效粒径
%Model.m     说明model的复折射率
%Model.struc 一个三维矩阵，用于标定哪些点有偶极子
%Model.N     说明model中偶极子的个数
%Model.d     说明model中偶极子之间的距离

%关于struc矩阵
%struc的长宽高为Nx , Ny , Nz
%struc的三个方向的节点数是Nx + 1 , Ny + 1 , Nz + 1
%所以 struc是一个size为( Nx + 1 , Ny + 1 , Nz + 1 ) 的零一矩阵
%而模型的长宽高为 Nx * d ，Ny * d ，Nz * d

%重要！为了加快fft运算，注意生成的struc矩阵，其长宽高的2倍加1 应该是一个能高度分解质因子的数
%例如，Nx = 40    2 * Nx + 1 = 81   81 = 3^4
%相反 如果 Nx = 41     2 * Nx + 1 = 83   83 是质数，会极大增加fft的计算量

%以下为生成Narray 此函数会从Narray中挑选可以包含模型的最小Nx Ny Nz来输出模型
twoNP1Array = [ 27 45 63 75 81 99 105 117 125 135 147 165 175 189 195 225 231 243 245 ...
    273 275 297 315 325 343 351 363 375 385 405 429 441 455 495 507 525 539 567 585 605 ...
    625 637 675 693 715 729 735 819 825 845];
Narray = ( twoNP1Array - 1 ) / 2 ;

%默认lambda 用于检验 |m|kd < 1
lambda = 0.532 ;
k = 2 * pi / lambda ;

%以下部分开始生成Model
Model.m = m ;
Model.shape = shape ;
Model.aeff = aeff ;

diax = diaArray(1) ;
diay = diaArray(2) ;
diaz = diaArray(3) ;

Nx = Narray(find(Narray > diax + 2 , 1)) ;
Ny = Narray(find(Narray > diay + 2 , 1)) ;
Nz = Narray(find(Narray > diaz + 2 , 1)) ;

%先生成 Model.struc 的零阵，指定其size
Model.struc = zeros( Nx + 1 , Ny + 1 , Nz + 1 ,'logical') ;

%球、椭球的情况
if strcmp( shape , 'ellipsoid')
    xc = round( (Nx + 2) / 2 ) ;
    yc = round( (Ny + 2) / 2 ) ;
    zc = round( (Nz + 2) / 2 ) ;
    
    for iterx = 1 : Nx + 1
        for itery = 1 : Ny + 1
            for iterz = 1 : Nz + 1
                Model.struc(iterx,itery,iterz) = ( xc - iterx)^2 / diax^2 + ( yc - itery )^2 / diay^2 ...
                    + ( zc - iterz )^2 / diaz^2 <= 0.25 ;
            end
        end
    end
    
    %计算偶极子个数
    pointNum = sum( Model.struc(:) ) ;
    
    %计算偶极子间距d
    d = ( 4/3*pi * aeff^3 / pointNum )^(1/3) ;
    Model.d = d ;
    Model.N = pointNum ;
    
     %Model.fileName
    dxName = num2str( d * diax );
    dyName = num2str( d * diay );
    dzName = num2str( d * diaz );
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

%圆柱的情况 注意此时的圆柱是一个有限的圆柱 而不是无限长圆柱
if strcmp( shape , 'cylinder')
    xc = round( (Nx + 2) / 2 ) ;
    yc = round( (Ny + 2) / 2 ) ;
    zc = round( (Nz + 2) / 2 ) ;
    
    for iterx = 1 : Nx + 1
        for itery = 1 : Ny + 1
            for iterz = 1 : Nz + 1
                Model.struc(iterx,itery,iterz) = ( xc - iterx)^2 / diax^2 + ( zc - iterz )^2 / diaz^2 <= 0.25 ...
                    && abs(yc - itery) <= abs( diay / 2 ) ;
            end
        end
    end
    
    %计算偶极子个数
    pointNum = sum( Model.struc(:) ) ;
    
    %计算偶极子间距d
    d = ( 4/3*pi * aeff^3 / pointNum )^(1/3) ;
    Model.d = d ;
    Model.N = pointNum ;
    
    %Model.fileName
    dxName = num2str( d * diax );
    dyName = num2str( d * diay );
    dzName = num2str( d * diaz );
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
    Model.struc(1:diax , 1: diay , 1:diaz) = 1;
    
    %计算偶极子个数
    pointNum = sum( Model.struc(:) ) ;
    
    %计算偶极子间距d
    d = ( 4/3*pi * aeff^3 / pointNum )^(1/3) ;
    Model.d = d ;
    Model.N = pointNum ;
    
    %Model.fileName
    dxName = num2str( d * diax );
    dyName = num2str( d * diay );
    dzName = num2str( d * diaz );
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


%输出|m|kd

mkd = abs(m) * k * d ;

if mkd < 0.5
    fprintf('|m|kd 的值为 %4.2f ，可以适当减小模型增加计算速度 \n' , mkd) ;
elseif mkd > 1
    fprintf('|m|kd 的值为 %4.2f ，请适当增大模型以增加计算精度 \n' , mkd) ;
else
    fprintf('|m|kd 的值为 %4.2f ，模型正常，在没有明显边界效应的情况下此模型良好 \n' , mkd) ;
end

strucFigureHandle = DDA_StrucPlot( Model ) ;
[ memoryRequired , GPUmemoryRequired ] = DDA_MemoryEstimate( size( Model.struc ) ) ;

save(['Model_', Model.fileName ,'.mat'],'Model') ;

end

