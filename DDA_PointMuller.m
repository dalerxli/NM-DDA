function pointMuller = DDA_PointMuller( tempMuller , L_theta , L_phi )
%DDA_POINTMULLER
%注：此函数含有静态变量，用查表的方式加快速度
%此函数实现从球面谱muller中找到空间某一点的muller矩阵
%此函数解决的核心问题是：生成tempMuller球面谱是空间中均匀分布的离散点
%如果有希望得到的点在原数据点列中没有，那么必然需要周围的点进行拟合得到此点数据

%输入部分
%函数接受一个 tempMuller（即球面谱的muller）
%L_theta L_phi为实验室坐标系下你想知道的点的极坐标 角度制

%输出部分
%pointMuller 一个4 X 4 muller矩阵

% tempCoe 为函数的静态变量， 用于储存为了得到某个特定实验室探测角而需要对原tempMuller操作的系数
% 比如对2586个点、  L_theta = 30 L_phi = 0 这样的情况， tempCoe会有mark = [2586 30 0]这样的标签
% 并在对应的coePosArray中给出系数值
% coePosArray的构造如下： 
% [ 第一个tempmuller的行数 第一个权重 ; 
%   ...         ...
%    第N个行数             第N个权重 ... ]
% coePosArray 最多4行
% coePosArray 的所有权重和为1
% coePosArray 的生成由 DDA_PointMullerCoeGener( tempMuller , L_theta , L_phi ) 实现

% 函数主体部分
% 查表，查找以前是否计算过对应角度的系数
% 运行完这段代码后一定可以生成 temp_coePosArray 

persistent tempCoe ;
if isempty( tempCoe )
    tempCoe(1).mark = [ size( tempMuller , 1) , L_theta , L_phi ] ;
    tempCoe(1).coePosArray = DDA_PointMullerCoeGener( tempMuller , L_theta , L_phi ) ;
    temp_coePosArray = tempCoe(1).coePosArray ;
else
    temp_coePosArray = [];
    for iter = 1 : size( tempCoe , 2 )
        if tempCoe( iter ).mark(1) == size( tempMuller , 1) && ...
                tempCoe( iter ).mark(2) == L_theta && ...
                tempCoe( iter ).mark(3) == L_phi 
            temp_coePosArray = tempCoe( iter ).coePosArray ;
        end
    end
    
    if isempty( temp_coePosArray )
        newIndice =  size( tempCoe , 2 ) + 1 ;
        tempCoe( newIndice ).mark = [ size( tempMuller , 1) , L_theta , L_phi ] ;
        tempCoe( newIndice ).coePosArray = DDA_PointMullerCoeGener( tempMuller , L_theta , L_phi ) ;
        temp_coePosArray = tempCoe( newIndice ).coePosArray ;
    end
end

%将周围点的数据叠加 得到此点数据
tempPointMullerRow = zeros(1,18) ;
for iter = 1 : size( temp_coePosArray , 1 )
    tempPointMullerRow = tempPointMullerRow + temp_coePosArray( iter , 2 ) * tempMuller( temp_coePosArray( iter , 1 ) , : ) ;
end

%将值赋值给pointMuller
pointMuller = zeros(4,4) ;
for iterRow = 1 : 4
    for iterColumn = 1 : 4
        pointMuller( iterRow , iterColumn ) = tempPointMullerRow( 4 * iterRow + iterColumn - 2 ) ;
    end
end

end



