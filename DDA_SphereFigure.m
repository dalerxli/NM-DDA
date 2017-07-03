function figHandle = DDA_SphereFigure( tempMuller , rowIter , columnIter )

%DDA_SPHEREFIGURE
%此函数为作图函数 返回的为图像句柄
%作出来的图为散射光的muller矩阵元在空间（球面）分布

%输入参量
%输入tempMuller 一个 N X 18大小的muller球面谱
%输入rowIter columnIter muller矩阵的行、列坐标

%输出
%会输出一幅图像，函数返回值为此图像的句柄

%设置参数
colorNum = 64 ;
az = 139.8400 ;
el = 34.3200 ;
cmap = colormap( jet ( colorNum )) ;

%函数主体部分
figHandle = figure ;
hold on;
xlabel('x') ;
ylabel('y') ;
zlabel('z') ;
view(az, el) ;
title('(激光方向沿y轴正向)') ;

%生成空间点的x y z坐标
angleArray = tempMuller(:,1:2) ;
angleArray = deg2rad( angleArray ) ;
x = sin(angleArray(:,1)) .* cos(angleArray(:,2)) ;
y = cos(angleArray(:,1)) ;
z = -sin(angleArray(:,1)) .* sin(angleArray(:,2)) ;

%将此矩阵元由大到小排列，最后由深到浅上色
%去掉重复的值
inputMullerEleArray = tempMuller(:,rowIter * 4 + columnIter - 2) ;
uniMuller =  unique( inputMullerEleArray ) ;
%每种颜色用于表示的数据点数
eleNumForEveryColor = ceil( length(uniMuller) / colorNum ) ;
%排序
[ sortedUni ,indiceMuller ]  = sort( uniMuller , 'descend') ;

for iter = 1 : size( angleArray , 1 )
    
    %找到此数据点在已排序数据中的位置
    indiceRank = find( sortedUni == inputMullerEleArray( iter ) ) ;
    %此位置应该对应的颜色
    colorRank = ceil( indiceRank / eleNumForEveryColor ) ;
    tempColorArray = cmap( colorRank , :) ;
    %画出此数据点
    plot3(x(iter),y(iter),z(iter),'.','MarkerSize',40,'Color',tempColorArray) ;
    
end

end

