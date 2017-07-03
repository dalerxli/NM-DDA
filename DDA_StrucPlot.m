function strucFigureHandle = DDA_StrucPlot( Model )
%函数接收 Model 并将其Model.struc绘制为图片
%函数只选择性地绘制该物体的表面，不绘制其内部

%图片初始部分
strucFigureHandle = figure ;
hold on;
az = 139.8400 ;
el = 34.3200 ;
hold on;
xlabel('x') ;
ylabel('y') ;
zlabel('z') ;
view(az, el) ;

%绘制图片边界点
Nx = size( Model.struc , 1) - 1 ;
Ny = size( Model.struc , 2) - 1 ;
Nz = size( Model.struc , 3) - 1 ;
d = Model.d ;
%得到边界点
strEdge = bwperim(Model.struc) ;
totalOriPlotNum = sum( strEdge(:) ) ;
%一般只plot 2000个点
ratio = ceil( (totalOriPlotNum / 2000)) ;

for i = 1 : Nx + 1
    for  j = 1 : Ny + 1
        for  k = 1 : Nz + 1
            if strEdge(i,j,k) == 1 && rand(1) < 1 / ratio          
                plot3(i*d ,j*d ,k*d ,'.') ;
            end
        end
    end
end

