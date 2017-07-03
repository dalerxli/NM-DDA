function strucFigureHandle = DDA_StrucPlot( Model )
%�������� Model ������Model.struc����ΪͼƬ
%����ֻѡ���Եػ��Ƹ�����ı��棬���������ڲ�

%ͼƬ��ʼ����
strucFigureHandle = figure ;
hold on;
az = 139.8400 ;
el = 34.3200 ;
hold on;
xlabel('x') ;
ylabel('y') ;
zlabel('z') ;
view(az, el) ;

%����ͼƬ�߽��
Nx = size( Model.struc , 1) - 1 ;
Ny = size( Model.struc , 2) - 1 ;
Nz = size( Model.struc , 3) - 1 ;
d = Model.d ;
%�õ��߽��
strEdge = bwperim(Model.struc) ;
totalOriPlotNum = sum( strEdge(:) ) ;
%һ��ֻplot 2000����
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

