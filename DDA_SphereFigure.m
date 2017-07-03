function figHandle = DDA_SphereFigure( tempMuller , rowIter , columnIter )

%DDA_SPHEREFIGURE
%�˺���Ϊ��ͼ���� ���ص�Ϊͼ����
%��������ͼΪɢ����muller����Ԫ�ڿռ䣨���棩�ֲ�

%�������
%����tempMuller һ�� N X 18��С��muller������
%����rowIter columnIter muller������С�������

%���
%�����һ��ͼ�񣬺�������ֵΪ��ͼ��ľ��

%���ò���
colorNum = 64 ;
az = 139.8400 ;
el = 34.3200 ;
cmap = colormap( jet ( colorNum )) ;

%�������岿��
figHandle = figure ;
hold on;
xlabel('x') ;
ylabel('y') ;
zlabel('z') ;
view(az, el) ;
title('(���ⷽ����y������)') ;

%���ɿռ���x y z����
angleArray = tempMuller(:,1:2) ;
angleArray = deg2rad( angleArray ) ;
x = sin(angleArray(:,1)) .* cos(angleArray(:,2)) ;
y = cos(angleArray(:,1)) ;
z = -sin(angleArray(:,1)) .* sin(angleArray(:,2)) ;

%���˾���Ԫ�ɴ�С���У�������ǳ��ɫ
%ȥ���ظ���ֵ
inputMullerEleArray = tempMuller(:,rowIter * 4 + columnIter - 2) ;
uniMuller =  unique( inputMullerEleArray ) ;
%ÿ����ɫ���ڱ�ʾ�����ݵ���
eleNumForEveryColor = ceil( length(uniMuller) / colorNum ) ;
%����
[ sortedUni ,indiceMuller ]  = sort( uniMuller , 'descend') ;

for iter = 1 : size( angleArray , 1 )
    
    %�ҵ������ݵ��������������е�λ��
    indiceRank = find( sortedUni == inputMullerEleArray( iter ) ) ;
    %��λ��Ӧ�ö�Ӧ����ɫ
    colorRank = ceil( indiceRank / eleNumForEveryColor ) ;
    tempColorArray = cmap( colorRank , :) ;
    %���������ݵ�
    plot3(x(iter),y(iter),z(iter),'.','MarkerSize',40,'Color',tempColorArray) ;
    
end

end

