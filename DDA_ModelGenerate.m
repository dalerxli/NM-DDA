function  Model  = DDA_ModelGenerate( shape , diaArray , m , accuracy)
%��������һ��Model
%Model�ṹ����ȫ��������һ�����ӵ���״����С�������ʡ��Լ���Ӧ�Ŀռ����

%����˵��
% shape = 'cylinder' 'ellipsoid' 'brick'�������������ʹ������������ȼ��ɣ�
% diaArray ���ӵ�X Y Z�����ϵ��߶ȣ���λΪum ����, 'ellipsoid' [ 1 2 1] ����һ��ֱ���ֱ�Ϊ 1 2 1
% um������ �� 'cylinder' [1 2 1] ����һ���� 2um ��ֱ�� 1 um ��Բ��
% m ���ӵĸ�������
% ע�⣡����תģ�������ΪY��������
% �����ȷ����DDA_TargetOrienGener������Ĭ�ϵ������йأ���ѡ�ô��������ȡ��ѡȡ��Ϊ��������ʱ��������
% |m|kd��ֵ��Ӱ����㾫�� ���ϣ���߾��ȼ��㣬����ͨ������ accuracy ���
% accuracy = 'Normal' 'High' 'Ultra' ��Ӧ��|m|kd �ֱ�Ϊ 0.8 0.5 0.25 ��ȱʡֵΪ 'Normal'

%���˵��
%Model �ṹ��
%Model.shape ˵��model����״���ַ�����ʽ
%Model.aeff  ˵��model�ĵ�Ч����
%Model.m     ˵��model�ĸ�������
%Model.struc һ����ά�������ڱ궨��Щ����ż����
%Model.N     ˵��model��ż���ӵĸ���
%Model.d     ˵��model��ż����֮��ľ���
%Model.fileName Ϊ�˷��㱣������ɵ�����

%����� DDA_ModelGenerate_aeff ���˺����������û������������ӿռ��߶ȣ�������Ӧѡ��d�����ɿռ����

%����struc����
%struc�ĳ����ΪNx , Ny , Nz
%struc����������Ľڵ�����Nx + 1 , Ny + 1 , Nz + 1
%���� struc��һ��sizeΪ( Nx + 1 , Ny + 1 , Nz + 1 ) ����һ����
%��ģ�͵ĳ����Ϊ Nx * d ��Ny * d ��Nz * d

%����Ϊ����Narray �˺������Narray����ѡ���԰���ģ�͵���СNx Ny Nz�����ģ��
twoNP1Array = [ 27 45 63 75 81 99 105 117 125 135 147 165 175 189 195 225 231 243 245 ...
    273 275 297 315 325 343 351 363 375 385 405 429 441 455 495 507 525 539 567 585 605 ...
    625 637 675 693 715 729 735 819 825 845];
Narray = ( twoNP1Array - 1 ) / 2 ;

% Ĭ��lambda
lambda = 0.532 ;
k = 2 * pi / lambda ;

% ���²��ֿ�ʼ����Model
Model.m = m ;
Model.shape = shape ;

if nargin < 4
    accuracy = 'Normal' ;
end

% �õ�d
switch accuracy
    case 'Normal'
        d = 0.8 / k / abs( m ) ;
    case 'High'
        d = 0.5 / k / abs( m ) ;
    case 'Ultra'
        d = 0.25 / k / abs( m ) ;
end

%��ʼ���� ...
diax = diaArray(1) ;
diay = diaArray(2) ;
diaz = diaArray(3) ;

%��������
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
    
    %����ż���Ӹ���
    pointNum = sum( Model.struc(:) ) ;
    
    %��������d �ԶԵ������
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

%Բ�������
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
    
    %����ż���Ӹ���
    pointNum = sum( Model.struc(:) ) ;
    
    %��������d �ԶԵ������
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

%��������
if strcmp( shape , 'brick')
    
    nx = ceil( diax / d ) ;
    ny = ceil( diay / d ) ;
    nz = ceil( diaz / d ) ;
    
    Model.struc(1:nx , 1: ny , 1:nz) = 1;
    
    %����ż���Ӹ���
    pointNum = sum( Model.struc(:) ) ;
    
    %��������d �ԶԵ������
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

%aeff��ż������N
Model.N = pointNum ;
Model.aeff = ( pointNum * d^3 / pi / 4 * 3 )^(1/3) ;
Model.d = d ;

strucFigureHandle = DDA_StrucPlot( Model ) ;
[ memoryRequired , GPUmemoryRequired ] = DDA_MemoryEstimate( size( Model.struc ) ) ;

save(['Model_', Model.fileName ,'.mat'],'Model') ;
end



