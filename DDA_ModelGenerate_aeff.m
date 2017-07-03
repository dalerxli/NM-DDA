function  Model = DDA_ModelGenerate_aeff( shape , diaArray , aeff ,m )
%�˺�����DDA_ModelGenerate��aeff��
%�������ӵĵ�Ч�����İ뾶aeff����������m �Լ����ӵ���״shape �Լ�����������״��Ӧ�ĸ����diaArray
%�������Զ������Ӧ�ľ�����d
%�������� |m|kd < 1ʱ�� �����ᷢ�����棨 Ĭ�ϵĲ���Ϊ lambda = 0.532 ��
%�˺���֧������ ��������
%20170418

%����˵��
% shape = 'cylinder' 'ellipsoid' 'brick'�������������ʹ������������ȼ��ɣ�
% diaArray һ�� 1 X 3 ���� ���ֱ��� x y z���������ڱ�ʾ������״�ĸ���� [40 40 40]��ʾһ��ֱ��Ϊ 40*d ����
% �� ע�⣡����תģ�������ΪY�������򡣱���������һ���������Ϊ1.5������ diaArray ����дΪ[40 60 40]
% �����ȷ����DDA_TargetOrienGener������Ĭ�ϵ������йأ���ѡ�ô��������ȡ��ѡȡ��Ϊ��������ʱ�������� ��
% aeff ������ĵ�Ч�뾶 �� ��λΪum
% m  ����ĸ�������

%���˵��
%Model �ṹ��
%Model.shape ˵��model����״���ַ�����ʽ
%Model.aeff  ˵��model�ĵ�Ч����
%Model.m     ˵��model�ĸ�������
%Model.struc һ����ά�������ڱ궨��Щ����ż����
%Model.N     ˵��model��ż���ӵĸ���
%Model.d     ˵��model��ż����֮��ľ���

%����struc����
%struc�ĳ����ΪNx , Ny , Nz
%struc����������Ľڵ�����Nx + 1 , Ny + 1 , Nz + 1
%���� struc��һ��sizeΪ( Nx + 1 , Ny + 1 , Nz + 1 ) ����һ����
%��ģ�͵ĳ����Ϊ Nx * d ��Ny * d ��Nz * d

%��Ҫ��Ϊ�˼ӿ�fft���㣬ע�����ɵ�struc�����䳤��ߵ�2����1 Ӧ����һ���ܸ߶ȷֽ������ӵ���
%���磬Nx = 40    2 * Nx + 1 = 81   81 = 3^4
%�෴ ��� Nx = 41     2 * Nx + 1 = 83   83 ���������Ἣ������fft�ļ�����

%����Ϊ����Narray �˺������Narray����ѡ���԰���ģ�͵���СNx Ny Nz�����ģ��
twoNP1Array = [ 27 45 63 75 81 99 105 117 125 135 147 165 175 189 195 225 231 243 245 ...
    273 275 297 315 325 343 351 363 375 385 405 429 441 455 495 507 525 539 567 585 605 ...
    625 637 675 693 715 729 735 819 825 845];
Narray = ( twoNP1Array - 1 ) / 2 ;

%Ĭ��lambda ���ڼ��� |m|kd < 1
lambda = 0.532 ;
k = 2 * pi / lambda ;

%���²��ֿ�ʼ����Model
Model.m = m ;
Model.shape = shape ;
Model.aeff = aeff ;

diax = diaArray(1) ;
diay = diaArray(2) ;
diaz = diaArray(3) ;

Nx = Narray(find(Narray > diax + 2 , 1)) ;
Ny = Narray(find(Narray > diay + 2 , 1)) ;
Nz = Narray(find(Narray > diaz + 2 , 1)) ;

%������ Model.struc ������ָ����size
Model.struc = zeros( Nx + 1 , Ny + 1 , Nz + 1 ,'logical') ;

%����������
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
    
    %����ż���Ӹ���
    pointNum = sum( Model.struc(:) ) ;
    
    %����ż���Ӽ��d
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

%Բ������� ע���ʱ��Բ����һ�����޵�Բ�� ���������޳�Բ��
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
    
    %����ż���Ӹ���
    pointNum = sum( Model.struc(:) ) ;
    
    %����ż���Ӽ��d
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

%��������
if strcmp( shape , 'brick')
    Model.struc(1:diax , 1: diay , 1:diaz) = 1;
    
    %����ż���Ӹ���
    pointNum = sum( Model.struc(:) ) ;
    
    %����ż���Ӽ��d
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


%���|m|kd

mkd = abs(m) * k * d ;

if mkd < 0.5
    fprintf('|m|kd ��ֵΪ %4.2f �������ʵ���Сģ�����Ӽ����ٶ� \n' , mkd) ;
elseif mkd > 1
    fprintf('|m|kd ��ֵΪ %4.2f �����ʵ�����ģ�������Ӽ��㾫�� \n' , mkd) ;
else
    fprintf('|m|kd ��ֵΪ %4.2f ��ģ����������û�����Ա߽�ЧӦ������´�ģ������ \n' , mkd) ;
end

strucFigureHandle = DDA_StrucPlot( Model ) ;
[ memoryRequired , GPUmemoryRequired ] = DDA_MemoryEstimate( size( Model.struc ) ) ;

save(['Model_', Model.fileName ,'.mat'],'Model') ;

end

