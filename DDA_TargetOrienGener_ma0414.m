function  scaData  = DDA_TargetOrienGener_ma0414( initialLog )
%DDA_TARGETORIENGENER
% �˺�����������ʦ0414���һ�������������ض�ɢ��ǵ���������뷨��
% ����scaData�ṹ��

%���� scaData �ṹ��
%-------.t_thetaMaxOrder
%-------.t_betaMaxOrder
% ------.orieNum
% ------.t_thetaOrder
% ------------------.phiArray
% ------------------.phiMaxNum
% ------------------.t_thetaValue
% ------------------.t_betaOrder
% -----------------------------.t_betaValue
% -----------------------------.horiP
%                              .vertiP
% ����theta phi beta �ĵ�λ��Ϊdeg

%psiΪ����xyƽ��֮�н� ����L_theta �� psi
L_theta = 30 ;
psiInter = 4;

%�õ�t_phi �� t_theta
L_thetaRad = deg2rad( L_theta ) ;
psiArray = [ 0 : psiInter : 180 ] ;
psiArrayRad = deg2rad( psiArray ) ;

%x y zΪ���弸��ȡ����㷽��ѡȡ���y�����ϵ�һ�����յ�
z = sin( psiArrayRad ) ;
x = cos( psiArrayRad ) * sin( L_thetaRad / 2 ) ;
y = cos( psiArrayRad ) * cos( L_thetaRad / 2 ) ;

t_thetaArrayRad = pi - acos( y ) ;
t_phiArrayRad = atan( x ./ z ) ;

% ��theta �� phi ��rad��Ϊdeg
t_thetaArray = rad2deg( t_thetaArrayRad ) ;
t_phiArray = rad2deg( t_phiArrayRad ) ;

scaData.t_thetaMaxOrder = length( t_thetaArray ) ;
scaData.t_betaMaxOrder = 1 ;

%ע��phiMaxNumΪ2 �� ��һ����phiΪ0����� resultOutput�ű��У� ���м���phi��0��������Ǵ�phiΪ0������·�ת����

for iter = 1 : scaData.t_thetaMaxOrder
    
    scaData.t_thetaOrder( iter ).t_thetaValue = t_thetaArray( iter ) ;
    scaData.t_thetaOrder( iter ).t_betaOrder(1).t_betaValue = 0 ;
    scaData.t_thetaOrder( iter ).phiMaxNum = 2 ;
    scaData.t_thetaOrder( iter ).phiArray = [ 0 t_phiArray( iter ) ] ;
    
end
    
    
