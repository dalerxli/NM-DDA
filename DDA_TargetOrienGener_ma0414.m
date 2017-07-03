function  scaData  = DDA_TargetOrienGener_ma0414( initialLog )
%DDA_TARGETORIENGENER
% 此函数用于马老师0414提的一个“柱加球在特定散射角等于椭球的想法”
% 生成scaData结构体

%生成 scaData 结构体
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
% 以上theta phi beta 的单位均为deg

%psi为轴与xy平面之夹角 设置L_theta 和 psi
L_theta = 30 ;
psiInter = 4;

%得到t_phi 和 t_theta
L_thetaRad = deg2rad( L_theta ) ;
psiArray = [ 0 : psiInter : 180 ] ;
psiArrayRad = deg2rad( psiArray ) ;

%x y z为物体几何取向计算方便选取其的y轴向上的一个参照点
z = sin( psiArrayRad ) ;
x = cos( psiArrayRad ) * sin( L_thetaRad / 2 ) ;
y = cos( psiArrayRad ) * cos( L_thetaRad / 2 ) ;

t_thetaArrayRad = pi - acos( y ) ;
t_phiArrayRad = atan( x ./ z ) ;

% 将theta 和 phi 由rad换为deg
t_thetaArray = rad2deg( t_thetaArrayRad ) ;
t_phiArray = rad2deg( t_phiArrayRad ) ;

scaData.t_thetaMaxOrder = length( t_thetaArray ) ;
scaData.t_betaMaxOrder = 1 ;

%注意phiMaxNum为2 ， 第一个是phi为0的情况 resultOutput脚本中， 所有计算phi非0的情况都是从phi为0的情况下翻转而来

for iter = 1 : scaData.t_thetaMaxOrder
    
    scaData.t_thetaOrder( iter ).t_thetaValue = t_thetaArray( iter ) ;
    scaData.t_thetaOrder( iter ).t_betaOrder(1).t_betaValue = 0 ;
    scaData.t_thetaOrder( iter ).phiMaxNum = 2 ;
    scaData.t_thetaOrder( iter ).phiArray = [ 0 t_phiArray( iter ) ] ;
    
end
    
    
