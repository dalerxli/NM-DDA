%clear all
clear all ;
profile on
%load
load('scaData.mat') ;
load('angleArray.mat') ;

%读取Model和initialLog
Model = scaData.Model ;
initialLog = scaData.initialLog ;

% 从Model中初始化参数
structure = Model.struc;
m = Model.m;
d = Model.d;

%得到Nx Ny Nz alpha
Nx = size(structure,1) - 1;
Ny = size(structure,2) - 1;
Nz = size(structure,3) - 1;

%lamda 为入射光在真空中的波长 k为波矢
lambda = initialLog.lambda ;
k =  2 * pi / lambda ;
kd = k * d ;

%single
k = single(k) ;
d = single(d) ;
kd = single(kd) ;

t_thetaMaxOrder = scaData.t_thetaMaxOrder ;
t_betaMaxOrder = scaData.t_betaMaxOrder ;
maxL_AangleNum = size(angleArray , 1) ;

%t_theta和t_beta的for循环计算在一种姿态下的散射图样，t_phi的变化可以靠变换来得到，会在下一步的计算中给出
for t_thetaIter = 1 : t_thetaMaxOrder
    for t_betaIter = 1 : t_betaMaxOrder
        
        fprintf('正在计算 thetaIter : %d ; betaIter : %d \n',t_thetaIter ,t_betaIter ) ;
        
        %获取此时的horiP和vertiP，将其暂时储存
        this_horiP = scaData.t_thetaOrder(t_thetaIter).t_betaOrder(t_betaIter).horiP ;
        this_vertiP = scaData.t_thetaOrder(t_thetaIter).t_betaOrder(t_betaIter).vertiP ;
        
        %获取此时t_theta和t_beta的值，并给出其弧度制的值，给出其正弦、余弦值减少计算量
        t_thetaValue = scaData.t_thetaOrder(t_thetaIter).t_thetaValue;
        t_betaValue = scaData.t_thetaOrder(t_thetaIter).t_betaOrder(t_betaIter).t_betaValue ;
        t_thetaValueRad = deg2rad( t_thetaValue ) ;
        t_betaValueRad = deg2rad( t_betaValue ) ;
        
        sinTtheta = sin(t_thetaValueRad) ;
        cosTtheta = cos(t_thetaValueRad) ;
        sinTbeta = sin(t_betaValueRad) ;
        cosTbeta = cos(t_betaValueRad) ;
        
        %初始化要存储的muller矩阵 ，前两列为此行对应的角度（角度制）
        tempMuller = zeros(maxL_AangleNum,18) ;
        tempMuller(:,1:2) = angleArray(:,3:4) ;
        
        %初始化对应的Jones矩阵 依次为S1 S2 S3 S4 对应的矩阵元为 a22 a11 a12 a21
        %见Absorption and Scattering of Light by Small Particles 书P73 电子版
        S1 = zeros(maxL_AangleNum,1) ;
        S2 = zeros(maxL_AangleNum,1) ;
        S3 = zeros(maxL_AangleNum,1) ;
        S4 = zeros(maxL_AangleNum,1) ;
        
        %计算部分
        %遍历探测器的方位进行计算
        parfor L_AngleIter = 1 : maxL_AangleNum
            
            %提取要计算的点的方位角（弧度制）
            L_thetaRad = angleArray(L_AngleIter , 1) ;
            L_phiRad = angleArray(L_AngleIter , 2) ;
            
            if mod(L_AngleIter,5) == 0
                fprintf('L_AngleIter now %d\n',L_AngleIter) ;
            end
            
            %出射光的水平、竖直方向，用于提取Jones矩阵中的S1 S2 S3 S4
            horiUni = [cos(L_thetaRad)*cos(L_phiRad) , -sin(L_thetaRad) , -sin(L_phiRad) * cos(L_thetaRad) ] ;
            vertiUni = [sin(L_phiRad),0,cos(L_phiRad)] ;
            
            %出射光的方向 给出其转置，加快计算速度
            rUni = [sin( L_thetaRad ) * cos( L_phiRad ) , cos( L_thetaRad ) , -sin( L_thetaRad ) * sin( L_phiRad ) ] ;
            rUniInver = rUni' ;
            %单位方向的叉乘矩阵，用于加快计算
            rcross = [0 -rUni(3) rUni(2) ; rUni(3) 0 -rUni(1) ; -rUni(2) rUni(1) 0 ] ;
            rcross2 = rcross * rcross ;
            
            %以下部分为循环计算
            horiE = 0 ;
            vertiE = 0 ;
            
            for nx = 0 : Nx
                for ny = 0 : Ny
                    for nz = 0 : Nz
                        
                        %正在计算的偶极子的位置 pos1,2,3分别代表偶极子的x y z位矢分量
                        pos = zeros(1,3) ;
                        pos(1) = nx * cosTbeta + nz * sinTbeta ;
                        pos(2) = -nx * sinTbeta * sinTtheta + ny * cosTtheta + nz * sinTtheta * cosTbeta ;
                        pos(3) = -nx * sinTbeta * cosTtheta - ny * sinTtheta + nz * cosTtheta * cosTbeta ;
                        
                        %获取此偶极子的极化向量
                        %                         tempHoriP(:) = squeeze(this_horiP(nx+1 , ny+1 , nz+1 ,: )) ;
                        %                         tempVeriP(:) = this_vertiP(nx+1 , ny+1 , nz+1 ,: ) ;
                        
                        %计算部分
                        dot_PosrUni = pos * rUniInver ;
                        tempCoeToMulti =  exp(-1i * kd * dot_PosrUni) * rcross2 ;
                        
                        horiE =   horiE - tempCoeToMulti * squeeze(this_horiP(nx+1 , ny+1 , nz+1 ,: ))  ;
                        vertiE = vertiE - tempCoeToMulti * squeeze(this_vertiP(nx+1 , ny+1 , nz+1 ,: ))  ;
                    end
                end
            end
            
            %Absorption and Scattering of Light by Small Particles书P75左右（电子版）
            S2( L_AngleIter ) = horiUni * horiE ;
            S4( L_AngleIter ) = vertiUni * horiE ;
            
            S3( L_AngleIter ) = horiUni * vertiE ;
            S1( L_AngleIter ) = vertiUni * vertiE ;
            
        end
        
        %将Muller矩阵信息存入tempMuller中，其顺序为S11 S12 S13 S14 S21 ...
        tempMuller(:,3) = ( S1.*conj(S1) + S2.*conj(S2) + S3.*conj(S3) + S4.*conj(S4) ) / 2;
        tempMuller(:,4) = ( S2.*conj(S2) - S1.*conj(S1) + S4.*conj(S4) - S3.*conj(S3) ) / 2;
        tempMuller(:,5) = real( S2.*conj(S3) + S1.*conj(S4) ) ;
        tempMuller(:,6) = imag( S2.*conj(S3) - S1.*conj(S4) ) ;
        
        tempMuller(:,7) = ( -S1.*conj(S1) + S2.*conj(S2) + S3.*conj(S3) - S4.*conj(S4) ) / 2;
        tempMuller(:,8) = ( S1.*conj(S1) + S2.*conj(S2) - S3.*conj(S3) - S4.*conj(S4) ) / 2;
        tempMuller(:,9) = real(S2.*conj(S3) - S1.*conj(S4) );
        tempMuller(:,10) = imag(S2.*conj(S3) + S1.*conj(S4) );
        
        tempMuller(:,11) = real(S2.*conj(S4) + S1.*conj(S3) );
        tempMuller(:,12) = real(S2.*conj(S4) - S1.*conj(S3) );
        tempMuller(:,13) = real(S1.*conj(S2) + S3.*conj(S4) );
        tempMuller(:,14) = imag(S2.*conj(S1) + S4.*conj(S3) );
        
        tempMuller(:,15) = imag(S4.*conj(S2) + S1.*conj(S3) );
        tempMuller(:,16) = imag(S4.*conj(S2) - S1.*conj(S3) );
        tempMuller(:,17) = imag(S1.*conj(S2) - S3.*conj(S4) );
        tempMuller(:,18) = real(S1.*conj(S2) - S3.*conj(S4) );
        tempMuller = single(tempMuller) ;
        outPutMuller.t_thetaOrder( t_thetaIter ).t_phiOrder( 1 ).t_betaOrder( t_betaIter ).Muller = tempMuller ;
        
        %补全phi
        %t_phi通过旋转矩阵来实现
        %以下代码将t_phiOrder补齐
        
        phiArray = scaData.t_thetaOrder( t_thetaIter ).phiArray ;
        
        %以下代码将t_phiOrder补齐
        phiMaxNum = scaData.t_thetaOrder( t_thetaIter ).phiMaxNum ;
        if phiMaxNum ~= 1
            for t_phiIter = 2 : phiMaxNum
                %得到需要旋转的phi角
                t_phiDeg = phiArray( t_phiIter ) ;
                outPutMuller.t_thetaOrder( t_thetaIter ).t_phiOrder( t_phiIter ).t_betaOrder( t_betaIter ).Muller ...
                    = DDA_MullerT_PhiFillUp( t_phiDeg , tempMuller ) ;
            end
        end
    end
end
%记录initialLog Model
outPutMuller.initialLog = initialLog ;
outPutMuller.Model = Model ;
save('outPutMuller.mat','outPutMuller');
%作图部分
%生成x y z
x = sin(angleArray(:,1)) .* cos(angleArray(:,2)) ;
y = cos(angleArray(:,1)) ;
z = -sin(angleArray(:,1)) .* sin(angleArray(:,2)) ;
tempMuller = outPutMuller.thetaOrder( 1 ).t_phiOrder( 1 ).t_betaOrder( 1 ).Muller ;

maxNum = max(tempMuller(:,3)) ;

%以下代码将t_phiOrder补齐

