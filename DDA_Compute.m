function [ scaData , residuOriArray ] = DDA_Compute( Model,initialLog )
%DDA_COMPUTE

%初始化部分

% 从Model中初始化参数
structure = Model.struc ;
m = Model.m ;
d = Model.d ;

%lamda 为入射光在真空中的波长 k为波矢
lambda = initialLog.lambda ;
k =  2 * pi / lambda ;

%得到Nx Ny Nz alpha
Nx = size(structure,1) - 1 ;
Ny = size(structure,2) - 1 ;
Nz = size(structure,3) - 1 ;

%alpha alphaw alphaLDR 根据Draine的论文 采取alphaLDR最为精确
kd = k * d ;
er = m^2 ;
n = 1 / d^3 ;
alpha0 = 3 / ( 4 * pi * n ) * ( er - 1 ) / ( er + 2 ) ;
b1 = -1.8915316 ;
b2 = 0.1648469 ;
b3 = -1.7700004 ;
S = 0.2 ;
alphanr = alpha0 / ( 1 + ( alpha0 / d^3 ) * ( b1 + m^2 * b2 + m^2 * b3 * S ) * kd^2 ) ;
alphaw = alpha0 / ( 1 - 2 / 3 * 1i * k^3 * alpha0 ) ;
alphaLDR = alphanr / ( 1 - 2 / 3 * 1i * k^3 * alphanr ) ;

%计算部分
%所有变量用single类型
er = single(er) ;
lambda = single(lambda) ;
k = single(k) ;
d = single(d) ;
kd = single(kd) ;
alphaw = single( alphaw ) ;
alphaLDR = single( alphaLDR ) ;


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

switch initialLog.targetOrienGener
    case 'default'
        scaData = DDA_TargetOrienGener( initialLog ) ;
    case 'ma0414'
        scaData = DDA_TargetOrienGener_ma0414( initialLog ) ;
end
%将Model和initialLog的信息存进scaData
scaData.Model = Model ;
scaData.initialLog = initialLog ;

t_thetaMaxOrder = scaData.t_thetaMaxOrder ;
t_betaMaxOrder = scaData.t_betaMaxOrder ;

%residuOriArray用于记载residuOri序列，以便debug
residuOriArray = [];

%定时模式代码
if initialLog.timeLimitFlag == 1
    timeLimitPerLoop = initialLog.desiredTimeCost / t_thetaMaxOrder / t_betaMaxOrder * 0.4 ;
    fprintf('timeLimitPerLoop is %5.4f \n',timeLimitPerLoop) ;
end

for j = 1 : t_thetaMaxOrder
    
    for betaIter = 1 : t_betaMaxOrder
        
        %获取此时的t_theta值，t_beta的值 ， 并给出其弧度值，正弦余弦值
        t_thetaValue = scaData.t_thetaOrder( j ).t_thetaValue ;
        t_betaValue = scaData.t_thetaOrder( j ).t_betaOrder( betaIter ).t_betaValue ;
        t_thetaValueRad = deg2rad( t_thetaValue ) ;
        t_betaValueRad = deg2rad( t_betaValue ) ;
        cosTheta = cos( t_thetaValueRad ) ;
        sinTheta = sin( t_thetaValueRad ) ;
        cosBeta = cos( t_betaValueRad ) ;
        sinBeta = sin( t_betaValueRad ) ;
        
        %生成A部分
        %初始化A 注意Axz 和 Azx完全相等 所以只用生成一个
        Axx = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        Ayy = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        Azz = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        
        Axy = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        Axz = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        Ayz = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        
        fprintf('正在生成A ...... \n') ;
        %三重循环 生成A的每个元
        
        %c1 c1和nx ny nz 无关 避免重复计算
        c1 = kd^2 ;
        
        for nx = -Nx : Nx
            for ny = -Ny : Ny
                for nz = -Nz : Nz
                    
                    %生成一些计算中会用的变量 rx ry rz
                    rx = nx * cosBeta + nz * sinBeta ;
                    ry = ny * cosTheta + ( nz * cosBeta - nx * sinBeta) * sinTheta ;
                    rz = ( nz * cosBeta - nx * sinBeta ) * cosTheta - ny * sinTheta ;
                    r = ( rx^2 + ry^2 + rz^2)^0.5  ;
                    
                    %生成A
                    if nx ~= 0 || ny ~= 0 || nz ~= 0
                        c3 = exp( 1i * kd * r) / d^3 / r^3 ;
                        %c2 c31
                        c2 = (1 - kd * 1i * r) / r^2;
                        
                        Axx(nx + Nx + 1, ny + Ny + 1 ,nz + Nz + 1) = c3 * (c1 * ( -ry^2 - rz^2 ) + c2 * ( ry^2 + rz^2 - 2*rx^2 )) ;
                        Ayy(nx + Nx + 1, ny + Ny + 1 ,nz + Nz + 1) = c3 * (c1 * ( -rx^2 - rz^2 ) + c2 * ( rx^2 + rz^2 - 2*ry^2 )) ;
                        Azz(nx + Nx + 1, ny + Ny + 1 ,nz + Nz + 1) = c3 * (c1 * ( -rx^2 - ry^2 ) + c2 * ( rx^2 + ry^2 - 2*rz^2)) ;
                        
                        Axy(nx + Nx + 1, ny + Ny + 1 ,nz + Nz + 1) = c3 * (c1 * ( ry*rx ) + c2 * ( -3 * rx * ry )) ;
                        Axz(nx + Nx + 1, ny + Ny + 1 ,nz + Nz + 1) = c3 * (c1 * ( rz*rx ) + c2 * ( -3 * rx * rz )) ;
                        Ayz(nx + Nx + 1, ny + Ny + 1 ,nz + Nz + 1) = c3 * (c1 * ( rz*ry ) + c2 * ( -3 * rz * ry )) ;
                        
                    end
                end
            end
        end
        fprintf('A生成完毕 ......\n') ;
        
        %运用卷积定理，使用fft加速计算
        %对A做fft 并将其储存入GPU之中
        Axxfg = gpuArray( fftn(Axx) ) ;
        Ayyfg = gpuArray( fftn(Ayy) ) ;
        Azzfg = gpuArray( fftn(Azz) ) ;
        
        Axyfg = gpuArray( fftn(Axy) ) ;
        Axzfg = gpuArray( fftn(Axz) ) ;
        Ayzfg = gpuArray( fftn(Ayz) ) ;
        
        %释放A所占的内存
        clear Axx Axy Axz Ayy Ayz Azz
        
        for eDirection = 0 : 1
            
            %定时模式代码
            if initialLog.timeLimitFlag == 1
                tempStartTime = clock ;
            end
            
            %初始化P
            Px = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            Py = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            Pz = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            %初始化E
            Ex = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            Ey = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            Ez = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            
            %三重循环 生成E 视eDirection的不同而生成不同的电场
            fprintf('正在生成电场E ...... \n') ;
            if eDirection == 0
                for nx = 0 : Nx
                    for ny = 0 : Ny
                        for nz = 0 : Nz
                            Ex( nx + 1 , ny + 1 , nz + 1 ) = exp( 1i * kd * ( ny * cosTheta + ...
                                ( nz * cosBeta - nx * sinBeta) * sinTheta)) ;
                        end
                    end
                end
                
            else
                for nx = 0 : Nx
                    for ny = 0 : Ny
                        for nz = 0 : Nz
                            Ez( nx + 1 , ny + 1 , nz + 1 ) = exp( 1i * kd * ( ny * cosTheta + ...
                                ( nz * cosBeta - nx * sinBeta) * sinTheta)) ;
                        end
                    end
                end
                
            end
            
            %去无效点
            Ex = Ex .* structure ;
            Ey = Ey .* structure ;
            Ez = Ez .* structure ;
            
            
            %迭代计算
            %迭代中需要的变量的初始化
            iterCount = 0;
            
            Zx = DDA_ConvAccelerate( Axxfg , Ex , 1) + DDA_ConvAccelerate( Axyfg , Ey , 1) + DDA_ConvAccelerate( Axzfg , Ez , 1) + Ex / conj(alphaLDR) ;
            Zy = DDA_ConvAccelerate( Axyfg , Ex , 1) + DDA_ConvAccelerate( Ayyfg , Ey , 1) + DDA_ConvAccelerate( Ayzfg , Ez , 1) + Ey / conj(alphaLDR) ;
            Zz = DDA_ConvAccelerate( Axzfg , Ex , 1) + DDA_ConvAccelerate( Ayzfg , Ey , 1) + DDA_ConvAccelerate( Azzfg , Ez , 1) + Ez / conj(alphaLDR) ;
            
            %去无效点
            Zx = Zx .* structure;
            Zy = Zy .* structure;
            Zx = Zx .* structure;
            
            gx = Zx;
            gy = Zy;
            gz = Zz;
            
            px = gx;
            py = gy;
            pz = gz;
            
            wx = zeros(size(Px),'single') ;
            wy = zeros(size(Py),'single') ;
            wz = zeros(size(Pz),'single') ;
            
            vx = DDA_ConvAccelerate( Axxfg , px , 0) + DDA_ConvAccelerate( Axyfg , py , 0) + DDA_ConvAccelerate( Axzfg , pz , 0) + px / alphaLDR ;
            vy = DDA_ConvAccelerate( Axyfg , px , 0) + DDA_ConvAccelerate( Ayyfg , py , 0) + DDA_ConvAccelerate( Ayzfg , pz , 0) + py / alphaLDR ;
            vz = DDA_ConvAccelerate( Axzfg , px , 0) + DDA_ConvAccelerate( Ayzfg , py , 0) + DDA_ConvAccelerate( Azzfg , pz , 0) + pz / alphaLDR ;
            
            %去无效点
            vx = vx .* structure ;
            vy = vy .* structure ;
            vz = vz .* structure ;
            
            EiEi0 = Ex(:)' * Ex(:) + Ey(:)' * Ey(:) + Ez(:)' *Ez(:) ;
            
            residuOri = 1;
            %数值初始化完毕，进入迭代
            fprintf('数值初始化完毕，进入迭代 \n');
            
            %记录 residuOri
            
            
            while  residuOri > initialLog.residu
                
                %定时模式代码
                if initialLog.timeLimitFlag == 1
                    tempCurrentTime = clock ;
                end
                
                %判断是否超时，导致跳出循环
                if TimePast( tempStartTime , tempCurrentTime ) > timeLimitPerLoop && residuOri < initialLog.timeModeLeastResidu
                    fprintf('超时，此时residu为 %5.4f \n',residuOri ) ;
                    break ;
                end
                
                iterCount = iterCount + 1 ;
                
                %ai 原文为alpha，与极化率冲突，故选用ai。 ai的意义是步长的系数。gigi 记为gi的内积和
                gigi = gx(:)' * gx(:) + gy(:)' * gy(:) + gz(:)' *gz(:) ;
                ai = gigi / (vx(:)' * vx(:) + vy(:)' * vy(:) + vz(:)' *vz(:) );
                
                %没有点的地方，P为0
                Px = Px + ai * px ;
                Py = Py + ai * py ;
                Pz = Pz + ai * pz ;
                
                wx = wx + ai * vx ;
                wy = wy + ai * vy ;
                wz = wz + ai * vz ;
                
                gx = Zx - DDA_ConvAccelerate( Axxfg , wx , 1) - DDA_ConvAccelerate( Axyfg , wy , 1) - DDA_ConvAccelerate( Axzfg , wz , 1) - wx / conj(alphaLDR) ;
                gy = Zy - DDA_ConvAccelerate( Axyfg , wx , 1) - DDA_ConvAccelerate( Ayyfg , wy , 1) - DDA_ConvAccelerate( Ayzfg , wz , 1) - wy / conj(alphaLDR) ;
                gz = Zz - DDA_ConvAccelerate( Axzfg , wx , 1) - DDA_ConvAccelerate( Ayzfg , wy , 1) - DDA_ConvAccelerate( Azzfg , wz , 1) - wz / conj(alphaLDR) ;
                
                %去无效点
                gx = gx .* structure ;
                gy = gy .* structure ;
                gz = gz .* structure ;
                
                beta = (gx(:)' * gx(:) + gy(:)' * gy(:) + gz(:)' *gz(:)) / gigi ;
                
                px = gx + beta * px ;
                py = gy + beta * py ;
                pz = gz + beta * pz ;
                
                vx = DDA_ConvAccelerate( Axxfg , gx , 0) + DDA_ConvAccelerate( Axyfg , gy , 0) + DDA_ConvAccelerate( Axzfg , gz , 0) + gx / alphaLDR + beta * vx ;
                vy = DDA_ConvAccelerate( Axyfg , gx , 0) + DDA_ConvAccelerate( Ayyfg , gy , 0) + DDA_ConvAccelerate( Ayzfg , gz , 0) + gy / alphaLDR + beta * vy ;
                vz = DDA_ConvAccelerate( Axzfg , gx , 0) + DDA_ConvAccelerate( Ayzfg , gy , 0) + DDA_ConvAccelerate( Azzfg , gz , 0) + gz / alphaLDR + beta * vz ;
                
                %去无效点
                vx = vx .* structure ;
                vy = vy .* structure ;
                vz = vz .* structure ;
                
                %每5次迭代更新一次residuOri
                %residuOri
                if mod( iterCount , 5 ) == 1
                    
                    residuOri = ((Ex(:) - wx(:))'*(Ex(:) - wx(:)) + (Ey(:) - wy(:))'*(Ey(:) - wy(:)) + (Ez(:) - wz(:))'*(Ez(:) - wz(:))) / EiEi0 ;
                    
                    fprintf('迭代次数 %i 当前原方程residu为 %d \n', iterCount , residuOri);
                    
                end
            end
            %记录此residuOri
            residuOriArray = [residuOriArray , residuOri] ;
            fprintf('完成了一次迭代计算。当前t_theta为 %d 序数为 %d ，当前电场方向为 %d ，最大t_theta序数为 %d \n',t_thetaValueRad,j,eDirection,t_thetaMaxOrder) ;
            
            %生成P，P的三分量为Px Py Pz
            P = zeros(Nx + 1, Ny + 1, Nz + 1,3,'single') ;
            P(:, :, : ,1) = Px ;
            P(:, :, : ,2) = Py ;
            P(:, :, : ,3) = Pz ;
            
            %将P的值存入scaData
            if eDirection == 0
                scaData.t_thetaOrder(j).t_betaOrder( betaIter ).horiP = P ;
                clear P ;
            end
            
            if eDirection == 1
                scaData.t_thetaOrder(j).t_betaOrder( betaIter ).vertiP = P ;
                clear P ;
            end
            
        end
    end
end

end


