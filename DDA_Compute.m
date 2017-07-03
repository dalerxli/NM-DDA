function [ scaData , residuOriArray ] = DDA_Compute( Model,initialLog )
%DDA_COMPUTE

%��ʼ������

% ��Model�г�ʼ������
structure = Model.struc ;
m = Model.m ;
d = Model.d ;

%lamda Ϊ�����������еĲ��� kΪ��ʸ
lambda = initialLog.lambda ;
k =  2 * pi / lambda ;

%�õ�Nx Ny Nz alpha
Nx = size(structure,1) - 1 ;
Ny = size(structure,2) - 1 ;
Nz = size(structure,3) - 1 ;

%alpha alphaw alphaLDR ����Draine������ ��ȡalphaLDR��Ϊ��ȷ
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

%���㲿��
%���б�����single����
er = single(er) ;
lambda = single(lambda) ;
k = single(k) ;
d = single(d) ;
kd = single(kd) ;
alphaw = single( alphaw ) ;
alphaLDR = single( alphaLDR ) ;


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

switch initialLog.targetOrienGener
    case 'default'
        scaData = DDA_TargetOrienGener( initialLog ) ;
    case 'ma0414'
        scaData = DDA_TargetOrienGener_ma0414( initialLog ) ;
end
%��Model��initialLog����Ϣ���scaData
scaData.Model = Model ;
scaData.initialLog = initialLog ;

t_thetaMaxOrder = scaData.t_thetaMaxOrder ;
t_betaMaxOrder = scaData.t_betaMaxOrder ;

%residuOriArray���ڼ���residuOri���У��Ա�debug
residuOriArray = [];

%��ʱģʽ����
if initialLog.timeLimitFlag == 1
    timeLimitPerLoop = initialLog.desiredTimeCost / t_thetaMaxOrder / t_betaMaxOrder * 0.4 ;
    fprintf('timeLimitPerLoop is %5.4f \n',timeLimitPerLoop) ;
end

for j = 1 : t_thetaMaxOrder
    
    for betaIter = 1 : t_betaMaxOrder
        
        %��ȡ��ʱ��t_thetaֵ��t_beta��ֵ �� �������仡��ֵ����������ֵ
        t_thetaValue = scaData.t_thetaOrder( j ).t_thetaValue ;
        t_betaValue = scaData.t_thetaOrder( j ).t_betaOrder( betaIter ).t_betaValue ;
        t_thetaValueRad = deg2rad( t_thetaValue ) ;
        t_betaValueRad = deg2rad( t_betaValue ) ;
        cosTheta = cos( t_thetaValueRad ) ;
        sinTheta = sin( t_thetaValueRad ) ;
        cosBeta = cos( t_betaValueRad ) ;
        sinBeta = sin( t_betaValueRad ) ;
        
        %����A����
        %��ʼ��A ע��Axz �� Azx��ȫ��� ����ֻ������һ��
        Axx = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        Ayy = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        Azz = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        
        Axy = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        Axz = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        Ayz = zeros(2 * Nx + 1, 2 * Ny + 1, 2 * Nz + 1,'single');
        
        fprintf('��������A ...... \n') ;
        %����ѭ�� ����A��ÿ��Ԫ
        
        %c1 c1��nx ny nz �޹� �����ظ�����
        c1 = kd^2 ;
        
        for nx = -Nx : Nx
            for ny = -Ny : Ny
                for nz = -Nz : Nz
                    
                    %����һЩ�����л��õı��� rx ry rz
                    rx = nx * cosBeta + nz * sinBeta ;
                    ry = ny * cosTheta + ( nz * cosBeta - nx * sinBeta) * sinTheta ;
                    rz = ( nz * cosBeta - nx * sinBeta ) * cosTheta - ny * sinTheta ;
                    r = ( rx^2 + ry^2 + rz^2)^0.5  ;
                    
                    %����A
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
        fprintf('A������� ......\n') ;
        
        %���þ������ʹ��fft���ټ���
        %��A��fft �����䴢����GPU֮��
        Axxfg = gpuArray( fftn(Axx) ) ;
        Ayyfg = gpuArray( fftn(Ayy) ) ;
        Azzfg = gpuArray( fftn(Azz) ) ;
        
        Axyfg = gpuArray( fftn(Axy) ) ;
        Axzfg = gpuArray( fftn(Axz) ) ;
        Ayzfg = gpuArray( fftn(Ayz) ) ;
        
        %�ͷ�A��ռ���ڴ�
        clear Axx Axy Axz Ayy Ayz Azz
        
        for eDirection = 0 : 1
            
            %��ʱģʽ����
            if initialLog.timeLimitFlag == 1
                tempStartTime = clock ;
            end
            
            %��ʼ��P
            Px = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            Py = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            Pz = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            %��ʼ��E
            Ex = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            Ey = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            Ez = zeros(Nx + 1, Ny + 1, Nz + 1,'single');
            
            %����ѭ�� ����E ��eDirection�Ĳ�ͬ�����ɲ�ͬ�ĵ糡
            fprintf('�������ɵ糡E ...... \n') ;
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
            
            %ȥ��Ч��
            Ex = Ex .* structure ;
            Ey = Ey .* structure ;
            Ez = Ez .* structure ;
            
            
            %��������
            %��������Ҫ�ı����ĳ�ʼ��
            iterCount = 0;
            
            Zx = DDA_ConvAccelerate( Axxfg , Ex , 1) + DDA_ConvAccelerate( Axyfg , Ey , 1) + DDA_ConvAccelerate( Axzfg , Ez , 1) + Ex / conj(alphaLDR) ;
            Zy = DDA_ConvAccelerate( Axyfg , Ex , 1) + DDA_ConvAccelerate( Ayyfg , Ey , 1) + DDA_ConvAccelerate( Ayzfg , Ez , 1) + Ey / conj(alphaLDR) ;
            Zz = DDA_ConvAccelerate( Axzfg , Ex , 1) + DDA_ConvAccelerate( Ayzfg , Ey , 1) + DDA_ConvAccelerate( Azzfg , Ez , 1) + Ez / conj(alphaLDR) ;
            
            %ȥ��Ч��
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
            
            %ȥ��Ч��
            vx = vx .* structure ;
            vy = vy .* structure ;
            vz = vz .* structure ;
            
            EiEi0 = Ex(:)' * Ex(:) + Ey(:)' * Ey(:) + Ez(:)' *Ez(:) ;
            
            residuOri = 1;
            %��ֵ��ʼ����ϣ��������
            fprintf('��ֵ��ʼ����ϣ�������� \n');
            
            %��¼ residuOri
            
            
            while  residuOri > initialLog.residu
                
                %��ʱģʽ����
                if initialLog.timeLimitFlag == 1
                    tempCurrentTime = clock ;
                end
                
                %�ж��Ƿ�ʱ����������ѭ��
                if TimePast( tempStartTime , tempCurrentTime ) > timeLimitPerLoop && residuOri < initialLog.timeModeLeastResidu
                    fprintf('��ʱ����ʱresiduΪ %5.4f \n',residuOri ) ;
                    break ;
                end
                
                iterCount = iterCount + 1 ;
                
                %ai ԭ��Ϊalpha���뼫���ʳ�ͻ����ѡ��ai�� ai�������ǲ�����ϵ����gigi ��Ϊgi���ڻ���
                gigi = gx(:)' * gx(:) + gy(:)' * gy(:) + gz(:)' *gz(:) ;
                ai = gigi / (vx(:)' * vx(:) + vy(:)' * vy(:) + vz(:)' *vz(:) );
                
                %û�е�ĵط���PΪ0
                Px = Px + ai * px ;
                Py = Py + ai * py ;
                Pz = Pz + ai * pz ;
                
                wx = wx + ai * vx ;
                wy = wy + ai * vy ;
                wz = wz + ai * vz ;
                
                gx = Zx - DDA_ConvAccelerate( Axxfg , wx , 1) - DDA_ConvAccelerate( Axyfg , wy , 1) - DDA_ConvAccelerate( Axzfg , wz , 1) - wx / conj(alphaLDR) ;
                gy = Zy - DDA_ConvAccelerate( Axyfg , wx , 1) - DDA_ConvAccelerate( Ayyfg , wy , 1) - DDA_ConvAccelerate( Ayzfg , wz , 1) - wy / conj(alphaLDR) ;
                gz = Zz - DDA_ConvAccelerate( Axzfg , wx , 1) - DDA_ConvAccelerate( Ayzfg , wy , 1) - DDA_ConvAccelerate( Azzfg , wz , 1) - wz / conj(alphaLDR) ;
                
                %ȥ��Ч��
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
                
                %ȥ��Ч��
                vx = vx .* structure ;
                vy = vy .* structure ;
                vz = vz .* structure ;
                
                %ÿ5�ε�������һ��residuOri
                %residuOri
                if mod( iterCount , 5 ) == 1
                    
                    residuOri = ((Ex(:) - wx(:))'*(Ex(:) - wx(:)) + (Ey(:) - wy(:))'*(Ey(:) - wy(:)) + (Ez(:) - wz(:))'*(Ez(:) - wz(:))) / EiEi0 ;
                    
                    fprintf('�������� %i ��ǰԭ����residuΪ %d \n', iterCount , residuOri);
                    
                end
            end
            %��¼��residuOri
            residuOriArray = [residuOriArray , residuOri] ;
            fprintf('�����һ�ε������㡣��ǰt_thetaΪ %d ����Ϊ %d ����ǰ�糡����Ϊ %d �����t_theta����Ϊ %d \n',t_thetaValueRad,j,eDirection,t_thetaMaxOrder) ;
            
            %����P��P��������ΪPx Py Pz
            P = zeros(Nx + 1, Ny + 1, Nz + 1,3,'single') ;
            P(:, :, : ,1) = Px ;
            P(:, :, : ,2) = Py ;
            P(:, :, : ,3) = Pz ;
            
            %��P��ֵ����scaData
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


