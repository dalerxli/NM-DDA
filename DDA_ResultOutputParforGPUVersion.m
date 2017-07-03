%此方法使用了如下加速手段：
%1.用矩阵乘法替代计算horiE时对偶极子的for循环
%2.使用了GPU加速计算
%3.使用了parfor
% 相对于 DDA_ResultOutputGPUVersion 此方法在遍历探测点的地方使用了parfor

%此方法在40*40*40对2586方位计算的时长是 2.9s 相对于DDA_ResultOutputGPUVersion 的17s运行时间是其五分之一
%此方法在63*63*63对10318方位的计算时长是 17.5s 相对于DDA_ResultOutputGPUVersion，运行时间是其四分之一
%（获得上述结果是在 i5 6500 平台上，所以此时应该是CPU的调度被充分利用了，并且显卡仍不存在瓶颈）
%此方法与DDA_ResultOutputGPUVersion比对，误差为0 故与原方法的误差在10e-6内

profile on
%文件路径部分
DDA_FILEDIRCTRL

%手动读取scaData ?
scaAutoLoadFlag = 0 ;

%load
if scaAutoLoadFlag == 1
    load('scaData.mat') ;
end
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
        %获取Px Py Pz
        horiPx = this_horiP(:,:,:,1) ;
        horiPy = this_horiP(:,:,:,2) ;
        horiPz = this_horiP(:,:,:,3) ;
        
        vertiPx = this_vertiP(:,:,:,1) ;
        vertiPy = this_vertiP(:,:,:,2) ;
        vertiPz = this_vertiP(:,:,:,3) ;
        
        %得到posx posy posz
        for nx = 1 : Nx + 1
            for ny = 1 : Ny + 1
                for nz = 1 : Nz + 1
                    posx( nx , ny , nz ) = nx * cosTbeta + nz * sinTbeta ;
                    posy( nx , ny , nz ) = -nx * sinTbeta * sinTtheta + ny * cosTtheta + nz * sinTtheta * cosTbeta ;
                    posz( nx , ny , nz ) = -nx * sinTbeta * cosTtheta - ny * sinTtheta + nz * cosTtheta * cosTbeta ;
                end
            end
        end
        
        %生成GPU上的矩阵
        %获取horiPx的size
        strucSize = size( horiPx ) ;
        horiPxg = gpuArray( horiPx ) ;
        horiPyg = gpuArray( horiPy ) ;
        horiPzg = gpuArray( horiPz ) ;
        
        vertiPxg = gpuArray( vertiPx ) ;
        vertiPyg = gpuArray( vertiPy ) ;
        vertiPzg = gpuArray( vertiPz ) ;
        
        %注意将posxg ...等设为single变量，节约显存以及加快计算
        posxg = gpuArray( single( posx ) ) ;
        posyg = gpuArray( single( posy ) ) ;
        poszg = gpuArray( single( posz ) ) ;
        
        %预先生成horiExg ...
        horiExg = gpuArray.zeros( maxL_AangleNum , 1 , 'single' ) ;
        horiEyg = gpuArray.zeros( maxL_AangleNum , 1 , 'single' ) ;
        horiEzg = gpuArray.zeros( maxL_AangleNum , 1 , 'single' ) ;
        
        vertiExg = gpuArray.zeros( maxL_AangleNum , 1 , 'single' ) ;
        vertiEyg = gpuArray.zeros( maxL_AangleNum , 1 , 'single' ) ;
        vertiEzg = gpuArray.zeros( maxL_AangleNum , 1 , 'single' ) ;
        
        %L_thetaRadArray L_phiRadArray horiUni vertiUni
        L_thetaRadArray = single(angleArray( : , 1)) ;
        L_phiRadArray = single(angleArray( : , 2)) ;
        horiUni = [ cos( L_thetaRadArray ) .* cos( L_phiRadArray ) , -sin( L_thetaRadArray ) , -sin( L_phiRadArray ) .* cos( L_thetaRadArray ) ] ;
        vertiUni = [sin( L_phiRadArray ) , zeros( length( L_phiRadArray ) , 1) , cos( L_phiRadArray )] ;
        
        %遍历探测器的方位进行计算
        parfor L_AngleIter = 1 : maxL_AangleNum
            
            %提取要计算的点的方位角（弧度制）
            L_thetaRad = L_thetaRadArray( L_AngleIter ) ;
            L_phiRad = L_phiRadArray( L_AngleIter ) ;
            
            if mod( L_AngleIter , 100 ) == 0
                fprintf('L_AngleIter now %d\n',L_AngleIter) ;
            end
            
            %出射光的方向 给出其转置，加快计算速度
            rUni = [sin( L_thetaRad ) * cos( L_phiRad ) , cos( L_thetaRad ) , -sin( L_thetaRad ) * sin( L_phiRad ) ] ;
            
            %单位方向的叉乘矩阵，用于加快计算
            rcross = [0 -rUni(3) rUni(2) ; rUni(3) 0 -rUni(1) ; -rUni(2) rUni(1) 0 ] ;
            rcross2 = rcross * rcross ;
            
            %以下部分为矩阵算法计算部分，去掉了原来的循环方法
            %循环方法核心是计算 exp(-1i * kd * pos（xyz） * rUni') * rcross2 * P（xyz）的连加
            %上式中pos 和 P 都是随着nx ny nz而改变
            %此时如果采用matlab中的 .*运算 并求和，可能有更好的速度，此即此matrixVersion的思想
            
            %需要生成的pos和P会在进入 parfor L_AngleIter = 1 : maxL_AangleNum 循环前给出
            %生成中间矩阵Tx Ty Tz
            %生成 rcross2 乘以P之后的矩阵 T
            horiTxg = rcross2(1,1) * horiPxg + rcross2(1,2) * horiPyg + rcross2(1,3) * horiPzg ;
            horiTyg = rcross2(2,1) * horiPxg + rcross2(2,2) * horiPyg + rcross2(2,3) * horiPzg ;
            horiTzg = rcross2(3,1) * horiPxg + rcross2(3,2) * horiPyg + rcross2(3,3) * horiPzg ;
            
            vertiTxg = rcross2(1,1) * vertiPxg + rcross2(1,2) * vertiPyg + rcross2(1,3) * vertiPzg ;
            vertiTyg = rcross2(2,1) * vertiPxg + rcross2(2,2) * vertiPyg + rcross2(2,3) * vertiPzg ;
            vertiTzg = rcross2(3,1) * vertiPxg + rcross2(3,2) * vertiPyg + rcross2(3,3) * vertiPzg ;
            
            %rUni 先乘以 -1i * kd 得到前面系数 rUniC
            rUniC = -1i * kd * rUni ;
            
            %生成exp（）后的矩阵 tempexpPosrUni
            tempexpPosrUnig = exp ( posxg * rUniC(1) + posyg * rUniC(2) + poszg * rUniC(3) ) ;
            
            % tempexpPosrUni 和 T矩阵相乘后输出结果
            %生成horiE vertiE
            tempSumMatrix_HXg = tempexpPosrUnig .* horiTxg ;
            horiExg( L_AngleIter ) = sum( tempSumMatrix_HXg(:) ) ;
            
            tempSumMatrix_HYg = tempexpPosrUnig .* horiTyg ;
            horiEyg( L_AngleIter ) = sum( tempSumMatrix_HYg(:) ) ;
            
            tempSumMatrix_HZg = tempexpPosrUnig .* horiTzg ;
            horiEzg( L_AngleIter ) = sum( tempSumMatrix_HZg(:) ) ;
            
            tempSumMatrix_VXg = tempexpPosrUnig .* vertiTxg ;
            vertiExg( L_AngleIter ) = sum( tempSumMatrix_VXg(:) ) ;
            
            tempSumMatrix_VYg = tempexpPosrUnig .* vertiTyg ;
            vertiEyg( L_AngleIter ) = sum( tempSumMatrix_VYg(:) ) ;
            
            tempSumMatrix_VZg = tempexpPosrUnig .* vertiTzg ;
            vertiEzg( L_AngleIter ) = sum( tempSumMatrix_VZg(:) ) ;
            
            %             transHoriE = horiE' ;
            %             transVertiE = vertiE' ;
            %矩阵算法计算部分结束
            
        end
        
        %先提取结果到内存
        horiEx = gather( horiExg ) ;
        horiEy = gather( horiEyg ) ;
        horiEz = gather( horiEzg ) ;
        
        vertiEx = gather( vertiExg ) ;
        vertiEy = gather( vertiEyg ) ;
        vertiEz = gather( vertiEzg ) ;
        
        %Absorption and Scattering of Light by Small Particles书P75左右（电子版）
        S2 = horiUni( : , 1 ) .* horiEx + horiUni( : , 2 ) .* horiEy + horiUni( : , 3 ) .* horiEz ;
        S4 = vertiUni( : , 1 ) .* horiEx + vertiUni( : , 2 ) .* horiEy + vertiUni( : , 3 ) .* horiEz ;
        
        S3 = horiUni( : , 1 ) .* vertiEx + horiUni( : , 2 ) .* vertiEy + horiUni( : , 3 ) .* vertiEz ;
        S1 = vertiUni( : , 1 ) .* vertiEx + vertiUni( : , 2 ) .* vertiEy +vertiUni( : , 3 ) .* vertiEz ;
        
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
save([outPutMullerPath ,'outPutMuller_',Model.fileName,'.mat'],'outPutMuller');


