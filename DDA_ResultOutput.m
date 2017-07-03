%clear all
clear all ;
profile on
%load
load('scaData.mat') ;
load('angleArray.mat') ;

%��ȡModel��initialLog
Model = scaData.Model ;
initialLog = scaData.initialLog ;

% ��Model�г�ʼ������
structure = Model.struc;
m = Model.m;
d = Model.d;

%�õ�Nx Ny Nz alpha
Nx = size(structure,1) - 1;
Ny = size(structure,2) - 1;
Nz = size(structure,3) - 1;

%lamda Ϊ�����������еĲ��� kΪ��ʸ
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

%t_theta��t_beta��forѭ��������һ����̬�µ�ɢ��ͼ����t_phi�ı仯���Կ��任���õ���������һ���ļ����и���
for t_thetaIter = 1 : t_thetaMaxOrder
    for t_betaIter = 1 : t_betaMaxOrder
        
        fprintf('���ڼ��� thetaIter : %d ; betaIter : %d \n',t_thetaIter ,t_betaIter ) ;
        
        %��ȡ��ʱ��horiP��vertiP��������ʱ����
        this_horiP = scaData.t_thetaOrder(t_thetaIter).t_betaOrder(t_betaIter).horiP ;
        this_vertiP = scaData.t_thetaOrder(t_thetaIter).t_betaOrder(t_betaIter).vertiP ;
        
        %��ȡ��ʱt_theta��t_beta��ֵ���������仡���Ƶ�ֵ�����������ҡ�����ֵ���ټ�����
        t_thetaValue = scaData.t_thetaOrder(t_thetaIter).t_thetaValue;
        t_betaValue = scaData.t_thetaOrder(t_thetaIter).t_betaOrder(t_betaIter).t_betaValue ;
        t_thetaValueRad = deg2rad( t_thetaValue ) ;
        t_betaValueRad = deg2rad( t_betaValue ) ;
        
        sinTtheta = sin(t_thetaValueRad) ;
        cosTtheta = cos(t_thetaValueRad) ;
        sinTbeta = sin(t_betaValueRad) ;
        cosTbeta = cos(t_betaValueRad) ;
        
        %��ʼ��Ҫ�洢��muller���� ��ǰ����Ϊ���ж�Ӧ�ĽǶȣ��Ƕ��ƣ�
        tempMuller = zeros(maxL_AangleNum,18) ;
        tempMuller(:,1:2) = angleArray(:,3:4) ;
        
        %��ʼ����Ӧ��Jones���� ����ΪS1 S2 S3 S4 ��Ӧ�ľ���ԪΪ a22 a11 a12 a21
        %��Absorption and Scattering of Light by Small Particles ��P73 ���Ӱ�
        S1 = zeros(maxL_AangleNum,1) ;
        S2 = zeros(maxL_AangleNum,1) ;
        S3 = zeros(maxL_AangleNum,1) ;
        S4 = zeros(maxL_AangleNum,1) ;
        
        %���㲿��
        %����̽�����ķ�λ���м���
        parfor L_AngleIter = 1 : maxL_AangleNum
            
            %��ȡҪ����ĵ�ķ�λ�ǣ������ƣ�
            L_thetaRad = angleArray(L_AngleIter , 1) ;
            L_phiRad = angleArray(L_AngleIter , 2) ;
            
            if mod(L_AngleIter,5) == 0
                fprintf('L_AngleIter now %d\n',L_AngleIter) ;
            end
            
            %������ˮƽ����ֱ����������ȡJones�����е�S1 S2 S3 S4
            horiUni = [cos(L_thetaRad)*cos(L_phiRad) , -sin(L_thetaRad) , -sin(L_phiRad) * cos(L_thetaRad) ] ;
            vertiUni = [sin(L_phiRad),0,cos(L_phiRad)] ;
            
            %�����ķ��� ������ת�ã��ӿ�����ٶ�
            rUni = [sin( L_thetaRad ) * cos( L_phiRad ) , cos( L_thetaRad ) , -sin( L_thetaRad ) * sin( L_phiRad ) ] ;
            rUniInver = rUni' ;
            %��λ����Ĳ�˾������ڼӿ����
            rcross = [0 -rUni(3) rUni(2) ; rUni(3) 0 -rUni(1) ; -rUni(2) rUni(1) 0 ] ;
            rcross2 = rcross * rcross ;
            
            %���²���Ϊѭ������
            horiE = 0 ;
            vertiE = 0 ;
            
            for nx = 0 : Nx
                for ny = 0 : Ny
                    for nz = 0 : Nz
                        
                        %���ڼ����ż���ӵ�λ�� pos1,2,3�ֱ����ż���ӵ�x y zλʸ����
                        pos = zeros(1,3) ;
                        pos(1) = nx * cosTbeta + nz * sinTbeta ;
                        pos(2) = -nx * sinTbeta * sinTtheta + ny * cosTtheta + nz * sinTtheta * cosTbeta ;
                        pos(3) = -nx * sinTbeta * cosTtheta - ny * sinTtheta + nz * cosTtheta * cosTbeta ;
                        
                        %��ȡ��ż���ӵļ�������
                        %                         tempHoriP(:) = squeeze(this_horiP(nx+1 , ny+1 , nz+1 ,: )) ;
                        %                         tempVeriP(:) = this_vertiP(nx+1 , ny+1 , nz+1 ,: ) ;
                        
                        %���㲿��
                        dot_PosrUni = pos * rUniInver ;
                        tempCoeToMulti =  exp(-1i * kd * dot_PosrUni) * rcross2 ;
                        
                        horiE =   horiE - tempCoeToMulti * squeeze(this_horiP(nx+1 , ny+1 , nz+1 ,: ))  ;
                        vertiE = vertiE - tempCoeToMulti * squeeze(this_vertiP(nx+1 , ny+1 , nz+1 ,: ))  ;
                    end
                end
            end
            
            %Absorption and Scattering of Light by Small Particles��P75���ң����Ӱ棩
            S2( L_AngleIter ) = horiUni * horiE ;
            S4( L_AngleIter ) = vertiUni * horiE ;
            
            S3( L_AngleIter ) = horiUni * vertiE ;
            S1( L_AngleIter ) = vertiUni * vertiE ;
            
        end
        
        %��Muller������Ϣ����tempMuller�У���˳��ΪS11 S12 S13 S14 S21 ...
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
        
        %��ȫphi
        %t_phiͨ����ת������ʵ��
        %���´��뽫t_phiOrder����
        
        phiArray = scaData.t_thetaOrder( t_thetaIter ).phiArray ;
        
        %���´��뽫t_phiOrder����
        phiMaxNum = scaData.t_thetaOrder( t_thetaIter ).phiMaxNum ;
        if phiMaxNum ~= 1
            for t_phiIter = 2 : phiMaxNum
                %�õ���Ҫ��ת��phi��
                t_phiDeg = phiArray( t_phiIter ) ;
                outPutMuller.t_thetaOrder( t_thetaIter ).t_phiOrder( t_phiIter ).t_betaOrder( t_betaIter ).Muller ...
                    = DDA_MullerT_PhiFillUp( t_phiDeg , tempMuller ) ;
            end
        end
    end
end
%��¼initialLog Model
outPutMuller.initialLog = initialLog ;
outPutMuller.Model = Model ;
save('outPutMuller.mat','outPutMuller');
%��ͼ����
%����x y z
x = sin(angleArray(:,1)) .* cos(angleArray(:,2)) ;
y = cos(angleArray(:,1)) ;
z = -sin(angleArray(:,1)) .* sin(angleArray(:,2)) ;
tempMuller = outPutMuller.thetaOrder( 1 ).t_phiOrder( 1 ).t_betaOrder( 1 ).Muller ;

maxNum = max(tempMuller(:,3)) ;

%���´��뽫t_phiOrder����

