%�˷���ʹ�������¼����ֶΣ�
%1.�þ���˷��������horiEʱ��ż���ӵ�forѭ��
%2.ʹ����GPU���ټ���
%3.ʹ����parfor
% ����� DDA_ResultOutputGPUVersion �˷����ڱ���̽���ĵط�ʹ����parfor

%�˷�����40*40*40��2586��λ�����ʱ���� 2.9s �����DDA_ResultOutputGPUVersion ��17s����ʱ���������֮һ
%�˷�����63*63*63��10318��λ�ļ���ʱ���� 17.5s �����DDA_ResultOutputGPUVersion������ʱ�������ķ�֮һ
%���������������� i5 6500 ƽ̨�ϣ����Դ�ʱӦ����CPU�ĵ��ȱ���������ˣ������Կ��Բ�����ƿ����
%�˷�����DDA_ResultOutputGPUVersion�ȶԣ����Ϊ0 ����ԭ�����������10e-6��

profile on
%�ļ�·������
DDA_FILEDIRCTRL

%�ֶ���ȡscaData ?
scaAutoLoadFlag = 0 ;

%load
if scaAutoLoadFlag == 1
    load('scaData.mat') ;
end
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
        %��ȡPx Py Pz
        horiPx = this_horiP(:,:,:,1) ;
        horiPy = this_horiP(:,:,:,2) ;
        horiPz = this_horiP(:,:,:,3) ;
        
        vertiPx = this_vertiP(:,:,:,1) ;
        vertiPy = this_vertiP(:,:,:,2) ;
        vertiPz = this_vertiP(:,:,:,3) ;
        
        %�õ�posx posy posz
        for nx = 1 : Nx + 1
            for ny = 1 : Ny + 1
                for nz = 1 : Nz + 1
                    posx( nx , ny , nz ) = nx * cosTbeta + nz * sinTbeta ;
                    posy( nx , ny , nz ) = -nx * sinTbeta * sinTtheta + ny * cosTtheta + nz * sinTtheta * cosTbeta ;
                    posz( nx , ny , nz ) = -nx * sinTbeta * cosTtheta - ny * sinTtheta + nz * cosTtheta * cosTbeta ;
                end
            end
        end
        
        %����GPU�ϵľ���
        %��ȡhoriPx��size
        strucSize = size( horiPx ) ;
        horiPxg = gpuArray( horiPx ) ;
        horiPyg = gpuArray( horiPy ) ;
        horiPzg = gpuArray( horiPz ) ;
        
        vertiPxg = gpuArray( vertiPx ) ;
        vertiPyg = gpuArray( vertiPy ) ;
        vertiPzg = gpuArray( vertiPz ) ;
        
        %ע�⽫posxg ...����Ϊsingle��������Լ�Դ��Լ��ӿ����
        posxg = gpuArray( single( posx ) ) ;
        posyg = gpuArray( single( posy ) ) ;
        poszg = gpuArray( single( posz ) ) ;
        
        %Ԥ������horiExg ...
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
        
        %����̽�����ķ�λ���м���
        parfor L_AngleIter = 1 : maxL_AangleNum
            
            %��ȡҪ����ĵ�ķ�λ�ǣ������ƣ�
            L_thetaRad = L_thetaRadArray( L_AngleIter ) ;
            L_phiRad = L_phiRadArray( L_AngleIter ) ;
            
            if mod( L_AngleIter , 100 ) == 0
                fprintf('L_AngleIter now %d\n',L_AngleIter) ;
            end
            
            %�����ķ��� ������ת�ã��ӿ�����ٶ�
            rUni = [sin( L_thetaRad ) * cos( L_phiRad ) , cos( L_thetaRad ) , -sin( L_thetaRad ) * sin( L_phiRad ) ] ;
            
            %��λ����Ĳ�˾������ڼӿ����
            rcross = [0 -rUni(3) rUni(2) ; rUni(3) 0 -rUni(1) ; -rUni(2) rUni(1) 0 ] ;
            rcross2 = rcross * rcross ;
            
            %���²���Ϊ�����㷨���㲿�֣�ȥ����ԭ����ѭ������
            %ѭ�����������Ǽ��� exp(-1i * kd * pos��xyz�� * rUni') * rcross2 * P��xyz��������
            %��ʽ��pos �� P ��������nx ny nz���ı�
            %��ʱ�������matlab�е� .*���� ����ͣ������и��õ��ٶȣ��˼���matrixVersion��˼��
            
            %��Ҫ���ɵ�pos��P���ڽ��� parfor L_AngleIter = 1 : maxL_AangleNum ѭ��ǰ����
            %�����м����Tx Ty Tz
            %���� rcross2 ����P֮��ľ��� T
            horiTxg = rcross2(1,1) * horiPxg + rcross2(1,2) * horiPyg + rcross2(1,3) * horiPzg ;
            horiTyg = rcross2(2,1) * horiPxg + rcross2(2,2) * horiPyg + rcross2(2,3) * horiPzg ;
            horiTzg = rcross2(3,1) * horiPxg + rcross2(3,2) * horiPyg + rcross2(3,3) * horiPzg ;
            
            vertiTxg = rcross2(1,1) * vertiPxg + rcross2(1,2) * vertiPyg + rcross2(1,3) * vertiPzg ;
            vertiTyg = rcross2(2,1) * vertiPxg + rcross2(2,2) * vertiPyg + rcross2(2,3) * vertiPzg ;
            vertiTzg = rcross2(3,1) * vertiPxg + rcross2(3,2) * vertiPyg + rcross2(3,3) * vertiPzg ;
            
            %rUni �ȳ��� -1i * kd �õ�ǰ��ϵ�� rUniC
            rUniC = -1i * kd * rUni ;
            
            %����exp������ľ��� tempexpPosrUni
            tempexpPosrUnig = exp ( posxg * rUniC(1) + posyg * rUniC(2) + poszg * rUniC(3) ) ;
            
            % tempexpPosrUni �� T������˺�������
            %����horiE vertiE
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
            %�����㷨���㲿�ֽ���
            
        end
        
        %����ȡ������ڴ�
        horiEx = gather( horiExg ) ;
        horiEy = gather( horiEyg ) ;
        horiEz = gather( horiEzg ) ;
        
        vertiEx = gather( vertiExg ) ;
        vertiEy = gather( vertiEyg ) ;
        vertiEz = gather( vertiEzg ) ;
        
        %Absorption and Scattering of Light by Small Particles��P75���ң����Ӱ棩
        S2 = horiUni( : , 1 ) .* horiEx + horiUni( : , 2 ) .* horiEy + horiUni( : , 3 ) .* horiEz ;
        S4 = vertiUni( : , 1 ) .* horiEx + vertiUni( : , 2 ) .* horiEy + vertiUni( : , 3 ) .* horiEz ;
        
        S3 = horiUni( : , 1 ) .* vertiEx + horiUni( : , 2 ) .* vertiEy + horiUni( : , 3 ) .* vertiEz ;
        S1 = vertiUni( : , 1 ) .* vertiEx + vertiUni( : , 2 ) .* vertiEy +vertiUni( : , 3 ) .* vertiEz ;
        
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
save([outPutMullerPath ,'outPutMuller_',Model.fileName,'.mat'],'outPutMuller');


