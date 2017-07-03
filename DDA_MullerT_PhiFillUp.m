function phiFilledMuller = DDA_MullerT_PhiFillUp( t_phiDeg , tempMuller )
% DDA_MULLERT_PHIFILLUP
% ��������һ����������t_phi = 0ʱ��ʱ��muller�����ף����������� t_phiDeg ʱ�������ײ���
% ���յ���tempMuller ���Ϊ����ͬ��t_theta��t_beta�µ�muller�����ף�����t_phi���Ϊ0

tempPhiFilledMuller = zeros( size( tempMuller ) ) ;
tempPhiFilledMuller(:,1:2) = tempMuller(:,1:2) ;
phiFilledMuller = zeros( size( tempMuller ) ) ;
phiFilledMuller(:,1:2) = tempMuller(:,1:2) ;

t_phiRad = deg2rad( t_phiDeg ) ;
cos2phi = cos( 2 * t_phiRad ) ;
sin2phi = sin( 2 * t_phiRad ) ;

%�Ƚ�tempMuller���� ��theta
%markPos�м�¼��ÿ��theta��ʼ��λ��
markPos= [1] ;
tempValue = tempMuller(1,1) ;

for iter = 1 : size( tempMuller , 1 )
    if tempMuller( iter , 1 ) ~= tempValue
        markPos = [ markPos , iter ] ;
        tempValue = tempMuller( iter , 1 ) ;
    end
end

%�õ�phiNum���� �����д洢���Ƕ�ÿ��theta�����к��е�phi����ĸ���
for iter = 1 : length( markPos ) - 1
    
    phiNum( iter )= markPos( iter + 1) - markPos( iter ) ;
    
end

phiNum( length( markPos ) ) = size( tempMuller , 1 ) - markPos( length( markPos ) ) + 1 ;

% ����ѭ�� ��ÿ��theta������ֱ���
%baseRowNum 
baseRowNum = 0 ;
for iter = 1 : length( markPos )
    if phiNum( iter ) == 1
        tempPhiFilledMuller( baseRowNum + 1 , 3:18 ) = tempMuller( baseRowNum + 1 , 3:18 ) ;
        baseRowNum = baseRowNum + 1 ;
        
    else
        phiInter = 360 / phiNum( iter ) ;
        rotateIndice = floor( t_phiDeg / phiInter ) ;
        c1 = t_phiDeg / phiInter - rotateIndice ;
        c2 = 1 - c1 ;
        
        %����ѭ�������ÿһ��
        for phiIter = 1 : phiNum( iter ) 
            tempPhiFilledMuller( baseRowNum + phiIter , 3:18 ) = ...
                c1 * tempMuller( baseRowNum + DDA_Mod( phiIter - rotateIndice - 1 , phiNum( iter ) ) , 3 : 18 ) ...
                + c2 * tempMuller( baseRowNum + DDA_Mod( phiIter - rotateIndice , phiNum( iter ) ) , 3 : 18 ) ;
        end
        
        baseRowNum = baseRowNum + phiNum( iter ) ;
        
    end
end

%��ת����
%�����������ת����ϵ�����ʹ�����һ����ת������Ӧ��muller����
%���²μ� Absorbtion and Scattering ���P61�����Ӱ棩
% rotateMatrix = [ 1 0 0 0 ; 0 cos2phi -sin2phi 0 ; 0 sin2phi cos2phi 0 ; 0 0 0 1] ;    
phiFilledMuller(:,3) = tempPhiFilledMuller(:,3) ;
phiFilledMuller(:,4) = cos2phi * tempPhiFilledMuller(:,4) + sin2phi * tempPhiFilledMuller(:,5) ;
phiFilledMuller(:,5) = -sin2phi * tempPhiFilledMuller(:,4) + cos2phi * tempPhiFilledMuller(:,5) ;
phiFilledMuller(:,6) = tempPhiFilledMuller(:,6) ;

phiFilledMuller(:,7) = tempPhiFilledMuller(:,7) ;
phiFilledMuller(:,8) = cos2phi * tempPhiFilledMuller(:,8) + sin2phi * tempPhiFilledMuller(:,9) ;
phiFilledMuller(:,9) = -sin2phi * tempPhiFilledMuller(:,8) + cos2phi * tempPhiFilledMuller(:,9) ;
phiFilledMuller(:,10) = tempPhiFilledMuller(:,10) ;

phiFilledMuller(:,11) = tempPhiFilledMuller(:,11) ;
phiFilledMuller(:,12) = cos2phi * tempPhiFilledMuller(:,12) + sin2phi * tempPhiFilledMuller(:,13) ;
phiFilledMuller(:,13) = -sin2phi * tempPhiFilledMuller(:,12) + cos2phi * tempPhiFilledMuller(:,13) ;
phiFilledMuller(:,14) = tempPhiFilledMuller(:,14) ;

phiFilledMuller(:,15) = tempPhiFilledMuller(:,15) ;
phiFilledMuller(:,16) = cos2phi * tempPhiFilledMuller(:,16) + sin2phi * tempPhiFilledMuller(:,17) ;
phiFilledMuller(:,17) = -sin2phi * tempPhiFilledMuller(:,16) + cos2phi * tempPhiFilledMuller(:,17) ;
phiFilledMuller(:,18) = tempPhiFilledMuller(:,18) ;
end

