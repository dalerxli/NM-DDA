function coePosArray = DDA_PointMullerCoeGener( tempMuller , L_theta , L_phi )
%DDA_POINTMULLERCOEGENER
% ������20170424�������������
% �˺���������ϵ�����˺�������DDA_PointMuller��ʹ��
% �˺������븴���ʺܵͣ�����һ��ʼ�ʹ���д����Щʵ�֣����ԾͲ������е�һЩ���̳���ɺ�����

% ������thetaArray �Լ�ÿ�� L_thetaOrder(iter).phiArray
thetaArray = unique( tempMuller(:,1) ) ;
%markNum ���ڱ�ʾ
markNum = 0 ;
for iter = 1 : length( thetaArray )
    phiIndex = find( tempMuller(:,1) == thetaArray( iter )) ;
    L_thetaOrder(iter).phiArray = tempMuller( phiIndex , 2 ) ;
    L_thetaOrder(iter).markNum = markNum ;
    markNum = markNum + length( L_thetaOrder(iter).phiArray ) ;
end

thetaInter = abs( thetaArray(2) - thetaArray(1) );
posTheta = find( thetaArray == L_theta ) ;

%�鿴�Ƿ���ֱ���ҵ�theta
if isempty( posTheta )
    
    %��ʱposTheta��������
    posTheta = find( abs(thetaArray - L_theta) < thetaInter) ;
    disThetaLA = abs(thetaArray - L_theta) ;
    coe1 = disThetaLA( posTheta(2) ) /( disThetaLA( posTheta(1) ) + disThetaLA( posTheta(2) )) ;
    coe2 = disThetaLA( posTheta(1) ) /( disThetaLA( posTheta(1) ) + disThetaLA( posTheta(2) )) ;
    
    %��posTheta(1)���д���
    posPhi = find( L_thetaOrder( posTheta(1) ).phiArray == L_phi ) ;
    if isempty( posPhi )
        %phi�ļ��
        phiInter = 360 / length( L_thetaOrder( posTheta(1) ).phiArray ) ;
        
        %����L_phi���phiArray��֮��ľ��룬ע��Բ���ϵľ�����������㣬Ҳ����������
        for iter = 1 : length( L_thetaOrder( posTheta(1) ).phiArray )
            if L_thetaOrder( posTheta(1) ).phiArray( iter ) >= L_phi
                arrayRight( iter ) = L_thetaOrder( posTheta(1) ).phiArray( iter ) - L_phi ;
                arrayLeft( iter ) = 360 - ( L_thetaOrder( posTheta(1) ).phiArray( iter ) - L_phi ) ;
            else
                arrayRight( iter ) = 360 + L_thetaOrder( posTheta(1) ).phiArray( iter ) - L_phi ;
                arrayLeft( iter ) = L_phi - L_thetaOrder( posTheta(1) ).phiArray( iter ) ;
            end
        end
        
        disLA = zeros(1, length(arrayRight)) ;
        for iter = 1 : length( L_thetaOrder( posTheta(1) ).phiArray )
            disLA( iter ) = min( arrayRight( iter ) , arrayLeft( iter )) ;
        end
        
        %�ҵ�����������������
        posPhi = find( disLA <  phiInter ) ;
        indice11 = posPhi(1) + L_thetaOrder( posTheta(1) ).markNum ;
        indice12 = posPhi(2) + L_thetaOrder( posTheta(1) ).markNum ;
        coe11 =  disLA(posPhi(2)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
        coe12 =  disLA(posPhi(1)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
       
        %ֱ���ҵ���phi�����
    else
        indice11 = posPhi + L_thetaOrder( posTheta(1) ).markNum ;
        indice12 = indice11 ;
        coe11 = 1 ;
        coe12 = 0 ;
    end
    
    %��posTheta(2)���д���
    posPhi = find( L_thetaOrder( posTheta(2) ).phiArray == L_phi ) ;
    if isempty( posPhi )
        %phi�ļ��
        phiInter = 360 / length( L_thetaOrder( posTheta(2) ).phiArray ) ;
        
        %����L_phi���phiArray��֮��ľ��룬ע��Բ���ϵľ�����������㣬Ҳ����������
        for iter = 1 : length( L_thetaOrder( posTheta(2) ).phiArray )
            if L_thetaOrder( posTheta(2) ).phiArray( iter ) >= L_phi
                arrayRight( iter ) = L_thetaOrder( posTheta(2) ).phiArray( iter ) - L_phi ;
                arrayLeft( iter ) = 360 - ( L_thetaOrder( posTheta(2) ).phiArray( iter ) - L_phi ) ;
            else
                arrayRight( iter ) = 360 + L_thetaOrder( posTheta(2) ).phiArray( iter ) - L_phi ;
                arrayLeft( iter ) = L_phi - L_thetaOrder( posTheta(2) ).phiArray( iter ) ;
            end
        end
        
        disLA = zeros(1, length(arrayRight)) ;
        for iter = 1 : length( L_thetaOrder( posTheta(2) ).phiArray )
            disLA( iter ) = min( arrayRight( iter ) , arrayLeft( iter )) ;
        end
        
        %�ҵ�����������������
        posPhi = find( disLA <  phiInter ) ;
        indice21 = posPhi(1) + L_thetaOrder( posTheta(2) ).markNum ;
        indice22 = posPhi(2) + L_thetaOrder( posTheta(2) ).markNum ;
        coe21 =  disLA(posPhi(2)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
        coe22 =  disLA(posPhi(1)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
       
        %ֱ���ҵ���phi�����
    else
        indice21 = posPhi + L_thetaOrder( posTheta(2) ).markNum ;
        indice22 = indice21 ;
        coe21 = 1 ;
        coe22 = 0 ;
    end
    
    %���� coePosArray 
    coePosArray = [indice11 coe11 * coe1 ;
        indice12 coe12 * coe1 ;
        indice21 coe21 * coe2 ;
        indice22 coe22 * coe2  ] ;
    
    % �ҵ���theta�ķ�֧
else
    %�鿴�Ƿ���ֱ���ҵ�phi
    posPhi = find( L_thetaOrder( posTheta ).phiArray == L_phi ) ;
    
    %�ҵ���theta ����û��ֱ���ҵ�phi�����
    if isempty( posPhi )
        %phi�ļ��
        phiInter = 360 / length( L_thetaOrder( posTheta ).phiArray ) ;
        
        %����L_phi���phiArray��֮��ľ��룬ע��Բ���ϵľ�����������㣬Ҳ����������
        for iter = 1 : length( L_thetaOrder( posTheta ).phiArray )
            if L_thetaOrder( posTheta ).phiArray( iter ) >= L_phi
                arrayRight( iter ) = L_thetaOrder( posTheta ).phiArray( iter ) - L_phi ;
                arrayLeft( iter ) = 360 - ( L_thetaOrder( posTheta ).phiArray( iter ) - L_phi ) ;
            else
                arrayRight( iter ) = 360 + L_thetaOrder( posTheta ).phiArray( iter ) - L_phi ;
                arrayLeft( iter ) = L_phi - L_thetaOrder( posTheta ).phiArray( iter ) ;
            end
        end
        
        disLA = zeros(1, length(arrayRight)) ;
        for iter = 1 : length( L_thetaOrder( posTheta ).phiArray )
            disLA( iter ) = min( arrayRight( iter ) , arrayLeft( iter )) ;
        end
        
        %�ҵ�����������������
        posPhi = find( disLA <  phiInter ) ;
        indice1 = posPhi(1) + L_thetaOrder( posTheta ).markNum ;
        indice2 = posPhi(2) + L_thetaOrder( posTheta ).markNum ;
        coe1 =  disLA(posPhi(2)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
        coe2 =  disLA(posPhi(1)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
        coePosArray = [ indice1 , coe1 ; indice2 , coe2 ] ;
        
        %ֱ���ҵ��� theta phi�����
    else
        indice1 = posPhi + L_thetaOrder( posTheta ).markNum ;
        coe1 = 1 ;
        coePosArray = [ indice1 , coe1 ];
    end
    
    
end

end

