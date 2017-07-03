function coePosArray = DDA_PointMullerCoeGener( tempMuller , L_theta , L_phi )
%DDA_POINTMULLERCOEGENER
% 函数于20170424测过，运行良好
% 此函数用于找系数，此函数仅在DDA_PointMuller中使用
% 此函数代码复用率很低，但是一开始就打算写死这些实现，所以就不对其中的一些过程抽象成函数了

% 先生成thetaArray 以及每个 L_thetaOrder(iter).phiArray
thetaArray = unique( tempMuller(:,1) ) ;
%markNum 用于表示
markNum = 0 ;
for iter = 1 : length( thetaArray )
    phiIndex = find( tempMuller(:,1) == thetaArray( iter )) ;
    L_thetaOrder(iter).phiArray = tempMuller( phiIndex , 2 ) ;
    L_thetaOrder(iter).markNum = markNum ;
    markNum = markNum + length( L_thetaOrder(iter).phiArray ) ;
end

thetaInter = abs( thetaArray(2) - thetaArray(1) );
posTheta = find( thetaArray == L_theta ) ;

%查看是否能直接找到theta
if isempty( posTheta )
    
    %此时posTheta有两个数
    posTheta = find( abs(thetaArray - L_theta) < thetaInter) ;
    disThetaLA = abs(thetaArray - L_theta) ;
    coe1 = disThetaLA( posTheta(2) ) /( disThetaLA( posTheta(1) ) + disThetaLA( posTheta(2) )) ;
    coe2 = disThetaLA( posTheta(1) ) /( disThetaLA( posTheta(1) ) + disThetaLA( posTheta(2) )) ;
    
    %对posTheta(1)进行处理
    posPhi = find( L_thetaOrder( posTheta(1) ).phiArray == L_phi ) ;
    if isempty( posPhi )
        %phi的间隔
        phiInter = 360 / length( L_thetaOrder( posTheta(1) ).phiArray ) ;
        
        %考虑L_phi点和phiArray点之间的距离，注意圆环上的距离可以右旋算，也可以左旋算
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
        
        %找到最近的两个点的坐标
        posPhi = find( disLA <  phiInter ) ;
        indice11 = posPhi(1) + L_thetaOrder( posTheta(1) ).markNum ;
        indice12 = posPhi(2) + L_thetaOrder( posTheta(1) ).markNum ;
        coe11 =  disLA(posPhi(2)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
        coe12 =  disLA(posPhi(1)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
       
        %直接找到了phi的情况
    else
        indice11 = posPhi + L_thetaOrder( posTheta(1) ).markNum ;
        indice12 = indice11 ;
        coe11 = 1 ;
        coe12 = 0 ;
    end
    
    %对posTheta(2)进行处理
    posPhi = find( L_thetaOrder( posTheta(2) ).phiArray == L_phi ) ;
    if isempty( posPhi )
        %phi的间隔
        phiInter = 360 / length( L_thetaOrder( posTheta(2) ).phiArray ) ;
        
        %考虑L_phi点和phiArray点之间的距离，注意圆环上的距离可以右旋算，也可以左旋算
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
        
        %找到最近的两个点的坐标
        posPhi = find( disLA <  phiInter ) ;
        indice21 = posPhi(1) + L_thetaOrder( posTheta(2) ).markNum ;
        indice22 = posPhi(2) + L_thetaOrder( posTheta(2) ).markNum ;
        coe21 =  disLA(posPhi(2)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
        coe22 =  disLA(posPhi(1)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
       
        %直接找到了phi的情况
    else
        indice21 = posPhi + L_thetaOrder( posTheta(2) ).markNum ;
        indice22 = indice21 ;
        coe21 = 1 ;
        coe22 = 0 ;
    end
    
    %生成 coePosArray 
    coePosArray = [indice11 coe11 * coe1 ;
        indice12 coe12 * coe1 ;
        indice21 coe21 * coe2 ;
        indice22 coe22 * coe2  ] ;
    
    % 找到了theta的分支
else
    %查看是否能直接找到phi
    posPhi = find( L_thetaOrder( posTheta ).phiArray == L_phi ) ;
    
    %找到了theta 但是没有直接找到phi的情况
    if isempty( posPhi )
        %phi的间隔
        phiInter = 360 / length( L_thetaOrder( posTheta ).phiArray ) ;
        
        %考虑L_phi点和phiArray点之间的距离，注意圆环上的距离可以右旋算，也可以左旋算
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
        
        %找到最近的两个点的坐标
        posPhi = find( disLA <  phiInter ) ;
        indice1 = posPhi(1) + L_thetaOrder( posTheta ).markNum ;
        indice2 = posPhi(2) + L_thetaOrder( posTheta ).markNum ;
        coe1 =  disLA(posPhi(2)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
        coe2 =  disLA(posPhi(1)) /( disLA(posPhi(1)) + disLA(posPhi(2)) )  ;
        coePosArray = [ indice1 , coe1 ; indice2 , coe2 ] ;
        
        %直接找到了 theta phi的情况
    else
        indice1 = posPhi + L_thetaOrder( posTheta ).markNum ;
        coe1 = 1 ;
        coePosArray = [ indice1 , coe1 ];
    end
    
    
end

end

