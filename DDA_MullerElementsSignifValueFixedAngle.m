function sumAbsMuller = DDA_MullerElementsSignifValueFixedAngle( outPutMuller ,L_theta , L_phi , normalizedFlag)
%DDA_MULLERELEMENTSSIGNIFVALUEFIXEDANGLE 
%此函数用于估算一个outPutMuller结构体中 各个muller矩阵元的幅度大小
%此函数属于固定一个探测角观察的版本


%4X4矩阵形式的总和
sumAbsMuller = zeros(4,4) ;

for iterTheta = 1 : size(outPutMuller.t_thetaOrder , 2)
    for iterPhi = 1 : size(outPutMuller.t_thetaOrder( iterTheta ).t_phiOrder , 2)
        for iterBeta = 1 : size(outPutMuller.t_thetaOrder( iterTheta ).t_phiOrder( iterPhi ).t_betaOrder , 2)
            
            %提取此时球面muller
            tempSphereMuller = outPutMuller.t_thetaOrder( iterTheta ).t_phiOrder( iterPhi ).t_betaOrder( iterBeta ).Muller ;
            tempPointMuller = DDA_PointMuller( tempSphereMuller ,L_theta , L_phi) ;
            tempAbsPointMuller = abs( tempPointMuller ) ;
            sumAbsMuller = sumAbsMuller + tempAbsPointMuller ;
            
        end
    end
end

if normalizedFlag == 1
    sumAbsMuller = sumAbsMuller / sumAbsMuller(1,1) ;
end

end

