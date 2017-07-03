function  sumAbsMuller = DDA_MullerElementsSignifValue( outPutMuller ,normalizedFlag)
%DDA_MULLERELEMENTSSIGNIFVALUE
%此函数用于估算一个outPutMuller结构体中 各个muller矩阵元的幅度大小

%行向量形式的求和 sumAbsMullerInRow
sumAbsMullerInRow = zeros(1,16) ;
%4X4矩阵形式的总和
sumAbsMuller = zeros(4,4) ;

%遍历所有取向
for iterTheta = 1 : size(outPutMuller.t_thetaOrder , 2)
    for iterPhi = 1 : size(outPutMuller.t_thetaOrder( iterTheta ).t_phiOrder , 2)
        for iterBeta = 1 : size(outPutMuller.t_thetaOrder( iterTheta ).t_phiOrder( iterPhi ).t_betaOrder , 2)
            
            %提取此时球面muller
            tempSphereMuller = outPutMuller.t_thetaOrder( iterTheta ).t_phiOrder( iterPhi ).t_betaOrder( iterBeta ).Muller(:,3:18) ;
            tempSphereAbsMullerInRow = sum( abs( tempSphereMuller )) ;
            sumAbsMullerInRow = tempSphereAbsMullerInRow + sumAbsMullerInRow ;
        end
    end
end

for rowIter = 1 : 4
    for columnIter = 1 : 4
        
        sumAbsMuller( rowIter , columnIter ) = sumAbsMullerInRow( rowIter * 4 + columnIter - 4 ) ;
        
    end
end

if normalizedFlag == 1
    sumAbsMuller = sumAbsMuller / sumAbsMuller(1,1) ;
end

end

