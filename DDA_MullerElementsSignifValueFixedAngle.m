function sumAbsMuller = DDA_MullerElementsSignifValueFixedAngle( outPutMuller ,L_theta , L_phi , normalizedFlag)
%DDA_MULLERELEMENTSSIGNIFVALUEFIXEDANGLE 
%�˺������ڹ���һ��outPutMuller�ṹ���� ����muller����Ԫ�ķ��ȴ�С
%�˺������ڹ̶�һ��̽��ǹ۲�İ汾


%4X4������ʽ���ܺ�
sumAbsMuller = zeros(4,4) ;

for iterTheta = 1 : size(outPutMuller.t_thetaOrder , 2)
    for iterPhi = 1 : size(outPutMuller.t_thetaOrder( iterTheta ).t_phiOrder , 2)
        for iterBeta = 1 : size(outPutMuller.t_thetaOrder( iterTheta ).t_phiOrder( iterPhi ).t_betaOrder , 2)
            
            %��ȡ��ʱ����muller
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

