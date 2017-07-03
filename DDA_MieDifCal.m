function resArray = DDA_MieDifCal( outPutMuller , mieTempMuller ,elePos)
%DDA_MIEDIFCAL 
%这个函数用于比对DDA的输出和 mie解的差距

sumMie = sum(abs(mieTempMuller(:,3))) ;
sumMieele = sum(abs(mieTempMuller(:,elePos))) ;
resArray =[] ;
for iterTheta = 1 : size(outPutMuller.t_thetaOrder , 2)
    for iterPhi = 1 : size(outPutMuller.t_thetaOrder( iterTheta ).t_phiOrder , 2)
        for iterBeta = 1 : size(outPutMuller.t_thetaOrder( iterTheta ).t_phiOrder( iterPhi ).t_betaOrder , 2)
            tempMuller = outPutMuller.t_thetaOrder( iterTheta ).t_phiOrder( iterPhi ).t_betaOrder( iterBeta ).Muller ;
            tempMuller(:,3:18) = tempMuller(:,3:18) * sumMie / sum( abs( tempMuller(:,3) ) ) ;
            difMatirx = tempMuller - mieTempMuller ;
            tempDif = sum( abs( difMatirx(:,elePos) )) /sumMieele;
            resArray = [resArray ; iterTheta iterPhi iterBeta tempDif] ;
        end
    end
end



end

