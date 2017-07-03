function [ sumTempHoriArray , sumTempVertiArray] = DDA_ScaDataVerify( scaDataOri , scaDataNew )
%DDA_SCADATAVERIFY
%这个函数用于检验改写的DDA_Compute 输出结果是否正确

%函数的输入为不同的scaData
%若用户修改了DDA_Compute内的内容，及其调用的函数，为了验证新的结果是否正确，
%对同一Model、initialLog 将旧版本DDA与新版本DDA生成的scaData进行对比即可

%函数的输出
%函数的输出为两组 N X 3 数组，第1、2、3列分别为 , 当前 iterTheta 的值，
% 当前 iterBeta 的值，当前角度下整体P阵相差的百分比
%两组数据分别是水平P阵 和 竖直P阵

t_thetaMaxOrder = scaDataOri.t_thetaMaxOrder ;
t_betaMaxOrder = scaDataOri.t_betaMaxOrder ;

sumTempHoriArray =[];
sumTempVertiArray = [];

for iterTheta = 1 : t_thetaMaxOrder
    for iterBeta = 1 : t_betaMaxOrder
        horiPori =  scaDataOri.t_thetaOrder( iterTheta ).t_betaOrder( iterBeta ).horiP ;
        vertiPori = scaDataOri.t_thetaOrder( iterTheta ).t_betaOrder( iterBeta ).vertiP ;
        
        horiPnew =  scaDataNew.t_thetaOrder( iterTheta ).t_betaOrder( iterBeta ).horiP ;
        vertiPnew = scaDataNew.t_thetaOrder( iterTheta ).t_betaOrder( iterBeta ).vertiP ;
        
        horiTemp = horiPori - horiPnew ;
        vertiTemp = vertiPori - vertiPnew ;
        
        sumHoriP = sum(abs(horiPori(:))) ;
        sumVertiP = sum(abs(vertiPori(:))) ;
        
        sumTempHori = sum( abs( horiTemp(:) ) ) ;
        sumTempVerti = sum( abs( vertiTemp(:) ) ) ;
        
        sumTempHoriArray = [sumTempHoriArray ; iterTheta , iterBeta ,sumTempHori / sumHoriP] ;
        sumTempVertiArray = [sumTempVertiArray ;iterTheta , iterBeta , sumTempVerti / sumVertiP] ;
    end
end
