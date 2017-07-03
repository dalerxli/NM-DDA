function [ sumTempHoriArray , sumTempVertiArray] = DDA_ScaDataVerify( scaDataOri , scaDataNew )
%DDA_SCADATAVERIFY
%����������ڼ����д��DDA_Compute �������Ƿ���ȷ

%����������Ϊ��ͬ��scaData
%���û��޸���DDA_Compute�ڵ����ݣ�������õĺ�����Ϊ����֤�µĽ���Ƿ���ȷ��
%��ͬһModel��initialLog ���ɰ汾DDA���°汾DDA���ɵ�scaData���жԱȼ���

%���������
%���������Ϊ���� N X 3 ���飬��1��2��3�зֱ�Ϊ , ��ǰ iterTheta ��ֵ��
% ��ǰ iterBeta ��ֵ����ǰ�Ƕ�������P�����İٷֱ�
%�������ݷֱ���ˮƽP�� �� ��ֱP��

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
