function [ scaData ] = DDA_TargetOrienGener( initialLog )
%DDA_TARGETORIENGENER
% �˺�����������scaData�ṹ��

%���� scaData �ṹ��
%-------.t_thetaMaxOrder
%-------.t_betaMaxOrder
% ------.orieNum
% ------.t_thetaOrder
% ------------------.phiArray 
% ------------------.phiMaxNum
% ------------------.t_thetaValue
% ------------------.t_betaOrder
% -----------------------------.t_betaValue
% -----------------------------.horiP
%                              .vertiP
% ����theta phi beta �ĵ�λ��Ϊdeg

%���ɵĽṹ�������ȡ����initialLog

%�����ǿ�����ȡ��ʱ����������ǿ��ǿ��������������ȡ��
if initialLog.rotationFlag == 0
    
    %������ȡ��ʱ��t_thetaValue Ϊ0
    scaData.t_thetaMaxOrder = 1 ;
    scaData.t_thetaOrder(1).t_thetaValue = 0 ;
    
    %phiMaxNum phiArray���ã�������ȡ��ʱphiΨһ��Ϊ0
    scaData.t_thetaOrder(1).phiMaxNum = 1 ;
    scaData.t_thetaOrder(1).phiArray = [ 0 ] ;
    
    %����������ת���
    %��������תʱ���� betaMaxOrder ��Ϊ1
    if initialLog.innerRotationFlag == 0
        scaData.t_betaMaxOrder = 1 ;
        scaData.t_thetaOrder(1).t_betaOrder(1).t_betaValue = 0 ;
    end
    
    %��������תʱ
    if initialLog.innerRotationFlag == 1
        
        scaData.t_betaMaxOrder = initialLog.innerRotationNum ;
        
        for i = 1 : scaData.t_betaMaxOrder
            scaData.t_thetaOrder(1).t_betaOrder(i).t_betaValue = 360 / scaData.t_betaMaxOrder * ( i - 1) ;
        end
        
    end
end

%����Ϊ���ǿ�����ȡ��ʱ��������Ѿ��ٶ��������ȡ��ֲ�Ϊ�������Ͼ��ȷֲ�
if initialLog.rotationFlag == 1
    
    %���� rotationNum ������������thetaת��ʱ�ļ�� thetaDelta
    thetaDivideNum = round(( pi / 4 * initialLog.rotationNum )^0.5 ) ;
    thetaDelta = 180 / thetaDivideNum ;
    
    %thetaMaxOrder Ϊ thetaDivideNum ��1
    scaData.t_thetaMaxOrder = thetaDivideNum + 1 ;
    
    %û��������ת�����
    if initialLog.innerRotationFlag == 0
        
        scaData.t_betaMaxOrder = 1 ;
        
        %���ڼ�ȡ���ĸ���
        tempN = 0 ;
        
        for i = 1 : scaData.t_thetaMaxOrder
            
            %��ʱ��thetaValue
            scaData.t_thetaOrder(i).t_thetaValue = thetaDelta * ( i - 1 ) ;
            
            %���ɴ�ʱ��phiArray
            %tempMaxNum�����theta���µ�phiȡ�������ע�����п�����Ϊtheta��Ϊ0��Ϊ0���ʶ���Ҫ����Ϊ1
            tempMaxNum = round(2 * pi * sin( deg2rad( scaData.t_thetaOrder(i).t_thetaValue ) ) / deg2rad( thetaDelta )) ;
            tempMaxNum = max(tempMaxNum , 1) ;
            scaData.t_thetaOrder(i).phiMaxNum = tempMaxNum ;
            
            %tempN���ڼ����ܹ���ȡ����
            tempN = tempN + tempMaxNum ;
            scaData.t_thetaOrder(i).phiArray = [ 0 : tempMaxNum - 1 ] * 360 / tempMaxNum ;
            scaData.t_thetaOrder(i).t_betaOrder(1).t_betaValue = 0 ;
        end
        %ȡ���ж��ٸ���...��ע��orieNum��һ��ʼ�趨��rotationNumֻ�ǽ������
        scaData.orieNum = tempN ;
    end
    
    %��������ת�������
    if initialLog.innerRotationFlag == 1
        
        scaData.t_betaMaxOrder = initialLog.innerRotationNum ;
        
        %���ڼ�ȡ���ĸ���
        tempN = 0 ;
        
        for i = 1 : thetaDivideNum + 1
            
            %��ʱ��thetaValue
            scaData.t_thetaOrder(i).t_thetaValue = thetaDelta * ( i - 1 ) ;
            %���ɴ�ʱ��phiArray
            tempMaxNum = round(2 * pi * sin( deg2rad( scaData.t_thetaOrder(i).t_thetaValue ) ) / deg2rad( thetaDelta )) ;
            tempMaxNum = max(tempMaxNum , 1) ;
            scaData.t_thetaOrder(i).phiMaxNum = tempMaxNum ;
            tempN = tempN + tempMaxNum ;
            scaData.t_thetaOrder(i).phiArray = [ 0 : tempMaxNum - 1 ] * 360 / tempMaxNum ;
            
            for j = 1 : scaData.t_betaMaxOrder
                scaData.t_thetaOrder(i).t_betaOrder(j).t_betaValue = 360 / scaData.t_betaMaxOrder * ( j - 1 ) ;
            end
            
        end
        %ȡ���ж��ٸ���...
        scaData.orieNum = tempN ;
        
    end
   
end



end

