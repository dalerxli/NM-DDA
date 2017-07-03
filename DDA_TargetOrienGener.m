function [ scaData ] = DDA_TargetOrienGener( initialLog )
%DDA_TARGETORIENGENER
% 此函数用于生成scaData结构体

%生成 scaData 结构体
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
% 以上theta phi beta 的单位均为deg

%生成的结构体的样子取决于initialLog

%不考虑颗粒物取向时的情况（但是考虑颗粒物可能有内在取向）
if initialLog.rotationFlag == 0
    
    %不考虑取向时，t_thetaValue 为0
    scaData.t_thetaMaxOrder = 1 ;
    scaData.t_thetaOrder(1).t_thetaValue = 0 ;
    
    %phiMaxNum phiArray设置，不考虑取向时phi唯一、为0
    scaData.t_thetaOrder(1).phiMaxNum = 1 ;
    scaData.t_thetaOrder(1).phiArray = [ 0 ] ;
    
    %考虑内在旋转情况
    %无内在旋转时，将 betaMaxOrder 设为1
    if initialLog.innerRotationFlag == 0
        scaData.t_betaMaxOrder = 1 ;
        scaData.t_thetaOrder(1).t_betaOrder(1).t_betaValue = 0 ;
    end
    
    %有内在旋转时
    if initialLog.innerRotationFlag == 1
        
        scaData.t_betaMaxOrder = initialLog.innerRotationNum ;
        
        for i = 1 : scaData.t_betaMaxOrder
            scaData.t_thetaOrder(1).t_betaOrder(i).t_betaValue = 360 / scaData.t_betaMaxOrder * ( i - 1) ;
        end
        
    end
end

%以下为考虑颗粒物取向时的情况，已经假定颗粒物的取向分布为在球面上均匀分布
if initialLog.rotationFlag == 1
    
    %根据 rotationNum 计算出颗粒物的theta转动时的间隔 thetaDelta
    thetaDivideNum = round(( pi / 4 * initialLog.rotationNum )^0.5 ) ;
    thetaDelta = 180 / thetaDivideNum ;
    
    %thetaMaxOrder 为 thetaDivideNum 加1
    scaData.t_thetaMaxOrder = thetaDivideNum + 1 ;
    
    %没有内在旋转的情况
    if initialLog.innerRotationFlag == 0
        
        scaData.t_betaMaxOrder = 1 ;
        
        %用于计取向点的个数
        tempN = 0 ;
        
        for i = 1 : scaData.t_thetaMaxOrder
            
            %此时的thetaValue
            scaData.t_thetaOrder(i).t_thetaValue = thetaDelta * ( i - 1 ) ;
            
            %生成此时的phiArray
            %tempMaxNum代表此theta角下的phi取向个数，注意其有可能因为theta角为0而为0，故而需要至少为1
            tempMaxNum = round(2 * pi * sin( deg2rad( scaData.t_thetaOrder(i).t_thetaValue ) ) / deg2rad( thetaDelta )) ;
            tempMaxNum = max(tempMaxNum , 1) ;
            scaData.t_thetaOrder(i).phiMaxNum = tempMaxNum ;
            
            %tempN用于计量总共的取向数
            tempN = tempN + tempMaxNum ;
            scaData.t_thetaOrder(i).phiArray = [ 0 : tempMaxNum - 1 ] * 360 / tempMaxNum ;
            scaData.t_thetaOrder(i).t_betaOrder(1).t_betaValue = 0 ;
        end
        %取向有多少个点...，注意orieNum和一开始设定的rotationNum只是近似相等
        scaData.orieNum = tempN ;
    end
    
    %有内在旋转的情况下
    if initialLog.innerRotationFlag == 1
        
        scaData.t_betaMaxOrder = initialLog.innerRotationNum ;
        
        %用于计取向点的个数
        tempN = 0 ;
        
        for i = 1 : thetaDivideNum + 1
            
            %此时的thetaValue
            scaData.t_thetaOrder(i).t_thetaValue = thetaDelta * ( i - 1 ) ;
            %生成此时的phiArray
            tempMaxNum = round(2 * pi * sin( deg2rad( scaData.t_thetaOrder(i).t_thetaValue ) ) / deg2rad( thetaDelta )) ;
            tempMaxNum = max(tempMaxNum , 1) ;
            scaData.t_thetaOrder(i).phiMaxNum = tempMaxNum ;
            tempN = tempN + tempMaxNum ;
            scaData.t_thetaOrder(i).phiArray = [ 0 : tempMaxNum - 1 ] * 360 / tempMaxNum ;
            
            for j = 1 : scaData.t_betaMaxOrder
                scaData.t_thetaOrder(i).t_betaOrder(j).t_betaValue = 360 / scaData.t_betaMaxOrder * ( j - 1 ) ;
            end
            
        end
        %取向有多少个点...
        scaData.orieNum = tempN ;
        
    end
   
end



end

