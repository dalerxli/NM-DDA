%配置文件
profile on

%文件路径部分
DDA_FILEDIRCTRL

%请手动载入一个Model ...

%程序设置部分
%迭代允许的最大残差
initialLog.residu = 1e-4;

% 是否启用定时模式？启用后程序会在desiredTimeCost的时间内运行完，但一般精度会有影响
% desiredTimeCost的单位为s
initialLog.timeLimitFlag = 1 ;
initialLog.desiredTimeCost = 2 * 3600 ;
% 此模式下至少要达到的精度
initialLog.timeModeLeastResidu = 0.01 ; 

%入射光波长
initialLog.lambda = 0.532;

%颗粒物取向部分
%initialLog.targetOrienGener 用于选取使用哪个取向生成函数
%现有如下选择 ： 'default' 'ma0414'
%'default' 一般使用的函数，可以生成球面上的均匀取向
%'ma0414' 马老师在0414提出用于采集特定角度下柱、椭球的信号，由于柱只在一个面上有信号，
%故此取向分布在此面内旋转
initialLog.targetOrienGener = 'default' ;

%颗粒物取向
initialLog.rotationFlag = 1;
initialLog.rotationNum = 2000;

%颗粒物内在旋转取向 
initialLog.innerRotationFlag = 0;
initialLog.innerRotationNum = 3;

% RUN 
% 参数设置到此为止，以下为计算部分
[scaData , residuOriArray ] = DDA_Compute( Model,initialLog ) ;

%保存scaData
save([scaDataPath,'scaData_',Model.fileName,'.mat'],'scaData') ;
