%�����ļ�
profile on

%�ļ�·������
DDA_FILEDIRCTRL

%���ֶ�����һ��Model ...

%�������ò���
%������������в�
initialLog.residu = 1e-4;

% �Ƿ����ö�ʱģʽ�����ú�������desiredTimeCost��ʱ���������꣬��һ�㾫�Ȼ���Ӱ��
% desiredTimeCost�ĵ�λΪs
initialLog.timeLimitFlag = 1 ;
initialLog.desiredTimeCost = 2 * 3600 ;
% ��ģʽ������Ҫ�ﵽ�ľ���
initialLog.timeModeLeastResidu = 0.01 ; 

%����Ⲩ��
initialLog.lambda = 0.532;

%������ȡ�򲿷�
%initialLog.targetOrienGener ����ѡȡʹ���ĸ�ȡ�����ɺ���
%��������ѡ�� �� 'default' 'ma0414'
%'default' һ��ʹ�õĺ������������������ϵľ���ȡ��
%'ma0414' ����ʦ��0414������ڲɼ��ض��Ƕ�������������źţ�������ֻ��һ���������źţ�
%�ʴ�ȡ��ֲ��ڴ�������ת
initialLog.targetOrienGener = 'default' ;

%������ȡ��
initialLog.rotationFlag = 1;
initialLog.rotationNum = 2000;

%������������תȡ�� 
initialLog.innerRotationFlag = 0;
initialLog.innerRotationNum = 3;

% RUN 
% �������õ���Ϊֹ������Ϊ���㲿��
[scaData , residuOriArray ] = DDA_Compute( Model,initialLog ) ;

%����scaData
save([scaDataPath,'scaData_',Model.fileName,'.mat'],'scaData') ;
