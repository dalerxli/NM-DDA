% ���ļ�Ϊ������־����¼���ָĶ�
% 20170405 �����DDA_Compute �� DDA_TargetOrienGener
% 
% DDA_TargetOrienGener ���Ծ�����������ϵ�ȡ��㣬�Լ����Ը����������������ת
% DDA_Compute ���ն�Ӧ�Ŀ�����ȡ������������P��PΪ����������

%0409 DDA_ResultOutput�ļ�顢�Ż�
%Ŀǰ������ ���¾�Ϊ�� 2586��̽���ļ�����
%����汾����Ҫ25860s���ҵ�ʱ���������
%�Ż�����
%1 ����vertiE �� horiEһ����ѭ���ڲ�������Ż�
%2 ���ɸ���temp���� �� ����ļ�����deg2rad squeeze �Ⱥ����ĵ��ô��� 
%��ʱ��ʱ�ӽ�3000s
%����parfor�Ż�
%����parfor�Ż� �������һ�ν����ʱ1000s���� ��cpu i5 6500�� 4�˻����£��ٶ�������������3�����ң���Ϊ����ͨ�ſ���

%0411 DDA_SphereFigure ������� ������ʾ����muller����Ԫ��������
%DDA_SphereFigure �����Ż�����֤����ͬ��С��Ԫ����ʾ����ɫΪͬһ��

%0412 DDA_ResultOutput���Ż�
%���þ�����㡢����˷��İ취�����Ż�
%��д�������µĽ������ű� DDA_ResultOutputMatrixVersion �� DDA_ResultOutputParforMatrixVersion
%�����ű��������˾�������ķ��������ټ���
%���ּ���������ʱ��
%DDA_ResultOutputMatrixVersion 14.0s
%DDA_ResultOutputParforMatrixVersion 12.5s
%�����۲췢�֣�����û��parfor��DDA_ResultOutputMatrixVersion���Զ����ö�ˣ�ʹ����Ч����parfor�汾һ��
%����������DDA_ResultOutput��10e-6������ɴ�ԭ��Ľ��������DDA_ResultOutput��forѭ���е�round off����

%0413 ResultOutput ����������GPU�㷨 ResultOutputGPUVersion �� ResultOutputParforGPUVersion 
%����Դ�����ʱ��ResultOutputGPUVersion ����ʱ��Ϊ DDA_ResultOutputMatrixVersion ������֮һ
%ResultOutputParforGPUVersion ����ʱ��Ϊ DDA_ResultOutputMatrixVersion ��ʮ����֮һ
%���Ͻű�������֤�� ����������DDA_ResultOutputһ��

% ��DDA_ResultOutput������Ч��Ϊ1 ���������ű�����Ч������
% DDA_ResultOutputMatrixVersion = DDA_ResultOutputParforMatrixVersion  = 80
% ResultOutputGPUVersion = 240
% ResultOutputParforGPUVersion = 960

%0416 ʹ�������mieɢ������֤����DDA����ȷ��
%���ɢ��������������� S11 S22 S23 S32 S33 S44Ϊ��Ҫֵ
%��Ҫֵ�ϣ�mieɢ��Ľ����DDA�Ľ�������Ƶ� 
%deltaPercent = MatrixCompare( mieTempMuller(:,indice), tempMuller(:,indice) ) ����
%���صĽ����ʾ��deltaPercentһ�����Ϊ0.2���ң���˵������Ľ���ǶԵģ��Լ�׫д��DDA��û������ġ�
%����0.2Ҳ�ǲ�С�����֣���˵�����DDA�Ľ����������ʵ�����һ������

%0417 ��������DDA_MullerT_PhiFillUp ������ɢ����֤����ȷ��
%��2586��̽��������£�������ת�����µ�ƫ��Ϊ 10e-8
%���⣬�������ӵ���תt_phi������̽��ǵ���תL_phi�������Լ������䷽����Ϊ��ʸ����Ϊ������

%�޸���DDA_Compute�� ��ת����deg��rad��������� ���������ȡ������ת��ת�Ľ�� ������Ϊ��������û������

%0419 ��ʼ�޸�DDA_Compute DDA_ConvAccelerate  �����ǰ� Ag ֱ���ȴ洢���Դ��� ���⽫A���ڴ洫�����Դ���
%�����޸�֮��DDA_ConvAccelerate  �Ľӿڴӽ���Af ��Ϊ����Afg �ʶ�DDA_Compute��Ҫ����Ӧ���޸�
%��outArrayg �� Yg д�����Դ���˺���Ч�ʵõ���2.5��������
%�޸ĺ�ĺ�����ԭ����e-4�Ĳ�𣬵�������ɢ����֤������ȷ��

%�� A�����ɷ��� eDirection ֮�� �� ��ʡ����A��ʱ�� �����߼�������

% 0420 ������DDA_MemoryEstimate��������������Ҫ���ĵ��ڴ桢�Դ�
% ������ 4.5 2.2 2.2um���ӵļ����ٶ� �� ��ģ�������DDSCAT��10������
% ��д�� DDA_ModelGenerate �����ڴ˺����Ĺ�����ֱ���������ӵļ����߶ȣ�Ȼ������Ӧ���ɵ���

% 0421
% ��д�� DDA_StrucPlot �˺������������汾��ModelGenerate�����е��ã���������Model֮�󣬽����ӵ�
% �ռ伸��ͼ����������

% 0424
% ����� DDA_PointMuller����������ͨ���������������������muller����