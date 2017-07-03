function  outArray = DDA_ConvAccelerate( Afg , Y , AHermiFlag)

%�˺�������  ����Afg ���� ����Y
% Afg �Ѿ���fft�������Ѿ��������Դ�֮��
% Y Ϊԭ����
% ���㷨�ĺ�������fft���پ������ �� ��ϸ���۲μ� Goodman+Draine+Flatau 1991 ������£��Լ�wiki�ϵľ������

% A�Ĵ�СΪ 2Nx + 1, 2Ny + 1 ,2Nz + 1
% Y�Ĵ�СΪ Nx + 1, Ny + 1 ,Nz + 1
% ��׺g������gpuArray

%����һ����Afgһ�����Yg�� outArrayg��
persistent Yg ;
persistent outArrayg ;
if isempty(Yg)
    Yg = gpuArray.zeros(size(Afg),'single') ;
end
if isempty(outArrayg)
    outArrayg = gpuArray.zeros(size(Afg),'single') ;
end
%��Y��ֵ��ֵ��Yg ��ʹYg�����ĵط�Ϊ0 ע��Yg������ 7/8 �ĵط�Ϊ0Ԫ��
n1 = size(Y,1) ;
n2 = size(Y,2) ;
n3 = size(Y,3) ;

Yg(1:n1 , 1:n2 , 1:n3 ) = Y;

%ʹ�þ������ ����ȡ���
%AHermiFlag ��Ӧ��������һ�����������A��Hermi�ģ����㷨����Ҫ�õ�A'�ĳ˷����ʶ���Aֱ�ӳ���Y' ���Խ��ȡ������
%��������ᵽ��paper

if AHermiFlag == 0
    outArrayg = ifftn(Afg .* fftn(Yg)) ;
    outArray = gather(outArrayg(n1 : 2*n1 - 1 , n2 : 2*n2 - 1 , n3 : 2*n3 - 1)) ;
end

if AHermiFlag == 1
    outArrayg = conj( ifftn( Afg .* fftn( conj( Yg ))) );
    outArray = gather(outArrayg(n1 : 2*n1 - 1 , n2 : 2*n2 - 1 , n3 : 2*n3 - 1)) ;
end


