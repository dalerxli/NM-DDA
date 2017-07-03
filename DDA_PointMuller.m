function pointMuller = DDA_PointMuller( tempMuller , L_theta , L_phi )
%DDA_POINTMULLER
%ע���˺������о�̬�������ò��ķ�ʽ�ӿ��ٶ�
%�˺���ʵ�ִ�������muller���ҵ��ռ�ĳһ���muller����
%�˺�������ĺ��������ǣ�����tempMuller�������ǿռ��о��ȷֲ�����ɢ��
%�����ϣ���õ��ĵ���ԭ���ݵ�����û�У���ô��Ȼ��Ҫ��Χ�ĵ������ϵõ��˵�����

%���벿��
%��������һ�� tempMuller���������׵�muller��
%L_theta L_phiΪʵ��������ϵ������֪���ĵ�ļ����� �Ƕ���

%�������
%pointMuller һ��4 X 4 muller����

% tempCoe Ϊ�����ľ�̬������ ���ڴ���Ϊ�˵õ�ĳ���ض�ʵ����̽��Ƕ���Ҫ��ԭtempMuller������ϵ��
% �����2586���㡢  L_theta = 30 L_phi = 0 ����������� tempCoe����mark = [2586 30 0]�����ı�ǩ
% ���ڶ�Ӧ��coePosArray�и���ϵ��ֵ
% coePosArray�Ĺ������£� 
% [ ��һ��tempmuller������ ��һ��Ȩ�� ; 
%   ...         ...
%    ��N������             ��N��Ȩ�� ... ]
% coePosArray ���4��
% coePosArray ������Ȩ�غ�Ϊ1
% coePosArray �������� DDA_PointMullerCoeGener( tempMuller , L_theta , L_phi ) ʵ��

% �������岿��
% ���������ǰ�Ƿ�������Ӧ�Ƕȵ�ϵ��
% ��������δ����һ���������� temp_coePosArray 

persistent tempCoe ;
if isempty( tempCoe )
    tempCoe(1).mark = [ size( tempMuller , 1) , L_theta , L_phi ] ;
    tempCoe(1).coePosArray = DDA_PointMullerCoeGener( tempMuller , L_theta , L_phi ) ;
    temp_coePosArray = tempCoe(1).coePosArray ;
else
    temp_coePosArray = [];
    for iter = 1 : size( tempCoe , 2 )
        if tempCoe( iter ).mark(1) == size( tempMuller , 1) && ...
                tempCoe( iter ).mark(2) == L_theta && ...
                tempCoe( iter ).mark(3) == L_phi 
            temp_coePosArray = tempCoe( iter ).coePosArray ;
        end
    end
    
    if isempty( temp_coePosArray )
        newIndice =  size( tempCoe , 2 ) + 1 ;
        tempCoe( newIndice ).mark = [ size( tempMuller , 1) , L_theta , L_phi ] ;
        tempCoe( newIndice ).coePosArray = DDA_PointMullerCoeGener( tempMuller , L_theta , L_phi ) ;
        temp_coePosArray = tempCoe( newIndice ).coePosArray ;
    end
end

%����Χ������ݵ��� �õ��˵�����
tempPointMullerRow = zeros(1,18) ;
for iter = 1 : size( temp_coePosArray , 1 )
    tempPointMullerRow = tempPointMullerRow + temp_coePosArray( iter , 2 ) * tempMuller( temp_coePosArray( iter , 1 ) , : ) ;
end

%��ֵ��ֵ��pointMuller
pointMuller = zeros(4,4) ;
for iterRow = 1 : 4
    for iterColumn = 1 : 4
        pointMuller( iterRow , iterColumn ) = tempPointMullerRow( 4 * iterRow + iterColumn - 2 ) ;
    end
end

end



