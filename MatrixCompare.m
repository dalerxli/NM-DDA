function deltaPercent = MatrixCompare( standardMatrix , comparedMatrix )
%MATRIXCOMPARE 
% �˺�����������size��С��ͬ�ľ��� ���ؾ������İٷֱ�

%�����鲿��
sizeS = size(standardMatrix) ;
sizeC = size(comparedMatrix) ;

if length(sizeS) == length(sizeC)
    sizeCompareArray = sizeS == sizeC ;
    
    tempMulti = 1;
    for iter = 1 : length(sizeCompareArray)
        tempMulti = tempMulti * sizeCompareArray( iter ) ;
    end
    
    if tempMulti == 0
        fprintf('�������ĳߴ粻ͬ \n') ;
        return 
    end
end

%��ʼ�Ƚ���������
 
differMatrix = comparedMatrix - standardMatrix ;
absDM = abs( differMatrix ) ;
absSM = abs( standardMatrix ) ;

deltaPercent = sum( absDM(:) ) / sum( absSM(:) ) ;

end

