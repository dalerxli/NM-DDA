function deltaPercent = MatrixCompare( standardMatrix , comparedMatrix )
%MATRIXCOMPARE 
% 此函数接收两个size大小相同的矩阵 返回矩阵相差的百分比

%输入检查部分
sizeS = size(standardMatrix) ;
sizeC = size(comparedMatrix) ;

if length(sizeS) == length(sizeC)
    sizeCompareArray = sizeS == sizeC ;
    
    tempMulti = 1;
    for iter = 1 : length(sizeCompareArray)
        tempMulti = tempMulti * sizeCompareArray( iter ) ;
    end
    
    if tempMulti == 0
        fprintf('输入矩阵的尺寸不同 \n') ;
        return 
    end
end

%开始比较两个矩阵
 
differMatrix = comparedMatrix - standardMatrix ;
absDM = abs( differMatrix ) ;
absSM = abs( standardMatrix ) ;

deltaPercent = sum( absDM(:) ) / sum( absSM(:) ) ;

end

