%此脚本用于测试DDA_MullerT_PhiFillUp函数是否正确
%此脚本的方法是用mieTempMuller作为输入，由于球在旋转时是不会变的，故而此时散射球面应该不会发生变化

%20170514对m=1.33 aeff=0.5 um lambda = 0.532um angleArray 10318个角的球粒子的测试
%总和偏差一般小于 2.5%
%0.0000    0.0026    0.0026       NaN    
%0.0000    0.0273    0.0251       NaN       
% NaN      0.0259    0.0281       0.0000       
% NaN      0.0032    0.0032       0.0000

randMaxNum = 100 ;

deltaArray = zeros(1,16) ;
deltaMatrix = [] ;
for randIter = 1 : randMaxNum
    
    t_phiDeg = rand(1) * 360 ;
    
    mieTempMuller_rotated = DDA_MullerT_PhiFillUp( t_phiDeg , mieTempMuller ) ;
    for eleIter = 3 : 18 
        deltaArray( eleIter - 2 ) = MatrixCompare( mieTempMuller(:,eleIter) , mieTempMuller_rotated(:,eleIter) ) ;
    end
    
    deltaMatrix = [ deltaMatrix ; deltaArray ] ;
end
for iter = 1 : 16
meanDeltaMatrix( iter ) = mean( deltaMatrix( : , iter ) ) ;
end