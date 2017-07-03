%�˽ű����ڲ���DDA_MullerT_PhiFillUp�����Ƿ���ȷ
%�˽ű��ķ�������mieTempMuller��Ϊ���룬����������תʱ�ǲ����ģ��ʶ���ʱɢ������Ӧ�ò��ᷢ���仯

%20170514��m=1.33 aeff=0.5 um lambda = 0.532um angleArray 10318���ǵ������ӵĲ���
%�ܺ�ƫ��һ��С�� 2.5%
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