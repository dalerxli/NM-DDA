%这个脚本用于生成球散射的mie解 mieTempMuller ，用于比对DDA的结果
%load
load('angleArray.mat') ;

%设置参数
lambda = 0.532 ;
radius = 0.5 ;
refrel = 1.33 ;

nang = 91 ;
%initial
x = 2 * pi * lambda / radius ;


%get S11...
[S11,S12,S33,S34]=mie(x,refrel,nang) ;

%get spheremuller
mullerNum = size(S11,2) ;
sphereMuller = zeros(4,4,mullerNum) ;

for i = 1 : mullerNum
    sphereMuller(1,1,i) = S11(i);
    sphereMuller(1,2,i) = S12(i);
    sphereMuller(2,1,i) = S12(i);
    sphereMuller(2,2,i) = S11(i);
    sphereMuller(3,3,i) = S33(i);
    sphereMuller(4,4,i) = S33(i);
    sphereMuller(3,4,i) = S34(i);
    sphereMuller(4,3,i) = -S34(i);
    
end

%生成带phi角的muller矩阵
for iter = 1 : size( angleArray , 1 )
    
    Ltheta = angleArray( iter , 3 ) ;
    Lphi = angleArray( iter , 4 ) ;
    
    %注意此处Lphi的正负号
    temp = sphereMuller(:,:,Ltheta+1) * [1 0 0 0; 0 cos( 2 * deg2rad(Lphi) ) -sin(2*deg2rad(Lphi)) 0; 0 sin(2*deg2rad(Lphi)) cos(2*deg2rad(Lphi)) 0; 0 0 0 1] ;
    
    mieTempMuller(iter,1) = angleArray(iter , 3) ;
    mieTempMuller(iter,2) = angleArray(iter , 4) ;
    
    mieTempMuller(iter,3) = temp( 1 , 1 ) ;
    mieTempMuller(iter,4) = temp( 1 , 2 ) ;
    mieTempMuller(iter,5) = temp( 1 , 3 ) ;
    mieTempMuller(iter,6) = temp( 1 , 4 ) ;
    
    mieTempMuller(iter,7) = temp( 2 , 1 ) ;
    mieTempMuller(iter,8) = temp( 2 , 2 ) ;
    mieTempMuller(iter,9) = temp( 2 , 3 ) ;
    mieTempMuller(iter,10) = temp( 2 , 4 ) ;
    
    mieTempMuller(iter,11) = temp( 3 , 1 ) ;
    mieTempMuller(iter,12) = temp( 3 , 2 ) ;
    mieTempMuller(iter,13) = temp( 3 , 3 ) ;
    mieTempMuller(iter,14) = temp( 3 , 4 ) ;
    
    mieTempMuller(iter,15) = temp( 4 , 1 ) ;
    mieTempMuller(iter,16) = temp( 4 , 2 ) ;
    mieTempMuller(iter,17) = temp( 4 , 3 ) ;
    mieTempMuller(iter,18) = temp( 4 , 4 ) ;

end

save('mieTempMuller.mat','mieTempMuller') ;
