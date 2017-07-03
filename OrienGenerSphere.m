total = 0 ;
inter = 4;
for i = 0 :inter:180
    total = total + 2*pi*sin(deg2rad(i)) / deg2rad(inter) ;
end

figure ;
x=[];
y=[];
z=[];
angleArray = [] ;

for i = 0 : inter : 180
        
    if i == 0 || i == 180
        phiMaxNum = 1 ;
        phiArray = [ 0 ] ;
    else
    
    
        phiMaxNum = round ( 2*pi*sin(deg2rad(i)) / deg2rad(inter) ) ;
        phiArray = [0 : phiMaxNum - 1] * 360 / phiMaxNum ;
    end
    
    for j = 1 : phiMaxNum
        newx = sin(deg2rad(i)) * sin(deg2rad(phiArray(j))) ;
        newy = -cos(deg2rad(i)) ;
        newz = sin(deg2rad(i)) * cos(deg2rad(phiArray(j))) ;
        
        x = [ x , newx];
        y = [ y , newy];
        z = [ z , newz] ;
        
        angleArray = [angleArray ; deg2rad(i) deg2rad(phiArray(j)) i phiArray(j) ] ;
    end
end

plot3(x,y,z,'.')

%输出的angleArray的格式 第1、2列是theta和phi的角度（弧度制），第3、4列是theta和phi的角度（角度制）

save('angleArray.mat','angleArray') ;