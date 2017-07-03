angleArray = [];
for iterTheta = 0 : 180
    angleArray = [angleArray ; deg2rad( iterTheta ) 0 iterTheta 0] ;
end

save('angleArray.mat','angleArray') ;