function pastTime = TimePast( startTime , currentTime )
%TIMEPAST 此函数计量已经经过了多少时间
%startTime , currentTime 都应该是clock函数的返回值
sec = 1 ; 
min = 60 ;
hou = 3600 ;
day = 3600 * 24 ;

if currentTime(3) < startTime(3)
    pastTime = day + ( currentTime(4) - startTime(4) ) * hou ...
        + ( currentTime(5) - startTime(5) ) * min ...
        + ( currentTime(6) - startTime(6) ) * sec  ;

else
    pastTime = ( currentTime(3) - startTime(3) ) * day + ( currentTime(4) - startTime(4) ) * hou ...
        + ( currentTime(5) - startTime(5) ) * min ...
        + ( currentTime(6) - startTime(6) ) * sec  ;
end


end

