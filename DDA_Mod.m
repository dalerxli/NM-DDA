function result = DDA_Mod( modedNum , modNum )
%DDA_MOD
% DDA_Mod(1,3) = 1
% DDA_Mod(2,3) = 2
% DDA_Mod(3,3) = 3
% DDA_Mod(4,3) = 1
% ...

if mod( modedNum , modNum ) == 0
    result = modNum ;
else result = mod( modedNum , modNum ) ;
end

end

