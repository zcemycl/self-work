% image processing_3
% ==================================================== %
% Erosion
% ==================================================== %
function [E] = erosion(S,x,y)
% E = S(-)SK
SK =[0 1 0; 1 1 1; 0 1 0];
% S: image; xxy result image

% Dimension
Dimen = size(S);
% Intermediate 3x3  image
I2 = zeros(3,3); % for erosion

% Erosion S(-)SF = E(5x5)
E = zeros(Dimen(1)+2,Dimen(2)+2);
dimen = size(E);
E(dimen(1)-Dimen(1):dimen(1)-1,dimen(2)-Dimen(2):dimen(2)-1) = S;
EE = zeros(dimen);
for i = 1:Dimen(1)
    for j = 1:Dimen(2)
        I2 = E(i:i+2,j:j+2);
        if I2 == double(or(I2,SK))
            EE(i+1,j+1) = 1;
        end
    end 
end
clear E
if x <= dimen(1) && y <= dimen(2)
    x = (x-1)/2; y = (y-1)/2;
    xbar = (dimen(1)+1)/2 ; ybar = (dimen(2)+1)/2 ;
    E = EE(xbar - x:xbar + x,ybar - y:ybar + y);
else
    E = zeros(x,y);
    E(x-dimen(1):dimen(1)+1 , y-dimen(2):dimen(2)+1) = EE;
    
end
end
