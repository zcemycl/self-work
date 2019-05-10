% image processing_4
% ==================================================== %
% Dilation
% ==================================================== %
function [D] = dilation(S,x,y)
% D = S(+)SK
SK = [0 1 0; 1 1 1; 0 1 0];
% S: image; xxy result

% Dimension
Dimen = size(S);
% Intermediate 3x3 image
I1 = zeros(3,3);

% Dilation S(+)SK = D(5x5)
D = zeros(Dimen(1)+2,Dimen(2)+2);
dimen = size(D);
D(dimen(1)-Dimen(1):dimen(1)-1,dimen(2)-Dimen(2):dimen(2)-1) = S;

for i = 1:Dimen(1)
    for j = 1:Dimen(2)
        if S(i,j) == 1
            I1 = D(i:i+2,j:j+2);
            I1 = double(or(I1,SK));
            D(i:i+2,j:j+2) = I1;
        end
    end 
end
EE = D;
clear D
if x <= dimen(1) && y <= dimen(2)
    x = (x-1)/2; y = (y-1)/2;
    xbar = (dimen(1)+1)/2 ; ybar = (dimen(2)+1)/2 ;
    D = EE(xbar - x:xbar + x,ybar - y:ybar + y);
else
    D = zeros(x,y);
    D(x-dimen(1):dimen(1)+1 , y-dimen(2):dimen(2)+1) = EE;
    

end