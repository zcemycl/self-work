% image processing_5
% Thresholding to binary images
function [T] = threshold(S)
    dim = size(S);
    for i = 1:dim(1)
        for j = 1:dim(2)
            if S(i,j) > 0 
                T(i,j) = 1;
            else 
                T(i,j) = 0;
            end
        end
    end
end