% edge, corner, blob
% Open the target file
fid = fopen('CT_vol.nii');
% move to a position 40 bytes from the beginning of the file.
fseek(fid,40,'bof') ;
% Read the 8 image dimensions
img_dims=fread(fid,8,'short');
% move to a position 70 bytes from the beginning of the file.
fseek(fid,70,'bof') ;
% Read the datatype code
datatype_code=fread(fid,1,'short');
switch datatype_code
    case 2
        datatype_string = 'uint8';
    case 4
        datatype_string = 'int16';
    case 8
        datatype_string = 'int32';
    case 16
        datatype_string = 'float';
    case 64
        datatype_string = 'double';
    otherwise
        error('This datatype is not supported');
end

details = dir('CT_vol.nii');
% move to the start of the image data.
fseek(fid,352,'bof');
% Read the image into a 1D vector
V= fread(fid,prod(img_dims),datatype_string);
% Reshape the image so it has the correct dimensions
V= reshape(V,img_dims(2:4)');
fseek(fid,76,'bof');
vox_dims = fread(fid,8,'float');
vox_dims = vox_dims(2:4)';
assignin('base','V',V);
fclose(fid);

n = 2; m = 4;

colormap gray;

H = squeeze(V(270,:,:));
sp1 = subplot(n,m,1);
imagesc(H');
daspect(1./[107,512,1]);


filter_x = [-1,0,1];
filtered_x = convolution(H, filter_x, 0);
sp2 = subplot(n,m,2);
imagesc(filtered_x');
daspect(1./[107,512,1]);

filter_y = [-1;0;1];
filtered_y = convolution(H, filter_y, 1);
sp3 = subplot(n,m,3);
imagesc(filtered_y');
daspect(1./[107,512,1]);


filtered_2 = filtered_x.^2 + filtered_y.^2;
sp4 = subplot(n,m,4);
imagesc(filtered_2');
daspect(1./[107,512,1]);

filtered_xy = filtered_x.*filtered_y;
sp5 = subplot(n,m,5);
imagesc(filtered_xy');
daspect(1./[107,512,1]);

filter_xx = [1,-2,1];
filtered_xx = convolution(H,filter_xx, 0);
sp6 = subplot(n,m,6);
imagesc(filtered_xx');
daspect(1./[107,512,1]);

function filtered = convolution(image, filter, direction)
    dim = size(image);
    filtered = zeros(dim);
    if direction == 1
        Intermediate_y = zeros(3,1);
        I = zeros(dim(1)+2, dim(2));
        Dim = size(I);
        I(Dim(1)-dim(1):Dim(1)-1,1:dim(2)) = image;
        for i = 1:dim(1)
            for j = 1:dim(2)
                Intermediate_y = I(i:i+2, j);
                filtered(i,j) = dot(Intermediate_y, filter);                
            end
        end
    else
        Intermediate_x = zeros(1,3);
        I = zeros(dim(1), dim(2)+2);
        Dim = size(I);
        I(1:dim(1),Dim(2)-dim(2):Dim(2)-1 ) = image;
        for i = 1:dim(1)
            for j = 1:dim(2)
                Intermediate_x = I(i, j:j+2);
                filtered(i,j) = dot(Intermediate_x, filter);                
            end
        end
    end
end
