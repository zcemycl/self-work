% image_processing_6
% finding edges in squeeze()

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

% display H(original),T(binary),E(erosion of binary)
% ,T-E (difference map)
H = squeeze(V(270,:,:));
sp1 = subplot(1,4,1);
imagesc(H');
daspect(1./[107,512,1]);

T = threshold(H);
sp2 = subplot(1,4,2);
imagesc(T');
daspect(1./[107,512,1]);

E = erosion(T,512,107);
sp3 = subplot(1,4,3);
imagesc(E');
daspect(1./[107,512,1]);

ER = T - E;
sp4 = subplot(1,4,4);
imagesc(ER');
daspect(1./[107,512,1]);