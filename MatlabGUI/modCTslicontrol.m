function modCTslicontrol(x)
f = figure('Name','2-D Slicing control','NumberTitle','off','Position',[5,5,900,720]); 
% Open the target file
fid=fopen(x);
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

details = dir(x);
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

% Sample profile zy,zx,yx
xc = 270; assignin('base','xc',270);
yc = 270; assignin('base','yc',270);
zc = 50; assignin('base','zc',50);
% -------------------------------------------------%
% z-y plane
% -------------------------------------------------%
xlabel = uicontrol('Parent',f,'Style','text',...
            'string','X-Axis slicing',...
            'Position',[30,40,80,15],'tag','xempty');
sliderx = uicontrol('Parent',f,'Style','slider',...
            'Value',xc,'min',1,'max',512,'SliderStep',[5/511 5/511],...
            'Position',[120,40,150,15],'Callback',{@sliderx_callback},'tag','xempty');
xvalue = uicontrol('Parent',f,'Style','edit',...
            'string',num2str(xc),...
            'Position',[280,40,25,15],'tag','xempty','Callback',{@xvalue_callback,sliderx});
plotimagex(xc);
    function plotimagex(xc)
        axx = axes('Units','Pixels','Position',[30,90,400,280]);
        imagesc(squeeze(V(xc,:,:))'); 
        daspect(1./[vox_dims(2) vox_dims(3) 1]);
        axis xy;
        colormap(gray);
    end
    function sliderx_callback(hObject,eventdata)
        xc = round(get(hObject,'Value'));
        assignin('base','xc',xc);
        set(xvalue,'string',xc);
        plotimagex(xc);
    end
    function xvalue_callback(source,eventdata, sliderx)
        xc = str2double(get(xvalue,'String'));
        assignin('base','xc',xc);
        set(sliderx,'Value',xc);
        plotimagex(xc);
    end 
% -------------------------------------------------%
% z-x plane
% -------------------------------------------------%
ylabel = uicontrol('Parent',f,'Style','text',...
            'string','Y-Axis slicing',...
            'Position',[460,40,80,15],'tag','xempty');
slidery = uicontrol('Parent',f,'Style','slider',...
            'Value',yc,'min',1,'max',512,'SliderStep',[5/511 5/511],...
            'Position',[550,40,150,15],'Callback',{@slidery_callback},'tag','xempty');
yvalue = uicontrol('Parent',f,'Style','edit',...
            'string',num2str(yc),...
            'Position',[710,40,25,15],'tag','xempty',...
            'Callback',{@yvalue_callback,slidery});
plotimagey(yc);
    function plotimagey(yc)
        axy = axes('Units','Pixels','Position',[460,90,400,280]);
        imagesc(squeeze(V(:,yc,:))'); 
        daspect(1./[vox_dims(1) vox_dims(3) 1]);
        axis xy;
        colormap(gray);
    end
    function slidery_callback(hObject,eventdata)
        yc = get(hObject,'Value');
        assignin('base','yc',yc);
        set(yvalue,'string',yc);
        plotimagey(yc);
    end
    function yvalue_callback(source,eventdata, slidery)
        yc = str2double(get(yvalue,'String'));
        assignin('base','yc',yc);
        set(slidery,'Value',yc);
        plotimagey(yc);
    end
% -------------------------------------------------%
% y-x plane
% -------------------------------------------------%
zlabel = uicontrol('Parent',f,'Style','text',...
            'string','Z-Axis slicing',...
            'Position',[210,360,80,15],'tag','xempty');
sliderz = uicontrol('Parent',f,'Style','slider',...
            'Value',zc,'min',1,'max',107,'SliderStep',[5/106 5/106],...
            'Position',[300,360,150,15],'Callback',{@sliderz_callback},'tag','xempty');
zvalue = uicontrol('Parent',f,'Style','edit',...
            'string',num2str(zc),...
            'Position',[460,360,25,15],'tag','xempty',...
            'Callback',{@zvalue_callback,sliderz});
plotimagez(zc);
    function plotimagez(zc)
        axz = axes('Units','Pixels','Position',[210,410,400,280]);
        imagesc(squeeze(V(:,:,zc))'); 
        daspect(1./[vox_dims(1) vox_dims(2) 1]);
        axis xy;
        colormap(gray);
    end
    function sliderz_callback(hObject,eventdata)
        zc = get(hObject,'Value');
        assignin('base','zc',zc);
        set(zvalue,'string',zc);
        plotimagez(zc);
    end
    function zvalue_callback(source, eventdata, sliderz)
        zc = str2double(get(zvalue,'String'));
        assignin('base','zc',zc);
        set(sliderz,'Value',zc);
        plotimagez(zc);
    end
% -------------------------------------------------%

end