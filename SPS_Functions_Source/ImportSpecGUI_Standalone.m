function []=ImportSpecGUI_Standalone(varargin)
hfig = figure('Visible','off','color','w','Position',[0,0, 400, 200]);
set(hfig, 'Name', 'Import Spectrum')
movegui(hfig,'center')

author = uicontrol('Parent', hfig, 'Style','text','String','Joanna A. Guse',...
    'Visible', 'on','BackgroundColor','white','units', 'normalized', 'position', [.74, 0., .3 .09]);
set(hfig, 'Visible', 'on')

specData=guidata(hfig);
specData.B0=3440e-4;
specData.fuw=9.7e9;
guidata(hfig, specData);
xx=[];
yy=[];

%% Select format
hBrowse= uicontrol('parent', hfig,'style', 'pushbutton',...
    'units','normalized', 'position',[.1, .8, .25, .15],...
    'string','Browse', 'callback', @BrowseInput);
txt_file = uicontrol('Parent', hfig, 'Style','text','String','filename',...
    'Visible', 'on','BackgroundColor',[.8, .8, .8],'units', 'normalized',...
'position', [.35, .8, .6 .15]);

%% parameter panel
hImportParams = uipanel('parent', hfig,'Title','Import parameters',...
    'BackgroundColor','white','Position',[.1, .3, .85, .5]);

[ txt_xcol, edit_xcol] = Make_TxtEditpair(hImportParams, 'xcol','x column index:', [.1, .7, .3, .2], 'white', '1',1, @xcol_callback, 'ON');
[ txt_ycol, edit_ycol] = Make_TxtEditpair(hImportParams, 'ycol','y column index :', [.1, .4, .3, .2], 'white', '2',1, @ycol_callback, 'ON');
[ txt_startrow, edit_startrow] = Make_TxtEditpair(hImportParams, 'startrow','start row index:', [.1, .1, .3, .2], 'white', '2',1, @startrow_callback, 'ON');
[ txt_xvar, edit_xvar] = Make_TxtEditpair(hImportParams, 'xvar','x variable name:', [.1, .8, .3, .2], 'white', 'B',1, @xvar_callback, 'OFF');
txt_xvar_descr = uicontrol('Parent',hImportParams, 'Style','text','String','The variable x should correspond to B (gauss)',...
    'Visible', 'on','BackgroundColor','white','units', 'normalized','tag', 'xdescr',...
'position',[.052, .55, .8, .2], 'Visible', 'OFF');
[ txt_yvar, edit_yvar] = Make_TxtEditpair(hImportParams, 'yvar','y variable name:', [.1, .25, .3, .2], 'white', 'Intensity',1, @yvar_callback, 'OFF');
 txt_yvar_descr = uicontrol('Parent',hImportParams, 'Style','text','String','The variable y should be the spectrum',...
    'Visible', 'on','BackgroundColor','white','units', 'normalized','tag', 'ydescr',...
'position',[.01, .0, .8, .2], 'Visible', 'OFF');            
guidata(hfig, specData);


hImportButton= uicontrol('parent',hfig,'style', 'pushbutton', ...
    'units','normalized', 'position',[.1, .1, .3, .15],...
    'string','Import Spectrum', 'callback', @ImportSpec_callback);
Note = uicontrol('Parent', hfig, 'Style','text','String','Note: figure will close automatically when import is done. Do not close figure',...
    'Visible', 'on','BackgroundColor','white','units', 'normalized', 'position', [.45, 0.1, .5 .15]);


%% callback funtions
   function [ htxt, hedit] = Make_TxtEditpair( parent, name,string, normpos, color, defaultVal,boxscale, callbackfunc, vis)
        
        htxt  = uicontrol('Parent', parent,'tag', name, 'Style','text','String',string,...
            'Visible', vis,'BackgroundColor',color,'units', 'normalized', 'position', normpos);
        pos=get(htxt, 'position');
        pos2=[pos(1)+pos(3), pos(2),0.4, pos(4)];
        hedit = uicontrol('Parent', parent,'tag', name,'Style','edit','units', 'normalized',...
            'Position',[pos(1)+pos(3), pos(2),0.4, boxscale*pos(4)],'String', defaultVal, 'BackgroundColor','white',...
            'Visible', vis,'Callback', callbackfunc);
        clear pos
        align([htxt, hedit],'VerticalAlignment','Top');
   end

%____________________________________________________
    function BrowseInput(Source, eventdata)
        specData=guidata(Source);
                [fname, path] = uigetfile( ...
            {  '*.DTA; *.mat; *.csv; *.xlsx; *.txt;*.ASCII','supported Files (*.DTA, *.mat, *.csv, *.xlsx, *.txt,*.ASCII)';...
            '*.*',  'All Files (*.*)'}, ...
            'Select Spertrum File');

        [pathstr,name,ext] = fileparts(fname);
        specData.fname=fname;
        specData.ImportExtenstion=ext;
        specData.filename=[path, fname];
        set(txt_file, 'String', fname);
        guidata(Source, specData);
        
       
        switch specData.ImportExtenstion
            case {'.DTA'}
                showhide({'xvar', 'yvar','xcol', 'ycol', 'startrow'},{'B0'})
          
            case {'.txt', '.ASCII'}
                set(edit_xcol, 'String', '2')
                set(edit_ycol, 'String', '3')
                showhide({'xvar', 'yvar'},{'xcol', 'ycol', 'startrow'})
            case '.csv'
                set(edit_xcol, 'String', '1')
                set(edit_ycol, 'String', '2')
                showhide( {'xvar', 'yvar'},{'xcol', 'ycol','startrow'})
            case '.xlsx'
                set(edit_xcol, 'String', '1')
                set(edit_ycol, 'String', '2')
                showhide({'xvar', 'yvar','startrow'},{'xcol', 'ycol'})
            case '.mat'
                showhide({'xcol', 'ycol','startrow'}, {'xvar', 'yvar', 'xdescr', 'ydescr'})
            otherwise
            showhide({'xvar', 'yvar','xcol', 'ycol', 'startrow'},{'B0'})        
        end
  end

%____________________________________________________
    function ImportSpec_callback(Source, ~)
        specData=guidata(Source);
        switch specData.ImportExtenstion
            case {'.DTA'}
                [xx, yy]= LoadData_Elexys(specData.filename);                
            case {'.txt', '.ASCII'}
                startrow_callback(edit_startrow,[]);
                xcol_callback(edit_xcol,[]);
                ycol_callback(edit_ycol,[]);
                [xx, yy]=LoadData_txt(specData.filename);
            case '.csv'
                startrow_callback(edit_startrow,[]);
                xcol_callback(edit_xcol,[]);
                ycol_callback(edit_ycol,[]);
                [xx, yy]=LoadData_csv(specData.filename);
            case '.xlsx'
                xcol_callback(edit_xcol,[]);
                ycol_callback(edit_ycol,[]);
                [xx, yy]=LoadData_Excel(specData.filename);
            case '.mat'
                xvar_callback(edit_xvar,[]);
                yvar_callback(edit_yvar,[]);
                [xx, yy]=LoadData_mat(specData.filename);
            otherwise
        end

        specData.x=xx;
        specData.y=yy;
        guidata(Source, specData);
        assignin('base', 'specData',specData)
        
        h2=figure;
        plot( specData.x, specData.y);
        pause(1)
        close(hfig)
    end


%____________________________________________________
    function [xx, yy]=LoadData_txt( filename)
        
        startRow = specData.startrow;
        delim='';
        formatSpec = '%21f%21f%21f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delim, 'WhiteSpace', '', 'ReturnOnError', false);
        fclose(fileID);
        x = dataArray{:, specData.xcol};
        y_real= dataArray{:, specData.ycol};
        xx=x*1e-4;
        yy=y_real./max(y_real);
        
    end
%____________________________________________________
    function [xx, yy]=LoadData_csv(filename)
        startRow = specData.startrow;
        delim=',';
        formatSpec = '%10f%21f%21f%21f%[^\n\r]';
        fileID = fopen(filename,'r');
        textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delim, 'WhiteSpace', '', 'ReturnOnError', false);
        fclose(fileID);
        x = dataArray{:, specData.xcol};
        y_real= dataArray{:, specData.ycol};
        
        xx=x*1e-4;
        yy=y_real./max(y_real);
        
    end
%____________________________________________________
    function [xx, yy]=LoadData_Excel(filename)
        [numbers, strings] = xlsread(filename);
        xx= numbers(:,specData.xcol)*1e-4;
        y=real(numbers(:,specData.ycol));
        yy=y./max(y);
        
    end
%____________________________________________________
    function [xx, yy]= LoadData_mat(filename)
        s=load(filename);  
        xx=getfield(s, specData.xvar)*1e-4;
        y=getfield(s, specData.yvar);       
        y=real(y);
        yy=y./max(y);       
    end
%____________________________________________________
    function [xx, y_real]= LoadData_Elexys( file)
        
        xx=[];
        y_real=[];
        M=1;
        
        file=strrep(file,'.DTA','');
        file=strrep(file,'.dta','');
        fid = fopen( [file '.DTA'], 'r','ieee-be.l64') ;% spc
        %         fid = fopen( [pathname file '.DTA'], 'r','ieee-be.l64') ;% spc
        %check existance of file
        if(fid == -1)
            return;
        end
        %load DTA file and read into DSC file.
        [Z,N]= fread(fid,inf,'float64');
        fclose(fid);
        fid = fopen([file '.DSC'], 'r');% par
        [S,~]= fscanf(fid,'%c');
        fclose(fid);
        
        % %Parse string S, the delimiter is ascii 10, this is the token delimiter
        Delim = setstr(10);
        VB=[];
        while(length(S) > 0 )
            [token,S] = strtok(S,Delim);
            VB = str2mat(VB,token);
        end
        
        GST = get_vb(VB,'XMIN');
        GSI = get_vb(VB,'XWID');
        N = get_vb(VB,'XPTS');
        dx = GSI/(N-1);
        x=GST + dx*(0:(N-1));
        dim2=get_vb(VB,'YTYP');
        if ~strncmp(dim2,'NODATA',3),
            ymin = get_vb(VB,'YMIN'); ywid = get_vb(VB,'YWID'); M = get_vb(VB,'YPTS');
            dy = ywid/(M-1);
            y=ymin + dy*[0:(M-1)];
        end;
        
        if length(Z)==2*N*M,
            z0=Z;
            Z=zeros(N,1);
            for k=1:N*M,
                Z(k)=z0(2*k-1)+1i*z0(2*k);
            end;
        end;
        Z=reshape(Z,N,M);
        Z=permute(Z,[2,1]);
        
        xx=x*1e-4;
        y_real=real(Z)./max(real(Z));
        specData=guidata(hfig);
        specData.fuw=double(get_vb(VB,'MWFQ'));
        specData.B0=double(get_vb(VB,'A1CT'));
        guidata(hfig, specData);
        
    end
%_______________________________________________
    function [ value , vstr] = get_vb(VB,mnem)      
        % Elexys file read
        % $Date: 2008/04/01 19:29:01 $
        % $Revision: 3.19 $   NOT Release Version        
        [m,n] = size(VB);
        value = [];
        vstr =[];
        k =1;
        while( (k<m) && ~ strcmp( strtok(VB(k,:)),mnem ) )
            k=k+1;
        end
        if (k<m) [t,vstr] = strtok(VB(k,:));
            value = str2num(vstr);
        end
        if(isempty(value)),
            vstr0=vstr;
            vstr=[];
            for k=1:length(vstr0)
                if ~isspace(vstr0(k))
                    vstr=[vstr vstr0(k)];
                end
            end
            value = vstr;
        end
       
    end
%_______________________________________________

%% Edit callbakcs
%______________________________________________________________
    function startrow_callback(Source, eventdata)
        specData=guidata(Source);
        
        startrow = str2double(get(Source,'String'));
        if isnan(startrow)
            errordlg('input numeric number')
            uicontrol(Source)
            return
        else
            specData.startrow=startrow;
            guidata(Source,specData);
        end
    end

%______________________________________________________________
    function xcol_callback(Source, eventdata)
        specData=guidata(Source);
        
        xcol = str2double(get(Source,'String'));
        if isnan(xcol)
            errordlg('input numeric number')
            uicontrol(Source)
            return
        else
            specData.xcol=xcol;
            guidata(Source,specData);
        end
    end

%______________________________________________________________
    function ycol_callback(Source, eventdata)
        specData=guidata(Source);
        
        ycol = str2double(get(Source,'String'));
        if isnan(ycol)
            errordlg('input numeric number')
            uicontrol(Source)
            return
        else
            specData.ycol=ycol;
            guidata(Source,specData);
        end
    end

%______________________________________________________________
    function yvar_callback(Source, eventdata)
        specData=guidata(Source);        
        yvar =get(Source,'String');
        specData.yvar=yvar;
        guidata(Source,specData);       
    end
%______________________________________________________________
    function xvar_callback(Source, eventdata)
        specData=guidata(Source);        
        xvar =get(Source,'String');
        specData.xvar=xvar;
        guidata(Source,specData);       
    end

%___________________________________________________
    function [  ] = showhide( hide, show )
        for jj = 1:length(hide)
            temp= findobj('Tag', hide{jj});
            set(temp, 'Visible', 'off');
        end
        
        for ii =1:length(show)
            temp= findobj('Tag', show{ii});
            set(temp, 'Visible', 'on');
        end
        
    end
end
