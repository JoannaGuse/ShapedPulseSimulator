function [ ] = Export_wavfileGUI_Standalone( varargin )
%This GUI exports pulse data to AWG compatible .wav files according to
%specified formatting inputs.
%
%examples:
%Export_wavfileGUI_Standalone(); % will result in a loading data window.
%Export_wavfileGUI_Standalone(pulse); %will load pulse 
%pulse must have the fields: pulse.Bx (I), pulse.By (Q), pulse.signal (Y),
%pulse.t

addpath(genpath(pwd))

hfig = figure('Visible','off','color','w','Position',[0,0, 800, 600]);
set(hfig, 'Name', 'Export Wav file')
movegui(hfig,'center')

%% _______________make panels_________________
hinputp = uipanel('Title','Input','Fontsize', 12,...
    'BackgroundColor','white',...
    'Position',[.05 .75,.3 .2]);

hOutputList = uipanel('Title','Output Selection','Fontsize', 12,...
    'BackgroundColor','white',...
    'Position',[.05 .5 .3 .3]);

hformat = uipanel('Title','Format Output','Fontsize', 12,...
    'BackgroundColor','white',...
    'Position',[.35 .5 .6 .45]);

hplots = uipanel('Title','Plot inputs/outputs','Fontsize', 12,...
    'BackgroundColor','white',...
    'Position',[.05 .03 .9 .48]);

%% ______________ make Input panel________________________

hImport=uicontrol('parent', hinputp,'Style','pushbutton',...
    'String','Import pulse',...
    'units', 'normalized','position',[.1,.4,.7,.4],...
    'Callback',@himport_Callback);

%% ______________ make Output panel________________________
hI = uicontrol('parent',hOutputList,'Style','Radio','String','I','backgroundColor', 'white',...
    'Value', 1, 'units', 'normalized', 'pos',[.1, .8, .8, .1],'tag', 'I', 'Callback', @hI_Callback);
hQ = uicontrol('parent',hOutputList,'Style','Radio','String','Q','backgroundColor', 'white',...
    'Value', 1,'units', 'normalized', 'pos',[.1, .6, .8, .1],'tag', 'Q', 'Callback', @hQ_Callback);
hY = uicontrol('parent',hOutputList,'Style','Radio','String','Y','backgroundColor', 'white',...
    'Value', 0,'units', 'normalized', 'pos',[.1, .4, .8, .1],'tag', 'Y', 'Callback', @hY_Callback);

%_________________make export button____________________
hExport=uicontrol('parent', hOutputList,'Style','pushbutton',...
    'String','Export Wav Files',...
    'units', 'normalized','position',[.1,.1,.7,.2],...
    'Callback',@hExport_Callback);


%% ______________________Format Panel_______________
[ txt_Scale, edit_Scale] = Make_TxtEditpair( hformat,'Scale','Amp Scale (%)', [.02, .8, .3, .1], 'white', '1',0.5,@Scale_callback, 'on');
[ txt_SRate, edit_SRate] = Make_TxtEditpair( hformat,'SampleRate','AWG SampleRate (Sa/s)', [.02, .7, .3, .1], 'white', '2.3e9',0.5,@SampleRate_callback, 'on');
[ txt_bits, edit_bits] = Make_TxtEditpair( hformat,' bits','AWG bits', [.02, .6, .3, .1], 'white', '14',0.5,@bits_callback, 'on');
[ txt_minLength, edit_minLength] = Make_TxtEditpair( hformat,'Min. Length','AWG Min Length(CLCK)', [.02, .5, .3, .1], 'white', '384',0.5,@minLength_callback, 'on');
[ txt_Lengthinc, edit_Lengthinc] = Make_TxtEditpair( hformat,'Length incr.','AWG Length inc(CLCK)', [.02, .4, .3, .1], 'white', '32',0.5,@Lengthinc_callback, 'on');
[ txt_maxLength, edit_maxLength] = Make_TxtEditpair( hformat,'Max Length','AWG Max Length(CLCK)', [.02, .3, .3, .1], 'white', '16e6',0.5,@maxLength_callback, 'on');

hLeadZero = uicontrol('parent',hformat,'Style','Radio','String','Add Leading Zero','backgroundColor', 'white',...
    'Value', 1,'units', 'normalized', 'pos',[.02, .2, .3, .1],'tag', 'Lead0', 'Callback', @LeadZero_Callback);

Make_Txt(hformat,'zero','Some oscilloscopes idle on the first point of the waveform while waiting for a trigger',...
    [.02, .05, .45, .15], 'white', 'on', 'Left', 8);

[ txt_nalias, edit_nalias] = Make_TxtEditpair( hformat,'nalias','anti-alias filter order', [.55, .6, .25, .1], 'white', '0',0.5,@nalias_callback, 'on');
[ txt_filter, edit_filter] = Make_TxtEditpair( hformat,'filter','Kaiser filter beta', [.55, .5, .25, .1], 'white', '0',0.5,@filter_callback, 'on');


hs=[txt_nalias, edit_nalias, txt_filter, edit_filter];
[ htxt, hEndedness] = Make_Txtpopuppair(hformat, 'Endedness','Endedness',{'little-endian','big-endian'}, [.55,.82,.2,.1], 'white',0.5, @Endedness_callback, 'on');
[ htxt, hresample] = Make_Txtpopuppair(hformat, 'Resample','Resample',{'ON','OFF'}, [.55,.72,.2,.1], 'white',0.5, @(hresample,eventdata)resample_callback(hresample, eventdata, hs), 'on');


%__________set defaults___________________
data=guidata(hfig);
resample_callback(hresample, [], hs);
Endedness_callback(hEndedness, []);
minLength_callback(edit_minLength, []);
maxLength_callback(edit_maxLength, []);
Lengthinc_callback(edit_Lengthinc, []);
Scale_callback(edit_Scale, []);
SampleRate_callback(edit_SRate, []);
bits_callback(edit_bits, []);
nalias_callback(edit_nalias, []);
filter_callback(edit_filter, []);
guidata(hfig, data);

%% ____________________load data_____________________
data=guidata(hfig);
if isempty(varargin)
      himport_Callback(hfig, [])
%     pulse.t=linspace(0, 230e-9, 1024)';
%     pulse.Bx=zeros(size(pulse.t));
%     pulse.By=zeros(size(pulse.t));
%     pulse.signal=zeros(size(pulse.t));
%     data.Bxscale=pulse.Bx;
%     data.Byscale=pulse.By;
%     data.Signal=pulse.signal;
else
    pulse=varargin{1};
    
    if max(abs((pulse.Bx)))==0
        data.Bxscale=zeros(length(pulse.t),1);
    else
        data.Bxscale=pulse.Bx./max(abs(pulse.Bx));
    end
    
    if  max(abs((pulse.By)))==0
        data.Byscale=zeros(length(pulse.t),1);
    else
        data.Byscale=pulse.By./max(abs(pulse.By));
    end
    
    if max(abs((pulse.signal)))==0
        data.Signal=zeros(length(pulse.t),1);
    else
        data.Signal=abs(pulse.signal)./max(abs(pulse.signal));
        
    end
    
end

guidata(hfig, data)
LeadZero_Callback(hLeadZero, []);
SelectivePlot;


%% _________________make preview button____________________
hpreview=uicontrol('parent', hformat,'Style','pushbutton',...
    'String','Preview Wav File Output',...
    'units', 'normalized','position',[.55,.05,.4,.2],...
    'Callback',@hpreview_Callback);

%% _________________make plot panel______________________________
SelectivePlot;

%% general functions
    function SelectivePlot(varargin)
        
        Ion=get(hI, 'Value');
        Qon=get(hQ, 'Value');
        Son=get(hY, 'Value');
        
        data=guidata(hfig);
        try
            delete(data.haxes)
        end
        
        data.haxes(1) = axes('Parent', hplots,'Position',[.1,.2,.25,.7], 'Visible', 'on');
        plot(data.tt/1e-9,data.Bx, 'Linewidth', 3);
        hold on
        title('I ')
        xlabel('Time (ns)')
        ylim([-1, 1])
        
        data.haxes(2) = axes('Parent', hplots,'Position',[.4,.2,.25,.7],  'Visible', 'on');
        plot(data.tt/1e-9,data.By, 'Linewidth', 3);
        hold on
        title('Q')
        xlabel('Time (ns)')
        ylim([-1, 1])
        
        data.haxes(3) = axes('Parent', hplots,'Position',[.7,.2,.25,.7],  'Visible', 'off');
        plot(data.tt/1e-9,data.sig, 'Linewidth', 3);
        hold on
        title('Y')
        xlabel('Time (ns)')
        ylim([-1, 1])
        
        guidata(hfig, data);
        
    end

    function [ htxt, hedit] = Make_TxtEditpair( parent, name,string, normpos, color, defaultVal,boxscale, callbackfunc, vis)
        
        htxt  = uicontrol('Parent', parent,'tag', name, 'Style','text','String',string,...
            'Visible', vis,'BackgroundColor',color,'units', 'normalized', 'position', normpos);
        pos=get(htxt, 'position');
        hedit = uicontrol('Parent', parent,'tag', name,'Style','edit','units', 'normalized',...
            'Position',[pos(1)+pos(3), pos(2),0.13, boxscale*pos(4)],'String', defaultVal, 'BackgroundColor','white',...
            'Visible', vis,'Callback', callbackfunc);
        clear pos
        align([htxt, hedit],'VerticalAlignment','Top');
    end

    function [ htxt, hpopup] = Make_Txtpopuppair( parent, name,string,popup, normpos, color,boxscale, callbackfunc, vis)
        
        htxt  = uicontrol('Parent', parent,'tag', name, 'Style','text','String',string,...
            'Visible', vis,'BackgroundColor',color,'units', 'normalized', 'position', normpos);
        pos=get(htxt, 'position');
        
        hpopup = uicontrol('Parent', parent,'tag', name,'Style','popupmenu','units', 'normalized',...
            'Position',[pos(1)+pos(3), pos(2),0.2, boxscale*pos(4)],'String', popup, 'BackgroundColor','white',...
            'Visible', vis,'Callback', callbackfunc);
        clear pos
        align([htxt, hpopup],'VerticalAlignment','Top');
        
    end

    function [ x,Val,I ] = hasfield( s, fname, varargin)
        % X = HASFIELD(S,FNAME) checks if the struct S has the fieldname FNAME. If
        % S is not a struct, X will be false.
        %
        % [X,L] = HASFIELD(S,FNAME,LEVEL) checks LEVEL number of levels into possible
        % nested structres for a field name. L returns the first level that the
        % field name was found at.
        
        L = 0;% level
        I = 0;% level iteratios
        x = 0;% default
        Val=0;
        
        if(~isstruct(s))
            warning('input is not a struct');
            return
        end
        
        if(length(varargin) > 0)
            if(~isnumeric(varargin{1}))
                error('LEVEL must be a numeric variable');
            end
            L = varargin{1};
        end
        recurseon({s});
        function recurseon(T)
            I = I+1;
            %return if max iterations
            if(L > 0 && I > L)
                I = I-1;
                return
            end
            
            %check for any matching fieldnames
            NEWT = {};
            for i=1:length(T)
                
                t = T{i};
                
                %detection
                x = isfield(t,fname);
                
                if(x)
                    Val=getfield(t, fname);
                    return
                end
                
                %breadth first elaboration
                N = fieldnames(t);
                for j = 1:length(N)
                    if(isstruct(t.(N{j})))
                        NEWT = vertcat(NEWT,t.(N{j}));
                    end
                end
                
            end
            
            %if theres actually anything to recurse on
            if(~isempty(NEWT))
                recurseon(NEWT);
            end
            
        end
    end

%________________________________________________________________
    function [ htxt] = Make_Txt( parent, tagname,string, normpos, color, vis, justification, fontsize)
        htxt  = uicontrol('Parent', parent,'tag', tagname, 'Style','text','String',string,...
            'Visible', vis,'BackgroundColor',color,'units', 'normalized','Fontsize', fontsize, ...
            'position', normpos, 'HorizontalAlignment', justification);
    end

%% Format callback functions

    function LeadZero_Callback(Source, eventdata)
        data=guidata(Source);
        pad=get(hLeadZero, 'Value');
        if pad==1
            data.tt=[pulse.t; (pulse.t(end)+pulse.t(2)-pulse.t(1))];
            data.Bx=[0; data.Bxscale];
            data.By=[0; data.Byscale];
            data.sig=[0; data.Signal];
        else
            data.tt=pulse.t;
            data.Bx=data.Bxscale;
            data.By=data.Byscale;
            data.sig=data.Signal;
        end
        
        
        guidata(Source, data);
        SelectivePlot;
        
    end

%_______________________________
    function Endedness_callback(Source, ~)
        data=guidata(Source);
        str = get(Source, 'String');
        val = get(Source,'Value');
        switch str{val}
            case 'little-endian'
                data.machinefmt='ieee-le';
            case 'big-endian'
                data.machinefmt='ieee-be';
        end
        guidata(Source, data);
    end

%_______________________________
    function resample_callback(Source, eventdata, h)
        data=guidata(Source);
        str = get(Source, 'String');
        val = get(Source,'Value');
        data.resample=str{val};
        
        if val==2
            set(h(1), 'Visible', 'off')
            set(h(2), 'Visible', 'off')
            set(h(3), 'Visible', 'off')
            set(h(4), 'Visible', 'off')
        else
            set(h(1), 'Visible', 'on')
            set(h(2), 'Visible', 'on')
            set(h(3), 'Visible', 'on')
            set(h(4), 'Visible', 'on')
        end
        guidata(Source, data);
    end

%_______________________________
    function nalias_callback(Source, ~)
        data=guidata(Source);
        nalias = str2double(get(Source,'String'));
        if isnan(nalias)
            errordlg('input numeric number')
            uicontrol(Source)
            return
        else
            data.nalias=nalias;
            guidata(Source,data);
        end
    end

%_______________________________
    function filter_callback(Source, ~)
        data=guidata(Source);
        filter= str2double(get(Source,'String'));
        if isnan( filter)
            errordlg('input numeric number)')
            uicontrol(Source)
            return
        else
            data.filter= filter;
            guidata(Source,data);
        end
    end

%_______________________________
    function SampleRate_callback(Source, ~)
        data=guidata(Source);
        SampleRate = str2double(get(Source,'String'));
        if isnan(SampleRate)
            errordlg('input numeric number (Sa/s)')
            uicontrol(Source)
            return
        else
            data.SampleRate=SampleRate;
            guidata(Source,data);
        end
    end

%_______________________________
    function Scale_callback(Source, ~)
        data=guidata(Source);
        ml = str2double(get(Source,'String'));
        if isnan(ml)
            errordlg('input numeric number')
            uicontrol(Source)
            return
        else
            data.Vscale=ml;
            guidata(Source,data);
        end
    end

%_______________________________
    function minLength_callback(Source, ~)
        data=guidata(Source);
        ml = str2double(get(Source,'String'));
        if isnan(ml)
            errordlg('input numeric number')
            uicontrol(Source)
            return
        else
            data.MinLength=ml;
            guidata(Source,data);
        end
    end

%_______________________________
    function maxLength_callback(Source, ~)
        data=guidata(Source);
        ml = str2double(get(Source,'String'));
        if isnan(ml)
            errordlg('input numeric number')
            uicontrol(Source)
            return
        else
            data.MaxLength=ml;
            guidata(Source,data);
        end
    end

%_______________________________
    function Lengthinc_callback(Source, ~)
        data=guidata(Source);
        ml = str2double(get(Source,'String'));
        if isnan(ml)
            errordlg('input numeric number')
            uicontrol(Source)
            return
        else
            data.Lengthinc=ml;
            guidata(Source,data);
        end
    end

%_______________________________
    function bits_callback(Source, ~)
        data=guidata(Source);
        bits= str2double(get(Source,'String'));
        if isnan(bits)
            errordlg('input numeric number')
            uicontrol(Source)
            return
        else
            data.bits=bits;
            guidata(Source,data);
        end
    end

%% Export callback functions
%_______________________________
    function hpreview_Callback(Source, eventdata)
        LeadZero_Callback(hLeadZero, []);
        data=guidata(Source);
        outputs=[get(hI, 'Value'),get(hQ, 'Value'),get(hY, 'Value')];
        sig=[data.Bx, data.By, data.sig];
        
        t1=data.tt;
        dt=t1(2)-t1(1);
        data.Ypad=[];
        
        for kk=1:3
            try
                delete(data.htemp(kk));
            end
        end
        
        for jj=1:3
            if outputs(jj)==1
                y1=sig(:,jj);
                
                %% resample data to make sure pulse length is correct given the SampleRate
                if strcmp(data.resample, 'ON')==1
                    N1=length(y1);
                    Tp1=t1(end);
                    FS=data.SampleRate;
                    NN=round(Tp1*FS);
                    ydata=resample(y1, NN, N1, data.nalias, data.filter);
                else
                    ydata=y1;
                end
                
                %% check AWG specifications and pad if necessary
                n=length(ydata);
                maxamp=2^data.bits/2;
                
                pad=data.Lengthinc-rem(n,data.Lengthinc);
                padarray=zeros(1, pad);
                %               padarray=ones(1, pad)*ydata(end);  %pad data with the last value of y, not with zeros
                output=[ydata', padarray(1:end)];
                if  length(output)<data.MinLength
                    output=[output, zeros(1, data.MinLength-length(output))];
                elseif rem(output, data.Lengthinc)==1
                    output=output;
                elseif length(output)>data.MaxLength
                    output=output(1:data.MaxLength);
                else
                end
                
                %make sure to output in units of AWG DAC resolution
                if max(abs(output))==0
                    Ypad=output;
                else
                    Ypad=output./max(abs(output))*data.Vscale*maxamp;
                end
                
                % replace maximum output with next best value to avoid overflow
                [~, ii]=find(round(Ypad)==maxamp);
                Ypad(ii)=maxamp-1;
                [~, ss]=find(round(Ypad)==-maxamp);
                Ypad(ss)=-maxamp+1;
                
                
                %% Write the waveform
                fid = fopen('temp.wav','w'); % Open .wav file
                bitout=nextpow2(data.bits);
                fwrite(fid,Ypad,['bit', num2str(2^bitout)]); % Write 16 bit binary data
                fclose(fid); %
                
                data.Ypad(:,jj)=Ypad;
                %read waveform
                fid= fopen('temp.wav','r'); % Open .wav file
                [yy, ~]=fread(fid,['bit', num2str(2^bitout)], data.machinefmt); % Write 16 bit binary data
                fclose(fid); % Close .wav file
                tt=(0: length(yy)-1)*1/data.SampleRate; %make t as if sampled by AWG
                
                %% plot
                if max(abs(yy))==0
                    yplot=yy;
                else
                    yplot=yy./max(abs(yy));
                end
                data.htemp(jj)=plot(data.haxes(jj), tt/1e-9, yplot, 'r.-');
                
                set(data.haxes(jj),'Yaxislocation','left','YLim',[-1, 1],'Ytick', [-1, 0, 1],...
                    'YTickLabel', {['-2^',num2str(data.bits-1)],'0',['2^',num2str(data.bits-1)]});
                delete('temp.wav')
                
            end
        end
        guidata(Source, data);
    end
%______________________
    function hExport_Callback(Source, eventdata)
        hpreview_Callback(hpreview)
        data=guidata(Source);
        outputs=[get(hI, 'Value'),get(hQ, 'Value'),get(hY, 'Value')];
        
        try
            [file,path] = uiputfile('.wav','Save Simulation As');
            filenames={[file(1:end-4), '_I.wav'], [file(1:end-4), '_Q.wav'], [file(1:end-4), '_Y.wav']};
        end
        
        for jj=1:3
            if outputs(jj)==1
                fid = fopen([path, filenames{jj}],'w'); % Open .wav file
                bitout=nextpow2(data.bits);
                fwrite(fid,data.Ypad(:,jj),['bit', num2str(2^bitout)]); % Write 16 bit binary data
                fclose(fid); %
            end
        end
        
    end

%_______________________________
    function himport_Callback(Source, eventdata)
        data=guidata(Source);
        
        [filename1,filepath1]=uigetfile({'*.*','All Files'},'Select Data File');
        rawdata1=load([filepath1, filename1]);
        
        [exists,pulse, ~] = hasfield(rawdata1, 'pulse', 10);
        if exists==1
            data.pulse=pulse;
        else
            pname=inputdlg('enter variabe pulse is stored in:');
           pulse=getfield(rawdata1, pname{1});
           data.pulse=pulse;
        end

        
        if max(abs((pulse.Bx)))==0
            data.Bxscale=zeros(length(pulse.t),1);
        else
            data.Bxscale=pulse.Bx./max(abs(pulse.Bx));
        end
        
        if max(abs((pulse.By)))==0
            data.Byscale=zeros(length(pulse.t),1);
        else
            data.Byscale=pulse.By./max(abs(pulse.By));
        end
        
        if max(abs((pulse.signal)))==0
            data.Signal=zeros(length(pulse.t),1);
        else
            data.Signal=abs(pulse.signal)./max(abs(pulse.signal));
        end
        
        %         try
        %         delete(data.haxes)
        %         end
        %         data.haxes=SelectivePlot;
        %
        guidata(Source, data)
        
    end
%_______________________________
    function hI_Callback(Source, eventdata)
        
    end
%_______________________________
    function hQ_Callback(Source, eventdata)
        
    end
%_______________________________
    function hY_Callback(Source, eventdata)
        
    end

%% visible on
author = uicontrol('Parent', hfig, 'Style','text','String','Author: Joanna A. Guse',...
    'Visible', 'on','BackgroundColor','white','units', 'normalized', 'position', [.7, 0.01, .3 .02]);
set(hfig, 'Visible', 'on')
end

