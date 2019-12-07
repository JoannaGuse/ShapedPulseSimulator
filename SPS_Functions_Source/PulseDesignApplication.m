
function [  ] = PulseDesignApplication( varargin )
%GUI for designing, simulating and exporting a single pulse. 
%Syntax:
%PulseDesignApplication();
close all
clear
clc

maincol=[207,224,250]/255;
maincol='white';
hfig = figure('Visible','off','color',maincol,'Position',[0,0, 1100, 900]);
set(hfig, 'Name', 'Square_Setup')
movegui(hfig,'center')
set(hfig, 'Resize', 'on')


author = uicontrol('Parent', hfig, 'Style','text','String','Author: Joanna A. Guse',...
    'Visible', 'off','BackgroundColor',maincol,'units', 'normalized', 'position', [.7, 0.98, .4 .02]);
% set Matlab Version specific functions
MATLABversion=version('-release');
MatVersion=str2double(MATLABversion(1:4));
    
%% set up default data
data.sim.M0=[0,0,1];
data.sim.df=1;%MHz
data.sim.fmax=100; %MHz

data.pulsename='square';
data.concatenation='None';
data.composite=@(pulse) pulse;
data.pdef=[];
data.coord='cartesian';
data.out.Mx=[];
data.out.My=[];
data.out.Mz=[];
data.out.theta=[];
data.out.phi=[];
data.out.f=[];

data.Spec.x=[];
data.Spec.y=[];
data.Spec.fuw=9.7e9;
data.Spec.B0=3440;

guidata(hfig,data)

%% make sim panel

% maincol=[.9, .9, .9];
% % hsimcolor=maincol;
% hsimcolor=maincol;
% hImportcolor=hsimcolor;
% importcol=hsimcolor;
% simoutcolor=hsimcolor;

maincol=[.9, .9, .9];
hsimcolor=maincol;
hImportcolor=maincol;
importcol=maincol;
simoutcolor=maincol;

hsim = uipanel('Title','Simulation Panel','Fontsize', 14,...
    'BackgroundColor',hsimcolor,...
    'Position',[.505 .01 .49 .84]);

%% make Import Panel
hImportPanel = uipanel('Parent', hfig,'Title','Spectrum Import Panel','Fontsize', 14', ...
    'BackgroundColor', hImportcolor,'Position',[.505, .84, .49, .15]);
hImportButton= uicontrol('parent', hImportPanel,'style', 'pushbutton', ...
    'units','normalized', 'position',[.01, .6, .3, .3],...
    'string','Import Spectrum','Fontsize', 12, 'callback', @hImportSpec_Callback);
txt_filename = uicontrol('Parent', hImportPanel, 'Style','text',...
    'String','','Visible', 'on','BackgroundColor',hImportcolor,...
    'units', 'normalized', 'position', [.31, .5, .65, .45]);

[ txt_fuw, edit_fuw] = Make_TxtEditpair2(hImportPanel, 'fuw','fuw (GHz):', [.01, .1, .2, .2], .15,hImportcolor, num2str(data.Spec.fuw/1e9),1, @fuw_Callback, 'ON');
[ txt_B0, edit_B0] = Make_TxtEditpair2(hImportPanel, 'B0','B0 (Gauss):', [.5, .1, .2, .2], .15,hImportcolor, num2str(data.Spec.B0),1, @B0_Callback, 'ON');


hShowSpec=uicontrol('parent',hImportPanel,'Style','radio',...
    'String','ShowSpec', 'units', 'normalized',...
    'position',[.01, .4, .3, .15], 'Visible', 'off',...
    'Callback',@ShowSpec_Callback);

%% make sim input_____________________________
pp=get(hImportPanel,'Position');
hsiminput = uipanel('Parent',hsim,'Title','Simulation inputs',...
    'BackgroundColor',importcol,'Position',[.0,.91, 1, .08]);
[txt_M0,edit_M0] = Make_TxtEditpair2( hsiminput, 'M0','M0:', [0, 0, .07, .85],.2, importcol, '[0,0,1]',0.8,@M0_Callback, 'on' );
[txt_df, edit_df] = Make_TxtEditpair2( hsiminput, 'df','Resolution (MHz):', [.3, .0, .2, .5], .1,importcol, '1',0.8,@df_Callback, 'on' );
[txt_fmax,edit_fmax] = Make_TxtEditpair2( hsiminput, 'fmax','Freq. Sweep (MHz):', [.3, .5, .2, .5],.1, importcol, '100',0.8,@fmax_Callback, 'on' );

[txt_T1, edit_T1] = Make_TxtEditpair2( hsiminput, 'T1','T1 (s):', [.65, .5, .2, .5], .1,importcol, 'inf',0.8,@T1_Callback, 'on' );
[txt_T2,edit_T2] = Make_TxtEditpair2( hsiminput, 'T2','T2 (s):', [.65, .0, .2, .5],.1, importcol, 'inf',0.8,@T2_Callback, 'on' );

M0_Callback(edit_M0, []);
df_Callback(edit_df,[]);
fmax_Callback(edit_fmax, []);
T1_Callback(edit_T1,[]);
T2_Callback(edit_T2,[]);
%%  make Sim Output selection________________________
pp=get(hsiminput,'Position');
hsimOutput = uipanel('Parent',hsim,'Title','Simulation Outputs',...
    'BackgroundColor',simoutcolor,'Position',[.0, 0, 1, .9]);

%__________ make plot axes______________
fa=12;
fl=14;
data.haxes(1) = axes('Parent', hsimOutput,'Position',[.2,.52,.7,.38]);
ylabel('\theta (\pi)', 'fontsize', fl)
ylim([0, 1])
set(gca, 'fontsize', fa)
data.haxes(2) = axes('Parent', hsimOutput,'Position',[.2,.1,.7,.38]);
ylabel('\phi (\pi)', 'fontsize', fl)
xlabel('Freq (MHz)', 'fontsize', fl)
ylim([-1, 1])
set(gca, 'fontsize', fa)
guidata(hfig, data);


hCoordButtons = uibuttongroup('Parent', hsimOutput,'visible','off',...
    'Position',[.01, 0.94, 0.5, .05], 'backgroundColor', simoutcolor, 'Tag', 'CoordSel');
hout1 = uicontrol('parent',hCoordButtons,'Style','Radio','String','cartesian','backgroundColor',simoutcolor,...
    'units', 'normalized', 'pos',[0.1, .1, .4, .8],'tag', '1','HandleVisibility','off');
hout2 = uicontrol('parent',hCoordButtons,'Style','Radio','String','polar','backgroundColor',simoutcolor,...
    'units', 'normalized', 'pos',[.5, .1, .4, .8],'tag', '2', 'HandleVisibility','off');
set(hCoordButtons,'SelectionChangeFcn', @Coordinate_selection);
set(hCoordButtons,'SelectedObject',hout1);
set(hCoordButtons,'Visible','on');

hNewPlot=uicontrol('parent',hsimOutput,'Style','pushbutton','String','Export Plot','Fontsize', 12,...
    'units', 'normalized','position',[.6, 0.94, 0.3, .05],...
    'Callback',@hNewPlot_Callback);


%%____________________________________________________________________________________________________________-
%% make pulse panel
% hpcolor=[.9, .8, .8];
hpcolor=maincol;
inputscol=maincol;
incol=maincol;
hpoutputcol=maincol;

hp = uipanel('Title','Pulse Panel','FontSize',14,...
    'BackgroundColor',hpcolor,...
    'Position',[.01 .01 .49 .98]);

hpinput= uipanel('Parent',hp,'Title','Pulse Inputs','Fontsize', 14,...
    'BackgroundColor',inputscol,...
    'Position',[.02 .081 0.5,.69]);

hInfoButton= uicontrol('parent', hp,'style', 'pushbutton', ...
    'units','normalized', 'position',[.8, .98, .2, .03],...
    'string','Software Info','Fontsize', 12, 'callback', @hInfo_Callback);


% hpoutput = uipanel('Parent',hp,'Title','Pulse Outputs','Fontsize', 14,...
%     'BackgroundColor',hpoutputcol,...
%     'Position',[.54 .081 .45 .69]);
hpoutput = uipanel('Parent',hp,'Title','Pulse Outputs','Fontsize', 14,...
    'BackgroundColor',hpoutputcol,...
    'Position',[.54 .081 .45 .54]);


hCompSel = uibuttongroup('Parent', hp,'visible','off','Title','Pulse Concatenation','Fontsize', 14,...
    'Position',[.54 .62 .45 .15], 'backgroundColor', hpoutputcol, 'Tag', 'CompositeSel');

houta = uicontrol('parent',hCompSel,'Style','Radio','String','None','backgroundColor',hpoutputcol,...
    'units', 'normalized', 'pos',[0.1, .8, .9, .2],'tag', 'None','HandleVisibility','off');

houtb = uicontrol('parent',hCompSel,'Style','Radio','String','PhaseTime','backgroundColor',hpoutputcol,...
    'units', 'normalized', 'pos',[.1, .45, .9, .2],'tag', 'PhaseTime', 'HandleVisibility','off');

houtc = uicontrol('parent',hCompSel,'Style','Radio','String','BIR4','backgroundColor',hpoutputcol,...
    'units', 'normalized', 'pos',[.1, .1, .9, .2],'tag', 'BIR4', 'HandleVisibility','off');

set(hCompSel,'SelectionChangeFcn', @Composite_selection);
set(hCompSel,'SelectedObject',houta);
set(hCompSel,'Visible','on');
[txt_BIRphase,edit_BIRphase] = Make_TxtEditpair(hCompSel, 'BIRphase','BIRPhase (pi):',[.45, .05, .32, .2], inputscol, '1',0.7,@BIRphase_Callback, 'off' );


hi1 = uipanel('Parent',hpinput,'title', 'Pulse independent parameters','fontweight', 'bold',...
    'BackgroundColor',incol,...
    'Position',[.0 .8 1 .2]);

hi2 = uipanel('Parent',hpinput,'title', 'Pulse shaping parameters','fontweight', 'bold',...
    'BackgroundColor',incol,...
    'Position',[.0 .7 1 .1]);

hi3 = uipanel('Parent',hpinput,'title', 'Pulse design constraints','fontweight', 'bold',...
    'BackgroundColor',incol,...
    'Position',[.0 .3 1 .4]);

hi4 = uipanel('Parent',hpinput,'title', 'FM pulse design parameters','fontweight', 'bold',...
    'BackgroundColor',incol,...
    'Position',[.0 .0 1 .3]);
%% input panel text
height1=0.4;
base=0.0;
dy=0.3;

[txt_npts,edit_npts] = Make_TxtEditpair( hi1, 'npts',' Number of points [npts]:',[.05, base+2*dy, .7, height1], incol, '1024',0.4,@npts_Callback, 'on' );
[txt_f0,edit_f0] = Make_TxtEditpair( hi1, 'f0','Center frequnecy [f_center] (MHz):', [.05,  base+dy, .7, height1], incol, '0',0.4,@f0_Callback, 'on' );
[txt_phi0,edit_phi0] = Make_TxtEditpair( hi1, 'phi0','Initial phase [phi0] (pi):', [.05, base, .7, height1], incol, '0',0.4,@phi0_Callback, 'on' );
% Make_Txt( hpinput,'divline','______________________________',[.0, .72, 1, .05] , inputscol, 'ON', 'Left', 10);

height1=0.5;
base=0.0;
dy=0.4;
%  Make_Txt( hpinput,'divline','______________________________',[.0, .2, 1, .01] , inputscol, 'ON', 'Left', 10);
[txt_n,edit_n] = Make_TxtEditpair( hi2, 'n','n:', [.05, .5, .7, height1], incol, '1',0.7,@n_Callback, 'on' );
[txt_sigma,edit_sigma] = Make_TxtEditpair( hi2, 'sigma','sigma:', [.05,  0, .7, height1], incol,'0.23',0.7,@sigma_Callback, 'on' );

height1=0.2;
base=0.01;
dy=0.18;

[txt_theta,edit_theta] = Make_TxtEditpair( hi3, 'theta','Rotation angle [theta] (pi):', [.05, base+4*dy, .7, height1], incol, '1',0.4,@theta_Callback, 'on' );
[txt_opBW,edit_opBW] = Make_TxtEditpair( hi3, 'opBW','Operation bandwidth [opBW] (MHz):', [.05, base+3*dy, .7, height1], incol, '50',0.4,@opBW_Callback, 'on' );

[txt_B1,edit_B1] = Make_TxtEditpair( hi3, 'B1','RF power [B1_max] (G):', [.05,  base+2*dy, .7, height1], incol, '5.6',0.4,@B1_Callback, 'on' );
[txt_Tp, edit_Tp]= Make_TxtEditpair( hi3, 'Tp','Pulse length [Tp] (ns):', [.05, base+dy, .7, height1], incol, '32',0.4,@Tp_Callback, 'on' );
[txt_wmax,edit_wmax] = Make_TxtEditpair( hi3, 'fmax','Frequency sweep limit [fmax] (MHz):', [.05,  base, .7, height1], incol, '40',0.4,@wmax_Callback, 'on' );


height1=0.3;
base=0.1;
dy=0.25;

[txt_R,edit_R] = Make_TxtEditpair( hi4, 'R','Time-bandwidth product [R]=Tp*fmax:', [.05,  base+2*dy, .7, height1], incol, '2',0.4,@R_Callback, 'on' );
[txt_K,edit_K] = Make_TxtEditpair( hi4, 'K','K=(gamma B1)/fmax:', [.05,  base+dy, .7, height1], incol, '0.5',0.4,@K_Callback, 'on' );
[txt_node,edit_node] = Make_TxtEditpair( hi4, 'node','chirp solution number [node]:', [.05,  base, .7, height1], incol, '2',0.4,@node_Callback, 'on' );


EvaluateAllInputs(hfig)
data=guidata(hfig);
data.AllInputHandles=[edit_npts, edit_phi0, edit_f0, edit_B1,edit_theta, edit_opBW, edit_Tp, edit_wmax, edit_n, edit_sigma, edit_R, edit_K, edit_node, edit_BIRphase];
guidata(hfig, data);

%% outputs panel text
%________ output text______________
out_theta = Make_Txt(hpoutput, 'out_theta','rotation angle theta (pi)=0', [.05, .8, .9, .2], hpoutputcol, 'ON', 'Left',10);
out_opBW = Make_Txt(hpoutput, 'out_opBW','operation BW (MHz)=0', [.05, .7, .9, .2], hpoutputcol, 'ON', 'Left',10);

out_B1 = Make_Txt(hpoutput, 'out_B1','B1max (Gauss)=0', [.05, .55, .9, .2], hpoutputcol, 'ON', 'Left',10);
out_Tp = Make_Txt(hpoutput, 'out_Tp','Tp pulse length (ns)=0', [.05, .45, .9, .2], hpoutputcol, 'ON', 'Left',10);
out_fmax = Make_Txt(hpoutput, 'out_fmax','Freq. sweep limit (MHz)=0', [.05, .35, .9, .2], hpoutputcol, 'ON', 'Left',10);

out_BWP = Make_Txt(hpoutput, 'out_BWP','Time-opBW product=0', [.05, .2, .9, .2], hpoutputcol, 'ON', 'Left',10);
out_R = Make_Txt(hpoutput, 'out_R','Time-BW product=NA', [.05, .1, .9, .2], hpoutputcol, 'ON', 'Left',10);
out_Q = Make_Txt(hpoutput, 'out_Q','Adiabaticity Q=NA', [.05, .005, .9, .2], hpoutputcol, 'ON', 'Left',10);

 
hPlotEnv= uicontrol('parent', hpoutput,'style', 'pushbutton', ...
    'units','normalized', 'position',[.01, .0, .98, .06],...
    'string','Plot Pulse Modulation Functions','Fontsize', 12, 'callback', @PlotEnv_Callback);
%% pulse Selection
PopH = uicontextmenu('Parent',hfig);
hPulseSelect = uicontrol('parent', hp, 'Style', 'PushButton', 'units', 'normalized', 'Position', [.05, .93, .25, .04], ...
    'String', 'Select Pulse','Fontsize', 14, 'UIContextMenu', PopH, 'Callback', {@showMyPopup, PopH});

pname = uicontrol('parent', hp, 'Style', 'text', 'units', 'normalized','Position',[.3, .93, .2, .03],...
    'String', 'No pulse selected','Fontsize', 12, 'backgroundcolor', hpcolor, 'UIContextMenu', PopH,'Callback',@Callback_pname);

    % Create submenu items
    m1 = uimenu('Parent',PopH,'Label','AM');
    m2 = uimenu('Parent',PopH,'Label','FM');
    
   pa = uimenu('Parent',m2,'Label','chirp','Callback',@Selectpulse);
   pb = uimenu('Parent',m2,'Label','HSn','Callback',@Selectpulse);
   pc = uimenu('Parent',m2,'Label','WURST','Callback',@Selectpulse);
  

   ps = uimenu('Parent',m1,'Label','square','Callback',@Selectpulse);
   p2 = uimenu('Parent',m1,'Label','gauss','Callback',@Selectpulse);
   p3 = uimenu('Parent',m1,'Label','EBURP','Callback',@Selectpulse);
   p4 = uimenu('Parent',m1,'Label','IBURP','Callback',@Selectpulse);
   p5 = uimenu('Parent',m1,'Label','UBURP','Callback',@Selectpulse);
   p6 = uimenu('Parent',m1,'Label','REBURP','Callback',@Selectpulse);
   p7 = uimenu('Parent',m1,'Label','G3','Callback',@Selectpulse);
   p8 = uimenu('Parent',m1,'Label','G4','Callback',@Selectpulse);
   p9 = uimenu('Parent',m1,'Label','Q3','Callback',@Selectpulse);
   p10 = uimenu('Parent',m1,'Label','Q5','Callback',@Selectpulse);
   p11 = uimenu('Parent',m1,'Label','Hermite90','Callback',@Selectpulse);
   p12 = uimenu('Parent',m1,'Label','Hermite180','Callback',@Selectpulse);

function showMyPopup(ButtonH, EventData, PopH)
% Pos = get(ButtonH, 'Position');
% set(PopH, 'Position', [Pos(1)+70, Pos(2)+855], 'Visible', 'on');
set(ButtonH, 'units', 'pixels');
Pos = get(ButtonH, 'Position');
set(ButtonH, 'units', 'normalized');
set(PopH, 'Position', [Pos(1)+Pos(3)/3, Pos(2)+Pos(4)/3], 'Visible', 'on');

end

function Selectpulse(Source, eventdata)
        data=guidata(Source);        
        val=get(Source, 'Label');
        data.pulsename=val;
        guidata(Source, data);     
        Callback_pname(pname, []);
        PulseSelect_Callback(hfig, [])
end

function Callback_pname(Source, eventdata)
        data=guidata(Source);        
        set(Source, 'string', data.pulsename)        
        guidata(Source, data);        
    end

%%  pulse Method Generation Panel
   
pulseDescr = Make_Txt(hp,'pulseDescr','Amplitude modulated pulse for universal operations.',...
    [.01, .82, .9, .1], hpcolor, 'on','Left', 8);

hregime = uicontrol('Parent', hp, 'Style','popupmenu',...
    'Visible','on','String',{'Adiabatic', 'Rapid Passage'},'Fontsize', 12,...
    'units', 'normalized','Position',[.65, .865, .3, .1],...
    'Tag', 'SweepRate','Callback',@hRegime_Callback);

htxt = Make_Txt(hp,'pulse','Select Pulse Generation Method:',[.01, .77, .55, .1], hpcolor, 'on','Left', 12);
hMethodSelect = uicontrol('Parent', hp, 'Style','popupmenu',...
    'String',{'Fixed operation BW', 'Fixed B1', 'Fixed Tp','Fully Specified'},'Fontsize', 12,...
    'units', 'normalized','Position',[.55, .77, .4, .1],...
    'Tag', 'SelectMethod','Callback',@MethodSelect_Callback);
methodDescr = Make_Txt(hp,'methodDescr','Pulse is optimized to perform a uniform theta operation within the specified operation bandwidth.',...
    [.01, .77, .9, .05], hpcolor, 'on','Left', 8);

% set(hMethodSelect, 'Value', 0);
PulseSelect_Callback(hPulseSelect, []);
MethodSelect_Callback(hMethodSelect, []);
hRegime_Callback(hregime, []);


hExport=uicontrol('parent', hp,'Style','pushbutton',...
    'String','Save to .mat','Fontsize', 12,...
    'units', 'normalized','position',[.7, .01, .3, 0.06],...
    'Callback',@hExport_Callback);

hExportWav=uicontrol('parent', hp,'Style','pushbutton',...
    'String','Export .wav file','Fontsize', 12,...
    'units', 'normalized','position',[.38, .01, .3, 0.06],...
    'Callback',@hExportWav_Callback);
hAnalyze=uicontrol('parent', hp,'Style','pushbutton',...
    'String','Simulate pulse','Fontsize', 14,...
    'units', 'normalized','position',[.0, .01, .36, 0.06],...
    'Callback',@hAnalyze_Callback);

%% %______________________________________________________________________________

%% Pulse Simulation Callbacks
%______________________________________________________
    function EvaluateAllInputs(Source)
        data=guidata(Source);
        npts_Callback(edit_npts, []);
        phi0_Callback(edit_phi0, []);
        f0_Callback(edit_f0, []);
        
        B1_Callback(edit_B0, []);
        Tp_Callback(edit_Tp, []);
        wmax_Callback(edit_wmax, []);
        opBW_Callback(edit_opBW, []);
        theta_Callback(edit_theta, []);
        R_Callback(edit_R, []);
        
        n_Callback(edit_n, []);
        sigma_Callback(edit_sigma, []);
        BIRphase_Callback(edit_BIRphase, []);
        
        guidata(Source, data);
    end

%____________________________________________________
    function PulseSelect_Callback(Source, evendtata)
       
        data=guidata(Source);
        Callback_pname(pname, []); 
        clear p1 pdef
        p1.name=data.pulsename;
           
        p1.npts=str2num(get(edit_npts, 'String'));
        p1=feval(['Def_', p1.name], p1);
        set(pulseDescr, 'String', p1.descr);
   
        pdef.name=p1.name;
        pdef.class=p1.class;
        pdef.DynamicVars=p1.DynamicVars;
        
        pdef.descr=p1.descr;
             
        data=rmfield(data, 'pdef');
        data=rmfield(data, 'p');
        data.pdef=pdef;
        data.p=pdef;
                         
       guidata(Source,data);
        
        set(edit_theta, 'Enable', 'on');
        set([edit_n, edit_sigma], 'visible', 'off');
        
  
        for jj=1:length(pdef.DynamicVars)
            field=pdef.DynamicVars{jj};
            val=p1.(field);
            a=eval(['edit_',field]);
            set(a,'String', num2str(val), 'visible', 'on');
        end
                
        val=p1.theta;
        if val==pi;
            set(edit_theta,'String','1', 'visible', 'on');
        elseif val==pi/2;
            set(edit_theta,'String','0.5', 'visible', 'on');
        else
        end
        
        if any(ismember(pdef.DynamicVars,'theta'))==1 && strcmp(data.concatenation, 'BIR4')==1
            set(edit_theta,'String','1', 'visible', 'on');
            set(edit_theta, 'Enable', 'off');
        end
        
        if any(ismember(pdef.DynamicVars,'theta'))==0
            set(edit_theta, 'Enable', 'off');
        end
           
        if strcmp(pdef.class, 'FM')==1
            set(hregime, 'visible', 'on'); 
            hRegime_Callback(hregime, []);
        else
            set(hregime, 'visible', 'off');
        end
           
        hRegime_Callback(hregime, [])
       
    end

% %__________________________________________________________
    function MethodSelect_Callback(Source, ~)
        data=guidata(Source);
        str = get(Source, 'String');
        val = get(Source,'Value');
        set([edit_R, edit_K, edit_node], 'Visible', 'off');
        
        switch str{val}
            case 'Fixed operation BW'
                method='Fixed_opBW';
                Descr='Pulse is optimized to perform a uniform theta operation within the specified operation bandwidth.';
                opBW_Callback(edit_opBW, []);
                theta_Callback(edit_theta, []);
                
                set([edit_opBW, edit_theta], 'Visible', 'on');
                set([edit_B1, edit_Tp, edit_wmax], 'Visible', 'off');
                
            case 'Fixed B1'
                method='Fixed_B1';
                Descr='Pulse is optimized to perform a desired theta operation for a given pulse power (B1max).';
                B1_Callback(edit_B1, []);
                theta_Callback(edit_theta, []);
                set([edit_B1, edit_theta], 'Visible', 'on');
                set([edit_opBW, edit_Tp, edit_wmax], 'Visible', 'off');
                
            case 'Fixed Tp'
                method='Fixed_Tp';
                Descr='Pulse is optimized to perform a desired theta operation for a given pulse length (Tp).';
                Tp_Callback(edit_Tp, []);
                theta_Callback(edit_theta, []);
                set([edit_Tp,edit_R, edit_theta], 'Visible', 'on');
                set([edit_B1,edit_opBW, edit_wmax], 'Visible', 'off');
                
            case 'Fully Specified'
                method='None';
                Descr='Pulse is fully specified with pulse power (B1max), length (Tp) and frequence sweep (fmax). The pulse rotation angle is not optimized.';
                B1_Callback(edit_B1, []);
                Tp_Callback(edit_Tp, []);
                set([edit_opBW, edit_theta], 'Visible', 'off');
                set([edit_B1, edit_Tp, edit_wmax], 'Visible', 'on');
                if strcmp(data.p.class, 'AM')==1
                    set( edit_wmax, 'Visible', 'off');
                end
            otherwise
                warndlg('undefined method');
        end
        
        set(methodDescr, 'String', Descr);
        data.p.method=method;
        
        guidata(Source,data);
        Select_DesignParams(Source)
    end

%________________________________________________
    function Select_DesignParams(Source)
        data=guidata(Source);
        try
            data.p=rmfield(data.p, 'K');
        end
        try
            data.p=rmfield(data.p, 'R');
        end
        try
            data.p=rmfield(data.p, 'node');
        end
        
        if strcmp(data.p.class, 'AM')==1
            set([ edit_R,  edit_node, edit_K], 'visible', 'off');
            try
                data.p=rmfield(data.p, 'K');
            end
            try
                data.p=rmfield(data.p, 'R');
            end
            try
                data.p=rmfield(data.p, 'node');
            end
        elseif strcmp(data.p.method, 'None')==1
            set([ edit_R,  edit_node, edit_K], 'visible', 'off');
        elseif strcmp(data.p.class, 'FM')==1 &&  strcmp(data.p.regime, 'Adiabatic')==1
            set([ edit_R,  edit_node, edit_K], 'visible', 'off');
            if strcmp(data.p.name, 'chirp')==1
                set(edit_node, 'visible', 'on')
            elseif strcmp(data.p.name, 'WURST')==1 || strcmp(data.p.name, 'HSn')==1
                if strcmp(data.p.method, 'Fixed_B1')
                    set(edit_K, 'visible', 'on')
                elseif strcmp(data.p.method, 'Fixed_opBW')==1||strcmp(data.p.method, 'Fixed_Tp')==1
                    set(edit_R, 'visible', 'on')
                end
            else
                disp('unsupported Adiabatic pulse name')
            end
        else
            set([ edit_R,  edit_node, edit_K], 'visible', 'off');
        end
        
        guidata(Source, data)
    end

%__________________________________________________________
    function hAnalyze_Callback(Source,eventdata)            
        data=guidata(Source);
        npts=round(2*data.sim.fmax/data.sim.df);
    
        Exp.npts=npts;
        Exp.MaxDetuning=data.sim.fmax*1e6;       
        Exp.T1=data.sim.T1;
        Exp.T2=data.sim.T2;
        Sys.M0=data.sim.M0;
        
        hcurr=findall(data.AllInputHandles, 'Visible', 'on');
        for jj=1:length(hcurr)
            hi=hcurr(jj);
            feval(get(hi, 'Callback'), hi);
        end
        
        pulse=data.p;
        pulse=Create_Optimized_Pulse(pulse);
        
        % check to see if there is a sufficient number of points
        dt=pulse.t(2)-pulse.t(1);
        df=1/dt;
        middle=round(length(pulse.esd)/2);
        cc=pulse.esd(middle:end);
        ii=find(abs(cc)<0.1);
        fmax=abs(pulse.esdf(middle+ii(1)));
        if df/fmax<10
        hwarn=warndlg('Insufficent number of timepoints.Errors may occur. Increase pulse.npts.');
        waitfor(hwarn);
        end
       
        
        %
        if strcmp(data.concatenation, 'BIR4')==1
            pulse.AFP=pulse;
        end
        
        pulse.concatenation=data.concatenation;
        pulse=data.composite(pulse);
        
        [ Sys, Exp, Mx, My, Mz ] = feval('EvolveM_Fast', pulse, Sys, Exp);
        [phi,theta,RR]=cart2sph(Mx, My, Mz);
        theta=-1*(theta-pi/2);
        
        
        %
        data.sim.Exp=Exp;
        data.sim.Sys=Sys;
        data.Simpulse=pulse;
        data.out.Mx=Mx;
        data.out.My=My;
        data.out.Mz=Mz;
        data.out.phi=phi/pi;
        data.out.theta=theta/pi;
        data.out.f=Exp.detuning/1e6;
        
        guidata(Source, data);
        
        plotdata([], []);
        try
            ShowSpec_Callback(Source, []);
        end

        op=theta(round(npts/2))/pi;
        [opbw, fpeak]=fwhm(Exp.detuning, theta);
        
        set(out_theta, 'String',['Rotation angle theta (pi)=', num2str(op, 2)]);
        set(out_opBW, 'String',['Operation BW (MHz)=',num2str(round(opbw/1e6))]);
        set(out_B1, 'String',['B1max (Gauss)=', num2str(pulse.B1max/1e-4, 2)]);
        set(out_Tp, 'String',['Tp pulse length (ns)=',num2str(round(pulse.t(end)/1e-9))]);
        
        if strcmp(pulse.class, 'AM')==1
            set(out_BWP, 'String',['Time-opBW product =',num2str(pulse.t(end)*opbw, 3)]);
            set(out_R, 'String',['Time-BW product R =', num2str(pulse.FWHM*pulse.t(end),3)]);
            set(out_Q, 'String','Adiabaticity Q = NA');
            set(out_fmax, 'String','Freq sweep limit (MHz)= 0');
        else
            set(out_BWP, 'String',['Time-opBW product =',num2str(pulse.t(end)*opbw, 3)]);
%             set(out_R, 'String',['Time-BW product R =',num2str(round(pulse.R))]);
            set(out_R, 'String',['Time-BW product R =',num2str(pulse.R, 3)]);
            set(out_Q, 'String',['Adiabaticity Q  =',num2str(pulse.Qcrit, 3)]);
            set(out_fmax, 'String',['Freq sweep limit (MHz)=',num2str(round(pulse.f_max/1e6))]);
       
        end        
        
            set(edit_theta, 'String',num2str(pulse.theta/pi, 3));
            try
            set(edit_R, 'String',num2str(pulse.R, 3));
            end
            try
            set(edit_node, 'String',num2str(pulse.node, 3));
            end
            try
            set(edit_K, 'String',num2str(pulse.K, 3));
            end
          
     end

%% ___________________composite _________________
    function Composite_selection(Source,eventdata)
        data=guidata(Source);
   
        choice=get(eventdata.NewValue,'tag');
        data.concatenation=choice;
        
        switch choice
            case 'None'
               data.composite=@(pulse) pulse;
               set([edit_BIRphase, txt_BIRphase],'visible', 'off')
                set(edit_theta,'Enable', 'on')
            
            case 'PhaseTime'
                 data.composite=@(pulse)Concatenate_PhaseTime(pulse);
%                  set(edit_BIRphase,'visible', 'off')
                 set([edit_BIRphase, txt_BIRphase],'visible', 'off')
                 set(edit_theta,'Enable', 'on')
            
            case 'BIR4'                           
                  set(edit_theta, 'String', '1')
                  set(edit_theta,'Enable', 'off') 
%                   set(edit_BIRphase,'visible', 'on')
                  set([edit_BIRphase, txt_BIRphase],'visible', 'on')
                  data.composite=@(pulse)Concatenate_BIR4(pulse); 
        end
       
       guidata(Source, data)              
    end

%% ______ pulse Input callbacks__________________________________________________________________--
% _______________________________________________________________
    function hRegime_Callback(hObject,eventdata)
        data=guidata(hObject);
        val = get(hObject,'Value');
        if val==1
             data.p.regime='Adiabatic';      
        else
            data.p.regime='RapidPassage';
           set([edit_R, edit_K, edit_node], 'visible', 'off')
        end
        guidata(hObject,data);
        MethodSelect_Callback(hMethodSelect, []);  
    end
% _______________________________________________________________
    function BIRphase_Callback(Source, eventdata)
        data=guidata(Source);
        
        BIRphase = str2num(get(Source,'String'));
        if isnan(BIRphase)
            errordlg('input must be numeric')
            return
        else
            data.p.BIRphase=BIRphase*pi;
            guidata(Source,data);
        end       
    end
%_________________________________________________________
    function fuw_Callback(Source, eventdata)
        data=guidata(Source);
        fuw= str2num(get(Source,'String'));
        if isnan(fuw)
            errordlg('input numeric number (GHz)')
            uicontrol(Source)
            return
        else
            data.Spec.fuw=fuw*1e9;
            guidata(Source,data);
            plotSpectrum(Source);
        end
    end
%_________________________________________________________
    function B0_Callback(Source, eventdata)
        data=guidata(Source);
        B0= str2num(get(Source,'String'));
        if isnan(B0)
            errordlg('input numeric number (Gauss)')
            uicontrol(Source)
            return
        else
            data.Spec.B0=B0*1e-4;
            guidata(Source,data);
            plotSpectrum(Source)
        end
    end
%___________________________________________________________---
    function opBW_Callback(Source, eventdata)
        data=guidata(Source);
        
        opBW = str2num(get(Source,'String'));
        if isnan(opBW)
            errordlg('input numeric number (MHz)')
            uicontrol(Source)
            return
        else
            data.p.opBW=opBW*1e6;
            guidata(Source,data);
        end
        
    end
%______________________________________________________________
    function npts_Callback(Source, eventdata)
        data=guidata(Source);
        
        npts = str2num(get(Source,'String'));
        if isnan(npts)
            errordlg('input numeric number')
            return
        elseif floor(npts)~=npts
            errordlg('npts must ne integer')
            npts=round(npts);
            set(Source, 'String', num2str(npts));
        else
            data.p.npts=npts;
            guidata(Source,data);
        end
    end
% _______________________________________________________________
    function theta_Callback(Source, eventdata)
        data=guidata(Source);        
        theta = str2num(get(Source,'String'));
        if isnan(theta)
            errordlg('input must be numeric')
            return
        else
            data.p.theta=theta*pi;
        end
        guidata(Source, data);
        
    end
% _______________________________________________________________
    function f0_Callback(Source, eventdata)
        data=guidata(Source);
        
        f0 = str2num(get(Source,'String'));
        if isnan(f0)
            errordlg('center frequency: input numeric number (MHz)')
            return
        else
            data.p.f_center=f0*1e6;
            guidata(Source,data);
        end
        
    end
% _______________________________________________________________
    function Tp_Callback(Source, eventdata)
        data=guidata(Source);
        
        Tp = str2num(get(Source,'String'));
        if isnan(Tp)
            errordlg('pulse length: input numeric number (ns)')
            return
        elseif Tp<=0
            errordlg('Tp must be >0')
            return
        else
            data.p.Tp=Tp*1e-9;
            guidata(Source,data);
        end
        
    end
% _______________________________________________________________
    function n_Callback(Source, eventdata)
        data=guidata(Source);
        
        n = str2num(get(Source,'String'));
        if isnan(n)
            errordlg('input numeric number for pulse.n')            
            return
        elseif floor(n)~=n
            errordlg('n must be an integer')
            n=round(n);
             set(edit_n, 'String',num2str(n));
        else
            data.p.n=n;
            
            try
            if strcmp(data.p.regime, 'Adiabatic')==1
                if strcmp(data.p.name, 'HSn')==1
                    if data.p.n>10
                        data.p.n=10;
                        set(edit_n, 'String',num2str(10))
                    end
                elseif strcmp(data.p.name, 'WURST')==1
                    if data.p.n==1
                        data.p.n=2;
                        set(edit_n, 'String',num2str(2))
                    end
                    if data.p.n>20
                        data.p.n=20;
                        set(edit_n, 'String',num2str(20))
                    end
                end
                
            end
            end 
            
            guidata(Source,data);
        end
        
    end
% _______________________________________________________________
    function B1_Callback(Source, eventdata)
        data=guidata(Source);
        B1max = str2num(get(Source,'String'));
        if isnan(B1max)
            errordlg('input numeric number for pulse.B1max')
            
            return
        else
            data.p.B1max=B1max*1e-4;
            guidata(Source,data);
        end
    end
% _______________________________________________________________
    function wmax_Callback(Source, eventdata)
        data=guidata(Source);
        fmax = str2num(get(Source,'String'));
        if isnan(fmax)
            errordlg('input numeric number for pulse.fmax')
            
            return
            
        else
            data.p.f_max=fmax*1e6;
            guidata(Source,data);
        end
    end
% _______________________________________________________________
    function sigma_Callback(Source, eventdata)
        data=guidata(Source);
        
        sigma = str2num(get(Source,'String'));
        if isnan(sigma)
            errordlg('input numeric number for pulse.sigma')           
            return
         elseif sigma<=0 || sigma>1
            errordlg('0<sigma<=1')           
            return
        else
            data.p.sigma=sigma;
            guidata(Source,data);
        end
        
    end
%_____________________________________________________________
    function K_Callback(Source, eventdata)
        data=guidata(Source);
        
        K = str2num(get(Source,'String'));
        if isnan(K)
            errordlg('input numeric number for K')
            
            return
        elseif K<=0
            errordlg('K cannot be <0')
            
            return
        else
            data.p.K=K;           
            guidata(Source,data);
        end
        
    end
%_____________________________________________________________
    function node_Callback(Source, eventdata)
        data=guidata(Source);
        
        R = str2num(get(Source,'String'));
        if isnan(R)
            errordlg('input numeric number for solution node integer')
            
            return
        elseif R<=0
            errordlg('node cannot be <0')
            
            return
        else
            data.p.node=R;           
            guidata(Source,data);
        end
        
    end
% _______________________________________________________________
    function R_Callback(Source, eventdata)
        data=guidata(Source);
        
        R = str2num(get(Source,'String'));
        if isnan(R)
            errordlg('input numeric number for time-bandiwdth product')
            
            return
        elseif R<=0
            errordlg('R cannot be <0')
            
            return
        else
            data.p.R=R;           
            guidata(Source,data);
        end
        
    end
% _______________________________________________________________
    function phi0_Callback(Source, eventdata)
        data=guidata(Source);
        
        phi0 = str2num(get(Source,'String'));
        if isnan(phi0)
            errordlg('input numeric number for phi0')
            
            return
        else
            data.p.phi0=phi0*pi;
            guidata(Source,data);
        end
        
    end
%__________________________________________________________________________
    function M0_Callback(Source, ~)
        data=guidata(Source);
        
        M0 = str2num(get(Source,'String'));
        if isnan(M0)
            errordlg('You must enter a vector like [0,0,1]')
            return
        elseif length(M0)~=3
            errordlg('You must enter a vector like [0,0,1]')
            return
        elseif norm(M0)>1
            errordlg('|M0|>1. This is unphysical')
            return 
        else 
            data.sim.M0=M0;
            guidata(Source,data);
        end
        
    end
%_____________________________________________________________
    function df_Callback(Source, ~)
        data=guidata(Source);
        df = str2num(get(Source,'String'));
        if isnan(df)
            errordlg('You must enter a numeric value')
            
            return
        elseif df<=0
            errordlg('df cannot be <0')
            
            return
        else
            data.sim.df=df;
            guidata(Source, data)
        end
    end
%__________________________________________________________--
    function fmax_Callback(Source,~)
        data=guidata(Source);
        fmax = str2num(get(Source,'String'));
        if isnan(fmax)
            errordlg('You must enter a numeric value')
            
            return
        elseif fmax<=0
            errordlg('fmax cannot be <0')
            
            return
        else
            data.sim.fmax=fmax;
            guidata(Source, data)
        end
    end
%__________________________________________________________--
    function T1_Callback(Source,~)
        data=guidata(Source);
        T = str2num(get(Source,'String'));
        if isnan(T)
            errordlg('You must enter a numeric value')
            
            return
        elseif T<=0
            errordlg('T cannot be <0')
            
            return
        else
            data.sim.T1=T;
            guidata(Source, data)
        end
    end
%__________________________________________________________--
    function T2_Callback(Source,~)
        data=guidata(Source);
        T = str2num(get(Source,'String'));
        if isnan(T)
            errordlg('You must enter a numeric value')
            
            return
        elseif T<=0
            errordlg('T cannot be <0')
            
            return
        else
            data.sim.T2=T;
            guidata(Source, data)
        end
    end

%% other callbacks
 function hInfo_Callback(Source, eventdata)
  figI=figure('color', 'white', 'units', 'norm', 'pos', [.3, .4, .2, .3], 'Name', 'Software Information');
  
  string={'Software Version: ShapedPulseSimulator 1.1', ...
       '',...
      'Release Date: 4/12/2017', ...
      '',...
      'Author: Joanna A. Guse',...
       '',...
      'Licence:If you use results obtained with the help of ShapedPulseSimulator in any scientific publication, cite the appropriate publication. Licence and citation information can be found at www.shapedpulsesimulator.org',...
       '',...
      'Documentation: Documentation and examples can be found at www.shapedpulsesimulator.org'};
  

  normpos=[.0, .0, .99, .99];
  infotxt  = uicontrol('Parent', figI,'tag', 'info', 'Style','text','String',string,...
            'Visible','on','BackgroundColor','white','units', 'normalized','Fontsize', 12, ...
            'position', normpos, 'HorizontalAlignment', 'Left');

     
 end

%_______________________________________________________________
    function PlotEnv_Callback(Source, eventdata)
        
        data=guidata(Source);
        hcurr=findall(data.AllInputHandles, 'Visible', 'on');
        for jj=1:length(hcurr)
            hi=hcurr(jj);
            feval(get(hi, 'Callback'), hi);
        end
        
        pulse=data.p;
        if isfield(pulse, 't')==0
            hAnalyze_Callback(hAnalyze,[])
            pulse=data.Simpulse;
        end
        
        
        fig2=figure('color',[1,1,1], 'units', 'normalized', 'position', [.1, .1, .7, .4]);
        subplot(2,2,1)
        plot(pulse.t/1e-9, pulse.env/1e-4, 'linewidth',2)
        xlabel('Time (ns)', 'Fontsize', 14)
        ylabel('AM (Gauss)', 'Fontsize', 14)
        set(gca, 'Fontsize', 14)
        subplot(2,2,2)
        plot(pulse.t/1e-9, pulse.f_mod/1e6, 'linewidth',2)
        xlabel('Time (ns)', 'Fontsize', 14)
        ylabel('FM (MHz)', 'Fontsize', 14)
        set(gca, 'Fontsize', 14)
        title([pulse.name ' modulation functions'],  'Fontsize', 16)
        subplot(2,2,3)
        plot(pulse.t/1e-9, pulse.phi/pi, 'linewidth',2)
        xlabel('Time (ns)', 'Fontsize', 14)
        ylabel('\phi (\pi)', 'Fontsize', 14)
        set(gca, 'Fontsize', 14)
        subplot(2,2,4)
        plot(pulse.esdf/1e6, pulse.esd, 'linewidth',2)
        xlabel('Freq(MHz)', 'Fontsize', 14)
        ylabel('|fft|^2', 'Fontsize', 14)
        set(gca, 'Fontsize', 14)
        try
            xlim(round([-5*pulse.FWHM, 5*pulse.FWHM]/1e6))
        end
    end
%_______________________________________________________________-
    function Coordinate_selection(Source,eventdata)
        data=guidata(Source);
        switch get(eventdata.NewValue,'tag')
            case '1'
                data.coord='cartesian';
             
            case '2'
                data.coord='polar';
        end
        
guidata(Source, data)  
plotdata
ShowSpec_Callback(hShowSpec, []);              
    end
%______________________________________________________________
    function hExportWav_Callback(Source, eventdata)
        data=guidata(Source);
        pulse=data.Simpulse;
        Export_wavfileGUI3(pulse)
    end
%_____________________________________________
    function hExport_Callback(Source, ~)
        data=guidata(Source);
        pulse=data.Simpulse;
        output=data.out;
        inputs=data.sim;
        inputs.p=data.p;
        [file,path] = uiputfile('*.mat','Save Simulation As');
        try
            save([path, file],'output', 'inputs', 'pulse');
        end
    end
% _______________________________________________________________
    function hNewPlot_Callback(hObject,~)
        GUI_fig_children=get(gcf,'children');
        Fig_Axes=findobj(GUI_fig_children,'type','Axes');
        fig=figure('color', 'white', 'units', 'normalized','position', [.1,.1,.5,.7]);
        new_handles=copyobj(Fig_Axes,fig);
        if length(new_handles)==2
            set(new_handles(2), 'Position',[.1, .1, .85, .4])
            set(new_handles(1), 'Position',[.1, .55, .85, .4])
        elseif length(new_handles)==3
            set(new_handles(1), 'Position',[.1, .1, .2, .1])
            set(new_handles(2), 'Position',[.1, .1, .85, .4])
            set(new_handles(3), 'Position',[.1, .55, .85, .4])
        end
        set(new_handles, 'Fontsize', 14)
    end


%% my GUI functions_______________________________________________________________

    function [FWHM,t0] = fwhm(t,y,range)
        % This function computes the intensity full-width at half maximum
        % (FWHM) of the periodic signal y.
        
        nt = length(y);
        dt = diff(t(1:2));
        k0 = (1:nt)';
        kl = [nt,1:nt-1]';
        kr = [2:nt,1]';
        
        if nargin < 3
            [~,ipeak] = max(y);
        else
            if length(range) == 2
                kv = find(range(1) < t & t < range(2));
                [~,ipeak] = max(y(kv));
                ipeak = k0(kv(ipeak));
            elseif length(range) == 1
                itr = interp1(t,k0,range);
                localmax = find(y(k0) >= y(kl) & y(k0) > y(kr));
                idiff = mod(localmax - itr + nt/2 - 1, nt) - nt/2 + 1;
                [~,i] = min(abs(idiff));
                ipeak = localmax(i);
            end
        end
        
        % 3 point quadratic fit to extract peak value
        
        pv = [kl(ipeak),k0(ipeak),kr(ipeak)];
        coefs = polyfit((-1:1),[y(pv(1)),y(pv(2)),y(pv(3))],2);
        coefs1 = polyder(coefs);
        idt = roots(coefs1);
        t0 = t(ipeak) + idt*dt;
        umax = polyval(coefs,idt);
        
        % find first rising half-max point to left of ipeak
        
        irise = find(y(k0) <= umax/2 & y(kr) > umax/2);
        irise = irise + (umax/2 - y(k0(irise)))./(y(kr(irise)) - y(k0(irise)));
        irise = min(mod(ipeak-irise,nt));
        
        % find first falling half-max point to right of ipeak
        
        ifall = find(y(k0) >= umax/2 & y(kr) < umax/2);
        ifall = ifall + (y(k0(ifall)) - umax/2)./(y(k0(ifall)) - y(kr(ifall)));
        ifall = min(mod(ifall-ipeak,nt));
        
        FWHM = dt*(irise+ifall);
    end
%__________________________________________________________
    function plotdata(varargin)
        % ____________delete all traces
        data=guidata(hfig);
        h1=get(data.haxes(1), 'Children');
        for jj=1:length(h1)
            delete(h1(jj));
        end
        h2=get(data.haxes(2), 'Children');
        for jj=1:length(h2)
            delete(h2(jj));
        end
        
        hLeg=findobj(gcf,'tag','legend');
        delete(hLeg);
        
        
fa=12;
fl=14;

        %______ plot new traces
        if strcmp(data.coord, 'cartesian')==1
            axes(data.haxes(1));
            plot( data.out.f, data.out.Mz,'b', 'Linewidth', 2);

            axis tight
            ylabel('M_z', 'fontsize', fl)
            ylim([-1, 1])
            set(gca, 'fontsize', fa)
            
            
            axes(data.haxes(2));
            plot(data.out.f, data.out.Mx, 'b', 'Linewidth', 2);
            hold on
            plot(data.out.f, data.out.My, 'r', 'Linewidth', 2);

            axis tight
            ylabel('M_x and M_y', 'fontsize', fl)
            legend('Mx', 'My','location', 'best');
            xlabel('Freq (MHz)', 'fontsize', fl)
            ylim([-1, 1]);
            yticklabel={'-1', '-0.5', '0','0.5','1'};
            yticks=[-1, -0.5, 0, 0.5, 1];
            set(gca,'YTick',yticks)
            set(gca,'YTickLabel',yticklabel)
            set(gca, 'fontsize', fa)
            
        else
            axes(data.haxes(1));
            plot(data.out.f, data.out.theta,'b', 'Linewidth', 2);
            axis tight
            ylabel('\theta (\pi)', 'fontsize', fl)
            ylim([0, 1]);
            set(gca, 'fontsize', fa)
            
            axes(data.haxes(2));
            plot(data.out.f, data.out.phi, 'b','Linewidth', 2);

            axis tight
            ylim([-1, 1])
            axis([data.out.f(1), data.out.f(end), -1.1, 1.1])
            yticks=[-1, -0.5, 0, 0.5, 1];
            yticklabel={'x', '-y', '-x','y'};
            set(gca,'YTick',yticks)
            set(gca,'YTickLabel',yticklabel);
            set(gca, 'fontsize', fa)
            
            xlabel('Freq (MHz)', 'fontsize', fl)
            ylabel('\phi (\pi)', 'fontsize', fl)
        end
        
        
%     
% inds=find(abs(data.out.f)<0.5*data.sim.fmax);
% figure
% A=phi(inds);
% B=RR(inds);
% cc = redblue(length(A));
% for ii = 1:length(A)
%     L(ii) =polar(A(ii),B(ii), 'o');
%     hold on
%     set(L(ii),'markerfacecolor',cc(ii,:),...
%         'markeredgecolor', cc(ii,:))
%     legendinfo{ii}=num2str(round(Exp.detuning(inds(ii))/1e6));
% end
% legend(L(1:round(length(L)/6):length(L)),legendinfo{1:round(length(L)/6):length(L)})
% 
% 
% hHiddenText = findall(gca,'type','text');
% set(hHiddenText, 'visible', 'off')
% fl=14;
% text(1.1*max(B),0,'+y', 'Fontsize', fl)
% text(0,1.1*max(B),'+x',  'Fontsize', fl)
% text(-1.2*max(B),0,'-y',  'Fontsize',fl)
% text(0,-1.1*max(B),'-x',  'Fontsize', fl)
%         
    guidata(hfig, data);
    end
%________________________________________________-
    function [ htxt, hedit] = Make_TxtEditpair( parent, name,string, normpos, color, defaultVal,boxscale, callbackfunc, vis)
        
        htxt  = uicontrol('Parent', parent,'tag', name, 'Style','text','String',string,...
            'Visible', vis,'BackgroundColor',color,'units', 'normalized', 'position', normpos);
        pos=get(htxt, 'position');
        pos2=[pos(1)+pos(3), pos(2),0.4, pos(4)];
        hedit = uicontrol('Parent', parent,'tag', name,'Style','edit','units', 'normalized',...
            'Position',[pos(1)+pos(3), pos(2),0.2, boxscale*pos(4)],'String', defaultVal, 'BackgroundColor','white',...
            'Visible', vis,'Callback', callbackfunc);
        clear pos
        align([htxt, hedit],'VerticalAlignment','Top');
        set([htxt, hedit],'HorizontalAlignment','Left');
        
    end
%________________________________________________-
    function [ htxt, hedit] = Make_TxtEditpair2( parent, name,string, normpos,editw, color, defaultVal,boxscale, callbackfunc, vis)
        
        htxt  = uicontrol('Parent', parent,'tag', name, 'Style','text','String',string,...
            'Visible', vis,'BackgroundColor',color,'units', 'normalized', 'position', normpos);
        pos=get(htxt, 'position');
        hedit = uicontrol('Parent', parent,'tag', name,'Style','edit','units', 'normalized',...
            'Position',[pos(1)+pos(3), pos(2),editw, boxscale*pos(4)],'String', defaultVal, 'BackgroundColor','white',...
            'Visible', vis,'Callback', callbackfunc);
        clear pos
        align([htxt, hedit],'VerticalAlignment','Top');

    end

%________________________________________________
function [ Phi ] =UnwrapPhase(x, phi)
%Unwraps Phase a from cart2sph(Mx, My)

phi=unwrap(phi);
phi=phi-min(phi);
diff(1)=0;
ind(1)=1;
 for jj=2:length(phi)
     diff(jj)=phi(jj)-phi(jj-1);
        if abs(diff(jj))>0.9
            ind(jj)=ind(jj-1)*-1;
        else
            ind(jj)=ind(jj-1);
        end
 end
ind=1-(ind+1)/2;
Phi=phi+pi*ind';
Phi=unwrap(Phi);

[~, ii]=min(abs(x));
Phi=Phi-Phi(ii);

end


%________________________________________________________________
    function [ htxt] = Make_Txt( parent, tagname,string, normpos, color, vis, justification, fontsize)
        htxt  = uicontrol('Parent', parent,'tag', tagname, 'Style','text','String',string,...
            'Visible', vis,'BackgroundColor',color,'units', 'normalized','Fontsize', fontsize, ...
            'position', normpos, 'HorizontalAlignment', justification);
    end
%________________________________________________________________
    function [ htxt] = Make_Txt2(tagname,string, normpos, color, vis, justification, fontsize)
        htxt  = text('position', [normpos(1), normpos(2)],'String',string);
        set(htxt,'HorizontalAlignment', justification, 'tag', tagname,...
            'Visible', vis,'BackgroundColor',color,'units', 'normalized','Fontsize', fontsize);
    end
%________________________________________
    function [  ] = showhide( show,hide)
        for jj = 1:length(hide)
            temp= findobj('Tag', hide{jj});
            set(temp, 'Visible', 'off');
        end
        
        for ii =1:length(show)
            temp= findobj('Tag', show{ii});
            set(temp, 'Visible', 'on');
        end
        
    end

%% spectrum callbacks
%___________________________________________________
    function hImportSpec_Callback(Source, ~)
        
        data=guidata(Source);
        hh=ImportSpecGUI5;
        uiwait(hh)
        a=load('tempSpec.mat');
        data.Spec=a.specData;
        
        guidata(Source, data);
        
        delete('tempSpec.mat')
        set(txt_filename, 'String', data.Spec.fname)
        
        set(edit_fuw, 'String',data.Spec.fuw/1e9)
        set(edit_B0, 'String',data.Spec.B0/1e-4)
           
        set(hShowSpec, 'visible', 'on', 'Value', 1)
        ShowSpec_Callback(hShowSpec, []);
        
        guidata(Source, data);
    end
%_______________________________________________
    function ShowSpec_Callback(Source, eventdata)
        data=guidata(Source);
        val=get(Source, 'Value');
        if val==0
            try
            delete( data.specax(1));
            delete( data.specax(2));
            end
        else    
         plotSpectrum(Source); 
        end
guidata(Source, data);
    end
%______________________________________________________________
    function plotSpectrum(Source)
        data=guidata(Source);            
             
        try
              %  delete old traces
              delete( data.specax(1));
              delete( data.specax(2));
        end   
              % rescale to correct coordinate
              spec=data.Spec.y-min(data.Spec.y);
              spec=spec./max(spec);
              
          if strcmp(data.coord, 'cartesian')==1
              spec1=2*spec-1;
              spec2=2*spec-1;
          else
               spec1=spec;
               spec2=2*spec-1;
          end
                    
          % rescale X axis
        mub=9.274009994*1e-24;
        h=6.626e-34;       
        gii=h*data.Spec.fuw./(mub*data.Spec.x);
        F=gii*mub*data.Spec.B0/h;
        F=(F-data.Spec.fuw)/1e6;      
        data.Spec.freq=F;        
        
        % plot         
                axes(data.haxes(1));
               
                hold on
                data.specax(1)=plot(data.Spec.freq,spec1, 'k','Linewidth', 2);   
                axis tight
                xlim([-data.sim.fmax, data.sim.fmax]);
                axes(data.haxes(2));
                hold on
                data.specax(2)=plot(data.Spec.freq,spec2, 'k','Linewidth', 2);  
                axis tight 
                xlim([-data.sim.fmax, data.sim.fmax]);
        guidata(Source, data)
    end

%% visible on
set(hfig, 'Visible', 'on')
end

%%______________________________________________________________
%%________________________________________________________________
function [ ] = Export_wavfileGUI3( varargin )
%this GUI exports pulse data to .wav files according to specified formatting inputs.
%Export_wavfileGUI() will result in a loading data window.
%Export_wavfileGUI(pulse)
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
        if exists==0
            disp('no structure called "pulse" in input file')
        end
        data.pulse=pulse;
        
        
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


%%
%%
%%
function [hfig1]=ImportSpecGUI5(varargin)
hfig1 = figure('Visible','off','color','w','Position',[0,0, 400, 200]);
set(hfig1, 'Name', 'Import Spectrum')
movegui(hfig1,'center')

author = uicontrol('Parent', hfig1, 'Style','text','String','Joanna A. Guse',...
    'Visible', 'on','BackgroundColor','white','units', 'normalized', 'position', [.74, 0., .3 .09]);
set(hfig1, 'Visible', 'on')

specData=guidata(hfig1);
specData.B0=3440e-4;
specData.fuw=9.7e9;
guidata(hfig1, specData);
xx=[];
yy=[];

%% Select format
hBrowse= uicontrol('parent', hfig1,'style', 'pushbutton',...
    'units','normalized', 'position',[.1, .8, .25, .15],...
    'string','Browse', 'callback', @BrowseInput);
txt_file = uicontrol('Parent', hfig1, 'Style','text','String','filename',...
    'Visible', 'on','BackgroundColor',[.8, .8, .8],'units', 'normalized',...
'position', [.35, .8, .6 .15]);

%% parameter panel
hImportParams = uipanel('parent', hfig1,'Title','Import parameters',...
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
guidata(hfig1, specData);


hImportButton= uicontrol('parent',hfig1,'style', 'pushbutton', ...
    'units','normalized', 'position',[.1, .1, .3, .15],...
    'string','Import Spectrum', 'callback', @ImportSpec_callback);
Note = uicontrol('Parent', hfig1, 'Style','text','String','Note: figure will close automatically when import is done. Do not close figure',...
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
        save('tempSpec.mat', 'specData')
        guidata(Source, specData);
         pause(1)
         close(hfig1)
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
        specData=guidata(hfig1);
        specData.fuw=double(get_vb(VB,'MWFQ'));
        specData.B0=double(get_vb(VB,'A1CT'));
        guidata(hfig1, specData);
        
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
