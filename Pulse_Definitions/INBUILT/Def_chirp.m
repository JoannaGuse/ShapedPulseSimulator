function pulse = Def_chirp(pulse)
% chirp definition function (FM). 
% pulse shaping parameters: -
% Methods='Fixed_opBW', 'Fixed_B1', 'Fixed_Tp','None'
% Regimes: 'Adiabatic', 'RapidPassage'
% Adiabatic design parameters : pulse.node valid for node<8 or <6 depending on
% method.
% theta: variable for 'RapidPassage', and [pi/2 or pi] for 'Adiabatic'
% regime
pulse.class='FM';
pulse.descr='Linear freqency adiabatic pulse. Best suited for Inversion';
pulse.DynamicVars={'theta'};
pulse.SupportedMethods={'Fixed_opBW', 'Fixed_B1', 'Fixed_Tp','None'};

%set up defaults
if isfield(pulse, 'theta')==0
    pulse.theta=pi;
end
if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

pulse.F1=F1_def(pulse);
pulse.F2=F2_def(pulse);

pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);

%_________ Add fucntion handles to do operation
pulse.SimMethod=@(pulse)SimMethod(pulse);

end

%________________________________________
function F1=F1_def(pulse)
tau=linspace(-1, 1, pulse.npts);
F1=ones(size(tau));
end
%________________________________________
function F2=F2_def(pulse)
tau=linspace(-1, 1, pulse.npts);
F2=tau;
end
%_____________________________________________________________________
function pulse= SimMethod( pulse)
gamma=1.760859e11;

if strcmp(pulse.regime, 'Adiabatic')==1
        Th=[pi, pi/2];
    if ismember(pulse.theta, Th)==0
        pulse.theta=pi;
    end
    
    switch pulse.method
       %_____________________________________________________________
        case 'Fixed_opBW'
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            if  isfield (pulse, 'opBW')==0
                pulse.opBW=40e6;
            end
            if  isfield (pulse, 'node')==0      
                      pulse.node=5;
            end            
            if pulse.node>8
                pulse.node=8;
            elseif pulse.node<1
                pulse.node=1;
            end
            
            % the following R values correspond to high fidelity regions.
            % May be adjusted to increase fidelity.
            Rs=[8.8, 12.4, 16, 20, 23.4, 27.4, 31, 35.4];
            %eta(R) is the banwidth correction factor. 
            if pulse.theta==pi/2
               eta=[1.19, 1.11, 1.075, 1.059, 1.054, 1.045, 1.06, 1.05];
            else
               eta=[1.14, 1.06, 1.03, 1.025, 1.005, 1.002, 0.996, 0.997];
            end
               
            nodes=1:length(eta);
            [~, ii]=min(abs(nodes-pulse.node));
            pulse.node=ii;
            pulse.eta=eta(ii);
            pulse.R=Rs(ii);

            pulse.f_max=pulse.eta*pulse.opBW/2;
            pulse.Tp=pulse.R/(2*pulse.f_max);
            Q=-2/pi*log(1-0.999*pulse.theta/pi);
            pulse.B1max=1/gamma*sqrt((4*pi*pulse.f_max*Q)/(pulse.Tp));
            pulse.t=linspace(0, pulse.Tp ,pulse.npts);
             
%_________________________________________________
        case 'Fixed_B1'
            if isfield(pulse, 'B1max')==0
                pulse.B1max=5e-4;
            end
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            
            if  isfield (pulse, 'node')==0   
                pulse.node=3;
            end
            if pulse.node>6
                pulse.node=6;
            elseif pulse.node<1
                pulse.node=1;
            end
            if pulse.theta==pi/2
                Ks=[0.87,1.2, 1.5, 2.1, 2.5, 2.9]; 
            else
                Ks=[0.5, 0.7, 1.5, 2.1, 2.5,2.9];
            end
            nodes=1:length(Ks);
            [~, ii]=min(abs(nodes-pulse.node));
            pulse.node=ii;
            pulse.K=Ks(ii);
            
            pulse.f_max=(gamma*pulse.B1max)*pulse.K;
            Q=-2/pi*log(1-0.999*pulse.theta/pi);
            pulse.Tp=2*pi*2*pulse.f_max*Q/(gamma*pulse.B1max).^2;
            pulse.t=linspace(0, pulse.Tp ,pulse.npts);
            pulse.R=pulse.Tp*pulse.f_max*2;

%____________________________________________-
               case 'Fixed_Tp'
            if isfield(pulse, 'Tp')==0
                pulse.Tp=200e-9;
            end
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            
            if  isfield (pulse, 'node')==0   
                pulse.node=10;
            end
            if pulse.node>8
                pulse.node=8;
            elseif pulse.node<1
                pulse.node=1;
            end
            Rs=[8.8, 12.4, 16, 20, 23.4, 27.4, 31, 35.4]; 
            nodes=1:length(Rs);
            [~, ii]=min(abs(nodes-pulse.node));
            pulse.node=ii;
            pulse.R=Rs(ii);
            
            pulse.f_max=pulse.R/(2*pulse.Tp);
            Q=-2/pi*log(1-0.999*pulse.theta/pi);
            pulse.B1max=1/gamma*sqrt((4*pi*pulse.f_max*Q)/(pulse.Tp));
            pulse.t=linspace(0, pulse.Tp ,pulse.npts);
                                     
        %___________________________________________________    
        case 'None'
            if  isfield (pulse, 'Tp')==0
                pulse.Tp=200e-9;
            end
            if isfield(pulse, 'B1max')==0
                pulse.B1max=5.6e-4;
            end
            if isfield(pulse, 'f_max')==0
                pulse.f_max=50e6;
            end
            
            pulse.t=linspace(0, pulse.Tp, pulse.npts);
            pulse.R=pulse.Tp*pulse.f_max*2;
        otherwise
            disp('unsupported method')
    end
    
    %_______________________________________________________
elseif strcmp(pulse.regime, 'RapidPassage')==1
    
    switch pulse.method
        %___________________________________________________
        case 'Fixed_opBW'
            %check defaults
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            if  isfield (pulse, 'opBW')==0
                pulse.opBW=40e6;
            end
             
            pulse.R=0.3;       
            if pulse.theta==pi
                 pulse.eta=0.374;  
            elseif pulse.theta==pi/2
                pulse.eta=0.268;
            else
                 pulse.eta=0.374;
            end 
            
            pulse.f_max=pulse.eta*pulse.opBW/2;
            pulse.Tp=pulse.R/(2*pulse.f_max);
            pulse.t=linspace(0, pulse.Tp, pulse.npts);
            pulse.B1max=pulse.theta/(gamma*pulse.Tp*pulse.RelArea);
            
            %___________________________________________________
        case 'Fixed_Tp'
            if isfield(pulse, 'Tp')==0
                pulse.Tp=40e-9;
            end
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end           
            R=0.5;
            pulse.B1max=pulse.theta/(gamma*pulse.Tp*pulse.RelArea);
            pulse.f_max=R/(pulse.Tp*2);
            pulse.t=linspace(0, pulse.Tp ,pulse.npts);
            %___________________________________________________
        case 'Fixed_B1'
            if isfield(pulse, 'B1max')==0
                pulse.B1max=2e-4;
            end
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            
            R=0.5;
            pulse.Tp=pulse.theta/(gamma*pulse.B1max*pulse.RelArea);
            pulse.f_max=R/(pulse.Tp*2);
            pulse.t=linspace(0, pulse.Tp ,pulse.npts);
       %___________________________________________________     
        case 'None'
            if  isfield (pulse, 'Tp')==0
                pulse.Tp=200e-9;
            end
            if isfield(pulse, 'B1max')==0
                pulse.B1max=5.6e-4;
            end
            if isfield(pulse, 'f_max')==0
                pulse.f_max=50e6;
            end
            
            pulse.t=linspace(0, pulse.Tp, pulse.npts);
            pulse.R=pulse.Tp*pulse.f_max*2;
        otherwise
            disp('unsupported method')
    end
   pulse.R=pulse.Tp*pulse.f_max*2;
else
    disp('unsupported regime')
end
end


