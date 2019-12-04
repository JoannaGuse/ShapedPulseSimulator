function pulse = Def_WURST(pulse)
% WURST definition function (FM).
% pulse shaping parameters: n Valid for n>2
% Methods='Fixed_opBW', 'Fixed_B1', 'Fixed_Tp','None'
% Regimes: 'Adiabatic', 'RapidPassage'
% Adiabatic design parameters: pulse.R ( for 'Fixed_opBW' and 'Fixed_Tp')
%                              pulse.K (for 'Fixed_B1')
% theta: variable for 'RapidPassage', and [pi/2 or pi] for 'Adiabatic'
% regime
pulse.class='FM';
pulse.descr='Wideband,uniform-rotaton,smooth-truncation adiabatic pulse. AM and FM. Best suited for Inversion';
pulse.DynamicVars={'n', 'theta'};
pulse.SupportedMethods={'Fixed_opBW', 'Fixed_B1', 'Fixed_Tp','None'};

%set up defaults
if isfield(pulse, 'theta')==0
    pulse.theta=pi;
end
if isfield(pulse, 'n')==0
    pulse.n=10;
end
    if pulse.n==1
        pulse.n=2;
    end

if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

% add F1 functions
pulse.F1=F1_def(pulse);
pulse.F2=F2_def(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);

%_________ Add fucntion handles to do operation
pulse.SimMethod=@(pulse)SimMethod(pulse);
end
%________________________________________
function F1=F1_def(pulse)
tau=linspace(-1, 1, pulse.npts);
F1=1-abs((sin(pi/2*tau)).^pulse.n);
end
%________________________________________
function F2=F2_def(pulse)
tau=linspace(-1, 1, pulse.npts);
F2=tau;
end

%________________________________________________________
function pulse=SimMethod(pulse)
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
            if  isfield (pulse, 'R')==0
                if pulse.n==1
                    pulse.R=50;
                else
                    pulse.R=20;
                end
            end
            
            if pulse.R<15 
                pulse.R=15;
            end
 % semi-analytical approximations. Can be modified by user to increase
 % fidelity
            if pulse.theta==pi
                pulse.eta=1.081*(1+pulse.n^-0.6748);
            elseif pulse.theta==pi/2
                pulse.eta=1.14+1.561*pulse.n^-0.8658;
            else
                pulse.eta=1.081*(1+pulse.n^-0.6748);
            end
            
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
            
            if pulse.theta==pi
                if isfield(pulse, 'K')==0
                    pulse.K=0.5;
                elseif pulse.K<0.5
                    pulse.K=0.5;
                end

 % semi-analytical approximations for xi. Can be modified by user to increase
 % fidelity
                pulse.xi=1;
            elseif pulse.theta==pi/2
                if isfield(pulse, 'K')==0
                    pulse.K=1.5;
                elseif pulse.K<1
                    pulse.K=1;
                end
                pulse.xi=1+0.4*pulse.n^(-8.6);
            end
            
            pulse.f_max=(gamma*pulse.B1max)*pulse.K;
            Q=-2/pi*log(1-0.999*pulse.theta/pi);
            pulse.Tp=pulse.xi*2*pi*2*pulse.f_max*Q/(gamma*pulse.B1max).^2;
            pulse.t=linspace(0, pulse.Tp ,pulse.npts);
            pulse.R=pulse.Tp*pulse.f_max*2;
            
        case 'Fixed_Tp'
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            if  isfield (pulse, 'Tp')==0
                pulse.Tp=200e-9;
            end
            if  isfield (pulse, 'R')==0
                if pulse.n==1
                    pulse.R=50;
                else
                    pulse.R=20;
                end
            end

            if pulse.R<15
                pulse.R=15;
                disp('WARNING: R must be >15 to achieve acceptable rotation angle fidelity. Defualt R=20;')
            end
            
            if pulse.n>20
                pulse.n=10;
            end
          
 % semi-analytical approximations. Can be modified by user to increase
 % fidelity
           if pulse.theta==pi
             pulse.kappa=1.08*(1+pulse.n^-0.68);
            elseif pulse.theta==pi/2                
                pulse.kappa=1;
            else
                pulse.kappa=1;
           end

            pulse.f_max=pulse.R/(2*pulse.Tp);
            Q=-2/pi*log(1-0.999*pulse.theta/pi);
            pulse.B1max=pulse.kappa*1/gamma*sqrt((4*pi*pulse.f_max*Q)/(pulse.Tp));
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
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            if  isfield (pulse, 'opBW')==0
                pulse.opBW=40e6;
            end
            
            pulse.R=0.3;
            eta0=(0.351*pulse.n^2+3.481*pulse.n+2.796)/(pulse.n^2+13.08*pulse.n+19.53);
            if pulse.theta==pi
                pulse.eta=eta0;
            elseif pulse.theta==pi/2
                pulse.eta=0.7166*eta0;
            else
               pulse.eta=eta0;
            end
            pulse.f_max=pulse.eta*pulse.opBW/2;
            pulse.Tp=pulse.R/(2*pulse.f_max);
            pulse.t=linspace(0, pulse.Tp, pulse.npts);
            pulse.B1max=pulse.theta/(gamma*pulse.Tp*pulse.RelArea);
            
            
            %___________________________________________________
        case 'Fixed_B1'
            if isfield(pulse, 'B1max')==0
                pulse.B1max=2e-4;
            end
            if isfield(pulse, 'Tp')==0
                pulse.Tp=40e-9;
            end
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            
            R=0.5;
            pulse.Tp=pulse.theta/(gamma*pulse.B1max*pulse.RelArea);
            pulse.f_max=R/(pulse.Tp*2);
            pulse.t=linspace(0, pulse.Tp ,pulse.npts);
            
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

