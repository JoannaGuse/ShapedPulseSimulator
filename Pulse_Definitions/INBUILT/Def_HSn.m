function pulse = Def_HSn(pulse)
% HSn definition function (FM). 
% pulse shaping parameters: n Valid for 1<=n<=10
% Methods='Fixed_opBW', 'Fixed_B1', 'Fixed_Tp','None'
% Regimes: 'Adiabatic', 'RapidPassage'
% Adiabatic design parameters : pulse.R ( for 'Fixed_opBW' and 'Fixed_Tp')
%                               pulse.K (for 'Fixed_B1')
% theta: variable for 'RapidPassage', and [pi/2 or pi] for 'Adiabatic'
% regime
pulse.class='FM';
pulse.descr='Hyperbolic secant adiabatic pulse. AM and FM. Best suited for Inversion';
pulse.DynamicVars={'n', 'theta'};
pulse.SupportedMethods={'Fixed_opBW', 'Fixed_B1','Fixed_Tp', 'None'};

%set up defaults
if isfield(pulse, 'theta')==0
    pulse.theta=pi;
end

    if isfield(pulse, 'n')==0
    pulse.n=1;
    end
    
    if pulse.n>10
    pulse.n=10;
    elseif pulse.n<1
    pulse.n=1;
    end
     
if isfield(pulse, 'npts')==0
    pulse.npts=1024;
end

pulse.beta=asech(0.01);

pulse.F1=F1_def(pulse);
pulse.F2=F2_def(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);

%_________ Add fucntion handles to do operation
pulse.SimMethod=@(pulse)SimMethod(pulse);

end
%________________________________________
function F1=F1_def(pulse)
tau=linspace(-1, 1, pulse.npts);
%________define F1 (AM function)____________
tau=linspace(-1, 1, pulse.npts);
F1=sech(pulse.beta*tau.^pulse.n);

end
%________________________________________
function F2=F2_def(pulse)
tau=linspace(-1, 1, pulse.npts);
%________define F2 (FM function)____________
F2=cumtrapz(tau, (sech(pulse.beta*tau.^pulse.n)).^2);
F2=F2/max(F2);
F2=1-(2*F2); %ensue F_max is in [-fmax, fmax].
F2=-F2;
end

%________________________________________________________--
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
                
            if pulse.theta==pi
                 Rmin=0.36*pulse.n^1.67+6.4;
                pulse.eta=0.975;
            elseif pulse.theta==pi/2
                 Rmin=20.29*sin(pulse.n*0.1562+0.013);
                pulse.eta=-0.25*pulse.n^-1.4+1.01;
            else
                 Rmin=0.36*pulse.n^1.67+6.4;
                pulse.eta=0.975;
            end
            if  isfield (pulse, 'R')==0
                pulse.R=Rmin;
            end
            if pulse.R<Rmin
                pulse.R=Rmin;
            end
            
           %Kappa cab ne adjusted to change fidelity. 
            kappa=(2.1*pulse.theta/pi);
            
            pulse.f_max=pulse.eta*pulse.opBW/2;
            pulse.Tp=pulse.R/(2*pulse.f_max);
            pulse.B1max=kappa*1/gamma*pulse.theta*pulse.beta^(1/(2*pulse.n))*pulse.f_max*2/sqrt(pulse.R);
            pulse.t=linspace(0, pulse.Tp ,pulse.npts);
            
            %_______________________________________________
        case 'Fixed_Tp'
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            if  isfield (pulse, 'Tp')==0
                pulse.Tp=200e-9;
            end
            
            if pulse.theta==pi
                Rmin=0.36*pulse.n^1.67+6.4;
            elseif pulse.theta==pi/2
                Rmin=20.29*sin(pulse.n*0.1562+0.013);
            else
                Rmin=0.36*pulse.n^1.67+6.4;
            end
           
            if  isfield (pulse, 'R')==0
                pulse.R=Rmin;
            end
            if pulse.R<Rmin
                pulse.R=Rmin;
            end
            kappa=(2.1*pulse.theta/pi);
            
            pulse.f_max=pulse.R/(2*pulse.Tp);
            pulse.B1max=kappa*1/gamma*pulse.theta*pulse.beta^(1/(2*pulse.n))*pulse.f_max*2/sqrt(pulse.R);
            pulse.t=linspace(0, pulse.Tp ,pulse.npts);

            %_____________________________________________________
    case 'Fixed_B1'
            if isfield(pulse, 'B1max')==0
                pulse.B1max=5e-4;
            end
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            if isfield(pulse, 'K')==0
                pulse.K=0.4;
            end     

            %xi can be changed to improve fidelity.
            if pulse.theta==pi/2
                xi=1.2;
                if pulse.K<1
                    pulse.K=1;
                end                
            else               
                xi=2.6;
               if pulse.K<0.4
                    pulse.K=0.4;
                end
            end 

            pulse.f_max=pulse.K*pulse.B1max*gamma;
            pulse.Tp=xi*2*pulse.f_max*(pulse.theta/(gamma*pulse.B1max*pulse.beta^(-1/(2.1*pulse.n))))^2;
            pulse.t=linspace(0, pulse.Tp ,pulse.npts);
            pulse.R=pulse.Tp*pulse.f_max*2;

            %____________________________________________-
        case 'Invariant_Tp'
            if isfield(pulse, 'Tp')==0
                pulse.Tp=200e-9;
            end
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            
            pulse.K=0.2;
            if pulse.theta==pi
                if pulse.n<5
                    Rmin=15;
                else
                    Rmin=1.91*pulse.n*2.33;
                end
                if isfield(pulse, 'R')==0
                    pulse.R=Rmin;
                end
                if pulse.R<Rmin
                    disp('R is too small. theta errors may occur')
                end
            elseif pulse.theta==pi/2
                pulse.R=1.85*pulse.n^-1.623+0.69;
            else
            end
            
            pulse.f_max=pulse.R/(pulse.Tp*2);
            pulse.B1max=pulse.f_max/(gamma*pulse.K);
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
            eta0=(-0.0002635*pulse.n.^2+0.3812*pulse.n+0.05669)/(pulse.n+1.806);
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




