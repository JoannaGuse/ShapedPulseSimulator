function pulse = Def_TanhTan(pulse)

pulse.class='FM';
pulse.descr='Tanh-tan adiabatic pulse. AM and FM. Best suited for Inversion';
pulse.DynamicVars={'n','beta', 'theta'};
pulse.SupportedMethods={'None', 'Fixed_B1', 'Fixed_Tp', 'None'};

%set up defaults
if isfield(pulse, 'theta')==0
    pulse.theta=pi;
end
if isfield(pulse, 'n')==0
    pulse.n=10;
end
if isfield(pulse, 'beta')==0
    pulse.beta=atan(18);
end

% add F1 functions
pulse.F1=F1_def(pulse);
pulse.F2=F2_def(pulse);
pulse.RelArea=trapz(linspace(0, 1, pulse.npts),pulse.F1);

%_________ Add fucntion handles to do operation
pulse.SimMethod=@(pulse)SimMethod(pulse);

end
% %________________________________________
% function F1=F1_def(pulse)
% tau=linspace(0, 1, round(pulse.npts/2));
% 
% A=tanh(pulse.n*tau);
% F1=[A, fliplr(A)];
% 
% if length(F1)~=pulse.npts
%     F1=[A(1:end-1), fliplr(A)];
% end
% end
% 
% %________________________________________
% function F2=F2_def(pulse)
% tau=linspace(-1, 1, pulse.npts);
% F2=tan(pulse.beta*tau)./tan(pulse.beta);
% 
% end

%________________________________________
function F1=F1_def(pulse)
tau=linspace(0, 1, round(pulse.npts/2));

A=tanh(pulse.n*tau);
F1=[A, fliplr(A)];

if length(F1)~=pulse.npts
    F1=[A(1:end-1), fliplr(A)];
end
end

%________________________________________
function F2=F2_def(pulse)
tau=linspace(-1, 1, pulse.npts);
F2=tan(pulse.beta*tau)./tan(pulse.beta);

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
            %_________________________________________________
        case 'Fixed_B1'
            if isfield(pulse, 'B1max')==0
                pulse.B1max=5e-4;
            end
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            
            pulse.R=30;
            pulse.K=1.6;
            
            pulse.f_max=(gamma*pulse.B1max)*pulse.K;
            pulse.Tp=pulse.R/(pulse.f_max*2);
            pulse.t=linspace(0, pulse.Tp ,pulse.npts);
         
            case 'Fixed_Tp'
            %check defaults
            if isfield(pulse, 'theta')==0
                pulse.theta=pi;
            end
            if  isfield (pulse, 'Tp')==0
                pulse.Tp=200e-9;
            end

            pulse.R=30;
            pulse.K=1.6;
            
            pulse.f_max=pulse.R/(2*pulse.Tp);
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
%         %___________________________________________________
%         case 'Fixed_opBW'
%             %check defaults
%             if isfield(pulse, 'theta')==0
%                 pulse.theta=pi;
%             end
%             if  isfield (pulse, 'opBW')==0
%                 pulse.opBW=40e6;
%             end
%             
%             pulse.R=0.3;
%             eta0=(0.351*pulse.n^2+3.481*pulse.n+2.796)/(pulse.n^2+13.08*pulse.n+19.53);
%             if pulse.theta==pi
%                 pulse.eta=eta0;
%             elseif pulse.theta==pi/2
%                 pulse.eta=0.7166*eta0;
%             end
%             pulse.f_max=pulse.eta*pulse.opBW/2;
%             pulse.Tp=pulse.R/(2*pulse.f_max);
%             pulse.t=linspace(0, pulse.Tp, pulse.npts);
%             pulse.B1max=pulse.theta/(gamma*pulse.Tp*pulse.RelArea);
%             
            
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

