function [ Phi ] =UnwrapPhase(x, phi)
%Unwraphs the phase of phi
%
%Syntax
%[ Phi ] =UnwrapPhase(detuning, phi)
[~, ii]=min(abs(x));
phiRes=phi(ii);

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
Phi=Phi-Phi(ii)+phiRes;

end

% function [ Phi ] =UnwrapPhase(x, phi)
% %Unwraps Phase a from cart2sph(Mx, My)
% 
% [~,ii]=min(abs(x));
% phi0=phi(ii);
% 
% phi=unwrap(phi);
% phi=phi-min(phi);
% figure
% plot(phi, 'g')
% hold on
% 
% diff(1)=0;
% ind(1)=1;
%  for jj=2:length(phi)
%      diff(jj)=phi(jj)-phi(jj-1);
%       if abs(diff(jj))>0.95
%             ind(jj)=ind(jj-1)*-1*sign(diff(jj));
%         else
%             ind(jj)=ind(jj-1);
%         end
%  end
% ind=1-(ind+1)/2;
% Phi=phi+pi*ind';
% plot(Phi, 'y')
% 
% Phi=unwrap(Phi);
% Phi=Phi-Phi(ii)+phi0;
% 
% plot(Phi, 'r')
% plot(diff, 'b')
% 
% end

