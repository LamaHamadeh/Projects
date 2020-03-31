


% %%%Damped Harmonic Oscillator
% %-----------------------------
% function rhs = Oscillators_rhs(t,y0,dummy,d,w)
%         %set the initial value for the position
%         y1 = y0(1);
%         %set the initial value for the velocity
%         y2 = y0(2);
%         %compute the equation
%         rhs = [y2 ; -w^2*y1-d*y2];
% 
% end
% 

% %%%Van Der Pol Equation
% %----------------------
% function rhs = Oscillators_rhs(tspan,y0,dummy,mu)
%         %set the initial value for the position
%         y1 = y0(1);
%         %set the initial value for the velocity
%         y2 = y0(2);
%         %compute the equation
%         rhs = [y2 ; (mu-mu*y1^2)*y2-y1];
% 
% end

