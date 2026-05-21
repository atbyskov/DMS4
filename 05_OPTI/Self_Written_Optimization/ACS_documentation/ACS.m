function T = ACS(IterNo,A,T)

% Function for ACS - Adaptive Constraint Scaling proposed by Le et al (2010).
% This implementation follows the algorithm suggested in
% Oest, J. and Lund, E. (2017): Topology optimization with finite-life fatigue constraints. 
% Structural and Multidisciplinary Optimization 56(5), 1045-1059.
%
% Input:
%   IterNo: Iteration number
%   T: structure containing aggregation info including
%      - T.ng: number of aggegated constraint functions, j=1,..,T.ng
%   A: structure containing analysis data including
%      - A.g(j)    : the aggregated value, j=1,..,T.ng from the previous iteration
%      - A.g_max(j): the maximum value of all included function values in g(j), j=1,..,T.ng from the previous iteration
% Output:
%   T.c_Le(j)     : The c parameter in Le et al (2010)
%   T.alpha_Le(j) : The alpha parameter in Le et al (2010)
%   T.c_Le_plot   : The c parameter for all aggregated functions for all iterations are stored for plotting
%   T.a_Le_plot   : The alpha parameter for all aggregated functions for all iterations are stored for plotting

% If IterNo = 0 we choose alpha = 1 and the scaling T.c_Le(j) is computed directly
if IterNo == 0
  for j = 1:T.ng
    T.c_Le(j) = A.g_max(j)/A.g(j); % We need to use values from this iteration as it is the first iteration
  end
  % Initialize alpha parameter to 1
  T.alpha_Le = ones(T.ng,1);

% If IterNo = 1 we still choose alpha = 1 and store old value in T.c_Le_old
elseif IterNo == 1
  T.c_Le_old = T.c_Le; % We store the old vector of c scaling vectors
  for j = 1:T.ng
    T.c_Le(j) = A.g_max(j)/A.g(j); % A.g_max(j) and A.g(j) are from the previous iteration
  end
  % Store T.c_Le_plot and T.a_Le_plot for plotting
  T.c_Le_plot(IterNo,:) = T.c_Le;
  T.a_Le_plot(IterNo,:) = ones(T.ng,1);

% If IterNo > 1 we can check for oscillations, as we have three values
elseif IterNo > 1
        
  % First, store c from previous iterations in T.c_Le_old and T.c_Le_oldold :-)
  T.c_Le_oldold   = T.c_Le_old; % We store the previous old vector of c scaling vectors
  T.c_Le_old      = T.c_Le;     % We store the old vector of c scaling vectors
        
  for j = 1:T.ng
    % Calculate c using current alpha value - again: A.g_max(j) and A.g(j) are from the previous iteration
    T.c_Le(j) = T.alpha_Le(j) * A.g_max(j)/A.g(j) + (1-T.alpha_Le(j)) * T.c_Le_old(j);

    % Do we oscillate?
    % We don't want to check oscillation on too small numbers
    % Thus, the parameter decimal decides the number of digits included in
    % the c parameters when checking (default decimal=2)
    
    % Scaling factors
    scale_fac_up = 1.20;
    scale_fac_down = 0.80;
    decimal = 2;


    if (round(T.c_Le_oldold(j),decimal) < round(T.c_Le_old(j),decimal) && ...
        round(T.c_Le_old(j),decimal)    > round(T.c_Le(j),decimal))    || ...
       (round(T.c_Le_oldold(j),decimal) > round(T.c_Le_old(j),decimal) && ...
        round(T.c_Le_old(j),decimal)    < round(T.c_Le(j),decimal))    ...
      % decrease alpha 
      T.alpha_Le(j) = max(0.5,T.alpha_Le(j)*scale_fac_down);
    else  % If its stable, then increase alpha
      % increase slowly
      T.alpha_Le(j) = min(1,T.alpha_Le(j)*scale_fac_up);
    end
            
    % Recalculate c using new (or maybe unchanged) alpha
    T.c_Le(j) = T.alpha_Le(j) * A.g_max(j)/A.g(j) + (1-T.alpha_Le(j)) * T.c_Le_old(j);
    
    % Testing: if limitation on the change of c is active (the logical variable T.c_Le_control 
    % is set true), then enforce a maximum change on T.c_Le(j) to be T.c_Le_control_delta (perhaps 0.02)
    % if abs(T.c_Le(j)-T.c_Le_old(j)) >= T.c_Le_control_delta && T.c_Le_control
    %   if T.c_Le > T.c_Le_old
    %     T.c_Le = T.c_Le_old + T.c_Le_control_delta;
    %   else
    %     T.c_Le = T.c_Le_old - T.c_Le_control_delta;
    %   end
    % end  

    % Store T.c_Le_plot and T.a_Le_plot for plotting
    T.c_Le_plot(IterNo,:) = T.c_Le;
    T.a_Le_plot(IterNo,:) = T.alpha_Le;
  end
end