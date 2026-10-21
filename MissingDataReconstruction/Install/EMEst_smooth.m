function obs = EMEst_smooth(Zt,Ht,F,alpha0,P0,R,Q,iters,direction)

% Description:
%   Expectation Maximization(EM) Algorithm with Kriged Kalman
%   Smooth to smooth observed data and interpolate missing data
%
% Input:
%   Zt - observed values, a n*m matrix, n is the length of 
%        observed time, m is the length of observed site.
%   Ht - Spatial Filed of Kriged Kalman Filter or a unit matrix
%   F - state transition matrix
%   alpha0 - initial state vector of Kriged Kalman Filter
%   P0 - covariance matrix of initial state vector
%   R - observation noise covariance matrix
%   Q - system state noise covariance matrix
%   iters - EM Algorithm iteration number
%
% Output:
%   obs - filtered and interpolated missing data

[n,m] = size(Zt);%size of observation matrix
p = size(Ht,2);%dimension of spatial filed
if direction == 1
    str_direciton = 'EM in N direciton:  ';
elseif direction == 2
    str_direciton = 'EM in E direciton:  ';
elseif direction == 3
    str_direciton = 'EM in U direciton:  ';
else
    error('Error in the direction!');
end
obs = zeros(n,m);
for j = 1:iters
    alpha = zeros(p,n+1);
    P = zeros(p,p,n+1);
    alpha(:,1) = alpha0;%initial state vector
    P(:,:,1) = P0;%covariance matrix of initial state vector
    
    %Kalman forward filter
    str_process = sprintf('the %dth Kalman forward Filter\n totally %d/%d',j,j,iters);
    str = [str_direciton,str_process];
    h1 = waitbar(0,str);
    
    alpha_ = zeros(p,n);
    P_ = zeros(p,p,n);
    for i=1:n
        if ~ishandle(h1)
            obs = [];
            return;
        end
        waitbar(i/n,h1);
        
        %change missing value and design matrix into 0
        II = any(isnan(Zt(i,:)'),2);
        Z = Zt(i,:)';  Z(II,:)=0;
        H = Ht;   H(II,:)=0;
        
        Rc = R;
        II_F = ~II;
        II = double(II);
        II_F = double(II_F);
        idex = II'*II + II_F'*II_F;
        idex = ~(logical(idex));
        Rc(idex) = 0;
        
        %Kalman Filter
        alpha_(:,i) = F*alpha(:,i);%one-step forecast state value
        P_(:,:,i) = F*P(:,:,i)*F' + Q;%covariance matrix of one-step forecast value
        K = P_(:,:,i)*H'/(H*P_(:,:,i)*H'+Rc);%gain matrix
        alpha(:,i+1) = alpha_(:,i) + K*(Z - H*alpha_(:,i));%filtered state value
        P(:,:,i+1) = P_(:,:,i) - K*H*P_(:,:,i);%covariance matrix of filtered state value
    end
    close(h1);
    
    %Kalman backward Smooth
    str_process = sprintf('the %dth Kalman backward Smooth\n totally %d/%d',j,j,iters);
    str = [str_direciton,str_process];
    h2 = waitbar(0,str);
    
    alpha_nt = zeros(p,n+1);
    alpha_nt(:,n+1) = alpha(:,n+1);
    obs(n,:) = (Ht*alpha_nt(:,n+1))';
    P_nt = zeros(p,p,n+1);
    P_nt(:,:,n+1) = P(:,:,n+1);
    J = zeros(p,p,n);%0~t-1
    for i=n:-1:1
        if ~ishandle(h2)
            obs = [];
            return;
        end
        waitbar((n-i+1)/n,h2);
        
        J(:,:,i) = P(:,:,i)*F'/(P_(:,:,i));
        alpha_nt(:,i) = alpha(:,i)+J(:,:,i)*(alpha_nt(:,i+1) - F*alpha(:,i));
        P_nt(:,:,i) = P(:,:,i)+J(:,:,i)*(P_nt(:,:,i+1) - P_(:,:,i))*J(:,:,i)';
        if i>1
            obs(i-1,:) = (Ht*alpha_nt(:,i))';
        end
    end
    close(h2);
    
    P_n_t = zeros(p,p,n);
    P_n_t(:,:,n) = (eye(p)-K*H)*F*P(:,:,n);
    for i=n:-1:2
        P_n_t(:,:,i-1) = P(:,:,i)*J(:,:,i-1)'+...
            J(:,:,i)*(P_n_t(:,:,i)-F*P(:,:,i))*J(:,:,i-1)';
    end
    
    %calculate the new parameter of Kalman filter
    str_process = sprintf('the %dth update parameter of Kalman filter\n totally %d/%d',j,j,iters);
    str = [str_direciton,str_process];
    h3 = waitbar(0,str);
    
    A = zeros(p,p);
    B = zeros(p,p);
    C = zeros(p,p);
    Rt = zeros(m,m);
    for i=1:n
        if ~ishandle(h3)
            obs = [];
            return;
        end
        waitbar(i/n,h3);
        
        A = A + P_nt(:,:,i) + alpha_nt(:,i)*alpha_nt(:,i)';
        B = B + P_n_t(:,:,i) + alpha_nt(:,i+1)*alpha_nt(:,i)';
        C = C + P_nt(:,:,i+1) + alpha_nt(:,i+1)*alpha_nt(:,i+1)';
        
        II = any(isnan(Zt(i,:)'),2);
        Z = Zt(i,:)';  Z(II,:)=0;
        H = Ht;   H(II,:)=0;
        
        %change R into new form
        Rc = zeros(m,m); 
        II = double(II); II = II*II'; II = logical(II);
        Rc(II) = R(II);
        
        Rt = Rt +...
            ((Z-H*alpha_nt(:,i))*(Z-H*alpha_nt(:,i))'+...
            H*P_nt(:,:,i)*H' + Rc)/n;
    end
    close(h3);
    F = B/A;
    Q = (C - F*B')/n;
    R = Rt;
	
    alpha0 = alpha_nt(:,1);
    P0 = P_nt(:,:,1);
end
