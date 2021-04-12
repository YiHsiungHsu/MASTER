function dy = calculate_ByTheta_passive(u_m,Theta,mu)
     dy = zeros(4,1); % a column vector
     dy(1) = Theta(2);
     dy(2) = -3/2/u_m*Theta(2)-sin(Theta(1));
     dy(3) = Theta(4);
     dy(4) = -3/2/u_m*Theta(4)-mu^2*sin(Theta(3));