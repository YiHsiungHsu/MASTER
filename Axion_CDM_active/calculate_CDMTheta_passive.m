function dy = calculate_CDMTheta_passive(u_m,Theta)
     dy = zeros(2,1); % a column vector
     dy(1) = Theta(2);
     dy(2) = -3/2/u_m*Theta(2)-sin(Theta(1));