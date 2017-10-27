function dydt = model(t, y, theta)
    % 
    m1 = y(1);
    p1 = y(2);
    m2 = y(3);
    p2 = y(4);
    m3 = y(5);
    p3 = y(6);
    
    % theta = [alpha0, n, beta, alpha]
    alpha0 = theta(1);
    n      = theta(2);
    beta   = theta(3);
    alpha  = theta(4);
    
    % define differential equation
    dm1_dt = alpha0 + alpha/(1+p3^n) - m1;
    dp1_dt = beta*(m1 - p1);
    dm2_dt = alpha0 + alpha/(1+p1^n) - m2;
    dp2_dt = beta*(m2 - p2);
    dm3_dt = alpha0 + alpha/(1+p2^n) - m3;
    dp3_dt = beta*(m3 - p3);
    
    % return vector
    dydt = [dm1_dt; dp1_dt; dm2_dt; dp2_dt; dm3_dt; dp3_dt];
end
