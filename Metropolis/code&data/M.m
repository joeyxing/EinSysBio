function z = M(theta)
    timespan = [0 30];
    y0 = [0, 2, 0, 1, 0, 3];
    diffeq = @(t,y) model(t,y,theta);
    [t, y] = ode45(diffeq, timespan,y0);
    
    t0 = 0:.01:30;
    yz = interp1(t,y,t0);
    % plot(t0, yz)
    % legend('m1','p1','m2','p2','m3','p3');
    z = zeros(7,3);
    for n = 1:7
        % tx = find(t0==4*n);
        % checkpoints(n,1) = 4*n;
        for m = 1:3
            z(n,m) = yz(t0==4*n,2*m-1);
        end
    end
end
