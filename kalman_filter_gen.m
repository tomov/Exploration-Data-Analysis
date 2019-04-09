function [data, latents] = kalman_filter(data, w)
   
    % same as kalman_filter but run generatively
    
    vars = [0.00001 16];   % reward variances for safe and risky options
    N = length(data.block);
    tau2 = zeros(N,2)+vars(1);
    tau2(data.cond==1|data.cond==3,1) = vars(2);
    tau2(data.cond==2|data.cond==3,2) = vars(2);
    
    for n = 1:N
        
        % initialization at the start of each block
        if n == 1 || data.block(n)~=data.block(n-1)
            m = [0 0];      % posterior mean
            s = [100 100];  % posterior variance
        end

        V = m(1) - m(2);
        RU = sqrt(s(1)) - sqrt(s(2));
        TU = sqrt(s(1) + s(2));

        DV = w(1) * V + w(2) * RU + w(3) * V/TU;
        p = normcdf(DV); % P(a == 1)
        c = binornd(1, p);
        if c == 0
            c = 2; % choice is 1 or 2 here
        end

        if c == 1
            r = data.r1(n);
        else
            assert(c == 2);
            r = data.r2(n);
        end
       
        data.choice(n) = c;
        data.reward(n) = r;
        %c = data.choice(n);
        %r = data.reward(n);
        
        % store latents
        latents.m(n,:) = m;
        latents.s(n,:) = s;
        
        % update, if choice was made
        if ~data.timeout(n)
            k = s(c)/(s(c)+tau2(n,c));    % Kalman gain
            err = r - m(c);            % prediction error
            m(c) = m(c) + k*err;       % posterior mean
            s(c) = s(c) - k*s(c);      % posterior variance
        end
        
    end
