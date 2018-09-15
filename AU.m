function latents = AU(data, x)
    
    % Actor with uncertainty (AU) model
    %
    % USAGE: latents = AU(data, x)
    %
    % INPUTS:
    %   data - single subject data
    %   x (optional) - parameters
  
    if nargin < 2
        alpha = 0.1;
        beta = 0.1;
        a = 2;
        b = 2;
    else
        alpha = x(1);
        beta = x(2);
        a = x(3);
        b = x(4);
    end

    % G + N = S = irreducible variance = 16 for risky, 0.00001 for safe (see Sam's kalman_filter.m)
    % G - N = Q = expected r = 0 for both
    G_0R = 8;
    N_0R = 8;
    G_0S = 0.000005;
    N_0S = 0.000005;

    N = length(data.block);

    latents.a = a;
    latents.b = b;
    
    for n = 1:N
        
        % initialization at the start of each block
        if n == 1 || data.block(n)~=data.block(n-1)
            switch data.cond(n)
                case 1 % RS
                    G = [G_0R G_0S];
                    N = [N_0R N_0S];
                case 2 % SR
                    G = [G_0S G_0R];
                    N = [N_0S N_0R];
                case 3 % RR 
                    G = [G_0R G_0R];
                    N = [N_0R N_0R];
                case 4 % SS
                    G = [G_0S G_0S];
                    N = [N_0S N_0S];
            end
        end

        loglik = a * G - b * N - logsumexp(a * G - b * N); % avoid numeric underflow

        c = data.choice(n);
        r = data.reward(n);

        % store latents
        latents.G(n,:) = G;
        latents.N(n,:) = N;
        latents.loglik(n,:) = loglik;
        
        % update, if choice was made
        if ~data.timeout(n)
            Q = G(c) - N(c);
            G(c) = G(c) + alpha * rect_pos(r - Q) - beta * G(c);
            N(c) = N(c) + alpha * rect_neg(r - Q) - beta * N(c);
        end
    end
end


function x = rect_pos(x)
    x = max(x,0);
end

function x = rect_neg(x)
    x = max(-x,0);
end
