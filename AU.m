function latents = AU(data, x)
    
    % Actor with uncertainty (AU) model
    %
    % USAGE: latents = AU(data, x)
    %
    % INPUTS:
    %   data - single subject data
    %   x (optional) - parameters
  
    if nargin < 2
        G_0 = [0.1 0.1];
        N_0 = [0.1 0.1];
        alpha = 0.1;
        beta = 0.1;
        a = 2;
        b = 2;
    else
        G_0 = [x(1) x(1)];
        N_0 = [x(2) x(2)];
        alpha = x(3);
        beta = x(4);
        a = x(5);
        b = x(6);
    end

    N = length(data.block);
    
    for n = 1:N
        
        % initialization at the start of each block
        if n == 1 || data.block(n)~=data.block(n-1)
            G = G_0;
            N = N_0;
        end

        loglik = a * G - b * N - logsumexp(a * G - b * N); % avoid numeric underflow

        c = data.choice(n);
        r = data.reward(n);

        % store latents
        latents.G(n,:) = G;
        latents.N(n,:) = N;
        latents.loglik(n,:) = loglik;
        
        % update
        Q = G(c) - N(c);
        G(c) = G(c) + alpha * rect_pos(r - Q) - beta * G(c);
        N(c) = N(c) + alpha * rect_neg(r - Q) - beta * N(c);
    end
end


function x = rect_pos(x)
    if x <= 0
        x = 0;
    end
end

function x = rect_neg(x)
    if x >= 0
        x = 0;
    end
end
