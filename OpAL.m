function latents = OpAL(data, x)
    
    % Opponent Actor Learning model (Collins & Frank 2014)
    %
    % USAGE: latents = OpAL(data, x)
    %
    % INPUTS:
    %   data - single subject data
    %   x (optional) - parameters
  
    if nargin < 2
        G_0 = [0.1 0.1];
        N_0 = [0.1 0.1];
        V_0 = 0.1;
        alpha = 0.1;
        a = 2;
        b = 2;
    else
        G_0 = [x(1) x(1)];
        N_0 = [x(2) x(2)];
        V_0 = x(3);
        alpha = x(4);
        a = x(5);
        b = x(6);
    end

    N = length(data.block);
    
    for n = 1:N
        
        % initialization at the start of each block
        if n == 1 || data.block(n)~=data.block(n-1)
            G = G_0;
            N = N_0;
            V = V_0;
        end

        loglik = a * G - b * N - logsumexp(a * G - b * N); % avoid numeric underflow

        c = data.choice(n);
        r = data.reward(n);

        % store latents
        latents.G(n,:) = G;
        latents.N(n,:) = N;
        latents.V(n,:) = V;
        latents.loglik(n,:) = loglik;
        
        % update
        V = V + alpha * (r - V);
        G(c) = G(c) + alpha * G(c) * (r - V);
        N(c) = N(c) - alpha * N(c) * (r - V);
    end
end

