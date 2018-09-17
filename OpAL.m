function latents = OpAL(data, x)
    
    % Opponent Actor Learning model (Collins & Frank 2014)
    %
    % USAGE: latents = OpAL(data, x)
    %
    % INPUTS:
    %   data - single subject data
    %   x (optional) - parameters
  
    if nargin < 2
        alpha = 0.1;
        a = 2;
        b = 2;
    else
        alpha = x(1);
        a = x(2);
        b = x(3);
    end
   
    % same as AU and ACU TODO look into it
    G_0R = 0.5 * var_to_S(alpha, beta, 16);
    N_0R = 0.5 * var_to_S(alpha, beta, 16); 
    G_0S = 0.5 * var_to_S(alpha, beta, 0.00001);
    N_0S = 0.5 * var_to_S(alpha, beta, 0.00001);
    V_0 = 0;

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
        
        % update, if choice was made
        if ~data.timeout(n)
            V = V + alpha * (r - V);
            G(c) = G(c) + alpha * G(c) * (r - V);
            N(c) = N(c) - alpha * N(c) * (r - V);
        end
    end
end

function S = var_to_S(alpha, beta, var)
    S = alpha / beta * sqrt(2 / pi) * sqrt(var);
end
