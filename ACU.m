function latents = ACU(data, x)
    
    % Actor-critic with uncertainty (ACU) model
    %
    % USAGE: latents = ACU(data, x)
    %
    % INPUTS:
    %   data - single subject data
    %   x (optional) - parameters
  
    if nargin < 2
        G_0R = 0.1;
        N_0R = 0.1;
        G_0S = 0.1;
        N_0S = 0.1;
        V_0 = 0;
        alpha = 0.1;
        beta = 0.1;
        a = 2;
        b = 2;
    else
        G_0R = x(1);
        N_0R = x(2);
        G_0S = x(3);
        N_0S = x(4);
        V_0 = x(5);
        alpha = x(6);
        beta = x(7);
        a = x(8);
        b = x(9);
    end

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
        
        % update
        V = V + alpha * (r - V);
        G(c) = G(c) + alpha * rect_pos(r - V) - beta * G(c);
        N(c) = N(c) + alpha * rect_neg(r - V) - beta * N(c);
    end
end


function x = rect_pos(x)
    x = max(x,0);
end

function x = rect_neg(x)
    x = max(-x,0);
end
