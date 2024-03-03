function [Particles_out, P_k_k, csi_k1_k1, e_k, w_posteriori_norm] = 
Particle_Filter(is_particle, Particles_in, tau_u, tau_q, y_alpha_m, z_m,...
    theta_m, t,z_boa_2, R, Q, S, csi_0, dt)
      
    m = S(1);
    I_y = S(2);
    Damping_matrix = diag([S(3),S(4),S(5)]);
    
    N = 20000;                      % Numero di particelle che generiamo 
    
    if(is_particle)
        w_posteriori = zeros(1,N);
        if(t == 0)               % Generazione particelle nell'istante t=0
    
            % Definiamo i range di incertezza
            
            range_x = 3;
            range_z = 3;
            range_theta = 0.01;
            range_u = 0.01;
            range_w = 0.01;
            range_q = 0.01;
    
            % Generiamo le particelle per ogni variabile del problema
    
            particles_x = rand(1,N)*range_x + (csi_0(1)-range_x/2);
            particles_z = rand(1,N)*range_z + (csi_0(2)-range_z/2);
            particles_theta = rand(1,N)*range_theta + 
                                (csi_0(3)-range_theta/2);
            particles_u = rand(1,N)*range_u + (csi_0(4)-range_u/2);
            particles_w = rand(1,N)*range_w + (csi_0(5)-range_w/2);
            particles_q = rand(1,N)*range_q + (csi_0(6)-range_q/2);
            Particles_k1_k = [particles_x; particles_z; particles_theta;
                particles_u; particles_w; particles_q];
            w_posteriori = ones(1,N)/N;
        
        else
        
            % Definizione del vettori vuoti da usare nell'algoritmo
            Particles_k1_k = zeros(size(Particles_in));
            e = zeros(3,N);
            w_posteriori = zeros(1,N);
            e_y_alpha_m = zeros(1,N);
            e_z_m = zeros(1,N);
            e_theta_m = zeros(1,N);
            is_w_zero = true;
            Particles_k_k = Particles_in;
            
            for i=1:N                % Per ogni i-esima particella viene 
                                     % fatta predizione e correzione
                x = Particles_k_k(1,i);
                z = Particles_k_k(2,i);
                theta = wrapTo2Pi(Particles_k_k(3,i));
                u = Particles_k_k(4,i);
                w = Particles_k_k(5,i);
                q = Particles_k_k(6,i);
                tau_du = randn()*Q(1,1);
                tau_dw = randn()*Q(2,2);
                tau_dq = randn()*Q(3,3);
        
                csi_dot = [cos(theta)*u + sin(theta)*w
                           cos(theta)*w - sin(theta)*u
                           q
                           (tau_u + tau_du - Damping_matrix(1,1)*u ...
                           - m*q*w)/m
                           (tau_dw - Damping_matrix(2,2)*w ...
                           + m*q*u)/(m)
                           (tau_q + tau_dq - Damping_matrix(3,3)*q)/I_y];
                
                % Predizione
                Particles_k1_k(:,i) = Particles_k_k(:,i) + dt*csi_dot;
                
                e_y_alpha_m(i) = y_alpha_m ...
                            - atan2(z_boa_2 - Particles_k1_k(2,i), ...
                            Particles_k1_k(1,i));
                e_z_m(i) = z_m - Particles_k1_k(2,i);
                e_theta_m(i) = theta_m - Particles_k1_k(3,i);
                e(:,i) = [e_y_alpha_m(i) e_z_m(i) e_theta_m(i)]';
        
                % Correzione  % Pesi a posteriori
                w_posteriori(i)=double(exp(-1/2*(e(:,i)'*inv(R)*e(:,i))));    
            end
            if(sum(w_posteriori) == 0)
                % Nel caso in cui la verosomiglianza sia nulla 
                % (per troncamentro) poniamo i pesi come di pdf uniforme
                w_posteriori = ones(1,N)/N;

            end
        end
        % Normalizzazione dei pesi
        w_posteriori_norm = w_posteriori./sum(w_posteriori);

        % Calcolo della stima (EAP/MMSE)
        csi_k1_k1 = zeros(6,1);
        for k=1:N
            csi_k1_k1 = csi_k1_k1 ...
                        + w_posteriori_norm(k).*Particles_k1_k(:,k);
        end

        csi_k1_k1(3) = wrapto2Pi(csi_k1_k1(3));
        
        % Calcolo dell'innovazione
        e_y_alpha = y_alpha_m ...
                    - atan2(z_boa_2 - csi_k1_k1(2), csi_k1_k1(1));
        e_z = z_m - csi_k1_k1(2);
        e_theta = theta_m - csi_k1_k1(3);
    
        e_k = [e_y_alpha, e_z, e_theta]';

        % Resampling (Algoritmo della roulette)
        Resample = zeros(size(Particles_in));                
        c = zeros(1,N);
        u = zeros(1,N);
        c(1) = w_posteriori_norm(1);
        for i=2:N
            c(i) = c(i-1) + w_posteriori_norm(i);
        end
        u(1) = rand()/N;
        i = 1;
        for j=1:N
            while(u(j)>c(i))
                i=i+1;
            end
            Resample(:,j) = Particles_k1_k(:,i);
            if(j<N)
                u(j+1) = u(j) + 1/N;
            end
        end
        
        % Calcolo della varianza di stima

        P_k_k = zeros(6,6);

        for i=1:N
            
           P_k_k = P_k_k ...
            + w_posteriori_norm(i).*...
           (Particles_k1_k(:,i) - csi_k1_k1)*...
           (Particles_k1_k(:,i) - csi_k1_k1)';

        end
        
        % Particelle in uscita (Sporcate con un rumore bianco)

        Particles_out = diag([0.01 0.01 0.001 0 0 0])*randn(6,N) +...
                        Resample;
        
    else

        % Nel caso in cui non usiamo il Particle Filter andiamo a dare in
        % uscita segnali nulli
        
        Particles_out = zeros(6,N);
        csi_k1_k1 = zeros(6,1);
        e_k = zeros(3,1);
        w_posteriori_norm = zeros(1,N);
        P_k_k = zeros(6,6);

    end
end