%% Animazioni e grafici per l'analisi dei filtri

if(~is_particle)
    Smoother;
end

system_states = out.system_states';
esteemated_system_states = out.esteemated_system_states';
tf = 60;
e_k = out.e_k;
x_t = system_states(1,:);
z_t = system_states(2,:);

x_hat_t = esteemated_system_states(1,:);
z_hat_t = esteemated_system_states(2,:);

% Plot della traiettoria reale e della sua stima (regolarizzata e non)
figure
title("Traiettoria del sistema e dello stimatore");
grid on
hold on
xlabel("x [m]");
ylabel("z [m]");
plot(x_boa_1,-z_boa_1,'o');
plot(x_boa_2,-z_boa_2,'o');
plot(x_t,-z_t);
plot(x_hat_t,-z_hat_t);

if(is_particle)
    legend(["Boa 1", "Boa 2", "Traiettoria reale", "Traiettoria stimata"]);
else
    plot(x_smooth,-z_smooth);
    legend(["Boa 1", "Boa 2", "Traiettoria reale", ...
        "Traiettoria stimata", "Stima regolarizzata"]);
end
hold off

% Plot dell'andamento degli errori di stima delle variabili di stato
figure
tl = tiledlayout(2,3,'TileSpacing','Compact');

%Posizione x
nexttile
plot(0:dt:tf,system_states(1,:) - esteemated_system_states(1,:));
title("Errore di Posizione x");
xlabel("tempo t [s]");

%Posizione z
nexttile
plot(0:dt:tf,-(system_states(2,:) - esteemated_system_states(2,:)));
title("Errore di Posizione z");
xlabel("tempo t [s]");

%Angolo theta
nexttile
plot(0:dt:tf,(system_states(3,:) - esteemated_system_states(3,:)));
title("Errore sull'angolo \theta");
xlabel("tempo t [s]");

%Velocità u
nexttile
plot(0:dt:tf,system_states(4,:) - esteemated_system_states(4,:));
title("Errore sulla velocità longitudinale u");
xlabel("tempo t [s]");

%Velocità w
nexttile
plot(0:dt:tf,system_states(5,:) - esteemated_system_states(5,:));
title("Errore sulla velocità laterale w");
xlabel("tempo t [s]");

%Velocità angolare q
nexttile
plot(0:dt:tf,system_states(6,:) - esteemated_system_states(6,:));
title("Errore sulla velocità angolare q");
xlabel("tempo t [s]");

title(tl,"Errori di stima");



if(~is_particle)
    % Plot dell'andamento degli errori di stima delle variabili di stato
    % (con regolarizzazione di Rauch Tung Striebel)
    figure
    tl = tiledlayout(2,3,'TileSpacing','Compact');
    
    %Posizione x
    nexttile
    plot(0:dt:tf,system_states(1,:) - x_smooth);
    title("Errore di Posizione x");
    xlabel("tempo t [s]");
    
    %Posizione z
    nexttile
    plot(0:dt:tf,-(system_states(2,:) - z_smooth));
    title("Errore di Posizione z");
    xlabel("tempo t [s]");
    
    %Angolo theta
    nexttile
    plot(0:dt:tf,(system_states(3,:) - theta_smooth));
    title("Errore sull'angolo \theta");
    xlabel("tempo t [s]");
    
    %Velocità u
    nexttile
    plot(0:dt:tf,system_states(4,:) - u_smooth);
    title("Errore sulla velocità longitudinale u");
    xlabel("tempo t [s]");
    
    %Velocità w
    nexttile
    plot(0:dt:tf,system_states(5,:) - w_smooth);
    title("Errore sulla velocità laterale w");
    xlabel("tempo t [s]");
    
    %Velocità angolare q
    nexttile
    plot(0:dt:tf,system_states(6,:) - q_smooth);
    title("Errore sulla velocità angolare q");
    xlabel("tempo t [s]");
    
    title(tl,"Errori di stima (Rauch Tung Striebel smoothing)");
end

%Plot delle componenti dell'innovazione e delle loro funzioni di
%autocovarianza
figure
tl_2 = tiledlayout(3,2,"TileSpacing","compact");

nexttile                  
plot(0:dt:tf,e_k(1,:));                % Prima componente dell'innovazione
title("Prima componente dell'innovazione e(1)");
xlabel("tempo t [s]");

nexttile                
% Autocovarianza della prima componente dell'innovazione
r1 = xcorr(e_k(1,:),e_k(1,:)); 
plot(0:dt:tf*2,r1);
title("Funzione di autocovarianza di e(1)");
xlabel("distanza temporale \Deltat [s]");

nexttile
plot(0:dt:tf,e_k(2,:));               % Seconda componente dell'innovazione
title("Seconda componente dell'innovazione e(2)");
xlabel("tempo t [s]");

nexttile
% Autocovarianza della seconda componente dell'innovazione
r2 = xcorr(e_k(2,:),e_k(2,:)); 
plot(0:dt:tf*2,r2);
title("Funzione di autocovarianza di e(2)");
xlabel("distanza temporale \Deltat [s]");

nexttile
plot(0:dt:tf,(e_k(3,:)));                 % Terza componente dell'innovazione
title("Terza componente dell'innovazione e(3)");
xlabel("tempo t [s]");

nexttile
% Autocovarianza della terza componente dell'innovazione
r3 = xcorr((e_k(3,:)),(e_k(3,:)));    
plot(0:dt:tf*2,r3);
title("Funzione di autocovarianza di e(3)");
xlabel("distanza temporale \Deltat [s]");

title(tl_2,"Innovazione");

% Animazione dell'andamento temporale della traiettoria e della sua stima
figure
hold on
grid on
title("Traiettoria del sistema e traiettoria stimata");
xlabel("x [m]");
ylabel("z [m]");
plot(x_boa_1,-z_boa_1,'o');
plot(x_boa_2,-z_boa_2,'o');
a_states = animatedline('Color',[0 0 1]);
a_esteem = animatedline('Color',[1 0 0]);
legend(["Boa 1", "Boa 2", "Traiettoria reale", "Traiettoria stimata"])
n = length(x_t);

%Ciclo for per la creazione dei frame
for t = 1:n
    % traiettoria del sistema
    addpoints(a_states,x_t(t),-z_t(t));

    % traiettoria della stima
    addpoints(a_esteem,x_hat_t(t),-z_hat_t(t));

    % update screen
    drawnow limitrate
    pause(0.0001)
end
hold off