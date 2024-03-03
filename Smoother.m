%% Regolarizzazione di Rauch, Tung, Striebel 

% in avanti = EKF
% all'indietro:
P_k1_k = out.P_k1_k;
P_k_k = out.P_k_k;
csi_k_k = out.csi_k_k;
csi_k1_k = out.csi_k1_k;
F_k = out.F_k;
csi_k_n = zeros(6,6001);
csi_k_n(:,end) = csi_k_k(:,:,end);
P_k_n = zeros(6,6,6001);
P_k_n(:,:,end) = P_k_k(:,:,end);

for k = 6000:-1:1
    Ck = P_k_k(:,:,k) * (F_k(:,:,k+1)') / (P_k1_k(:,:,k+1));
    csi_k_n(:,k) = csi_k_k(:,:,k) + Ck * (csi_k_n(:,k+1) - csi_k1_k(:,:,k+1));
    P_k_n(:,:,k) = P_k_k(:,:,k) + Ck * (P_k_n(:,:,k+1) - P_k1_k(:,:,k+1)) * (Ck');
end


x_smooth = csi_k_n(1,:);
z_smooth = csi_k_n(2,:);
theta_smooth = csi_k_n(3,:);
u_smooth = csi_k_n(4,:);
w_smooth = csi_k_n(5,:);
q_smooth = csi_k_n(6,:);
