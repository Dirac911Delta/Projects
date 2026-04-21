function post = post_processing(material, settings, results)
% post_processing.m
% Generates derived metrics and visualizations.

r = results.r;
z = results.z;
t = results.t;
T = results.T;

nt = numel(t);
peakT = zeros(nt,1);
for k = 1:nt
    peakT(k) = max(T(:,:,k), [], 'all');
end

% Radial gradient at mid-plane vs time (max absolute)
[~, jmid] = min(abs(z - z(end)/2));
radial_grad_max = zeros(nt,1);
for k = 1:nt
    Tr = T(:,jmid,k);
    dTdr = gradient(Tr, r);
    radial_grad_max(k) = max(abs(dTdr));
end

% Final-time axial profile on axis
axial_profile_final = squeeze(T(1,:,end));

% Optional thermal lens proxy: dn/dT * max radial gradient (final time)
if isnan(material.dn_dT)
    lens_proxy = NaN;
else
    dTdr_final = gradient(T(:,jmid,end), r);
    lens_proxy = material.dn_dT * max(abs(dTdr_final));
end

% Plots
figure('Name','Temperature contour (final time)','Color','w');
contourf(r*1e3, z*1e3, T(:,:,end).', 40, 'LineColor', 'none');
colorbar;
xlabel('r [mm]'); ylabel('z [mm]');
title(sprintf('T(r,z,t(end)) at t = %.3f s', t(end)));

figure('Name','Peak temperature vs time','Color','w');
plot(t, peakT, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]'); ylabel('Peak Temperature [K]');

figure('Name','Radial gradient vs time','Color','w');
plot(t, radial_grad_max, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]'); ylabel('max |dT/dr| [K/m]');

figure('Name','Axial profile on axis (final)','Color','w');
plot(z*1e3, axial_profile_final, 'LineWidth', 1.5);
grid on;
xlabel('z [mm]'); ylabel('T(r=0,z,t_{end}) [K]');

post = struct();
post.peak_temperature_K = peakT;
post.radial_gradient_max_K_per_m = radial_grad_max;
post.axial_profile_final_K = axial_profile_final;
post.thermal_lens_proxy_per_m = lens_proxy;
end
