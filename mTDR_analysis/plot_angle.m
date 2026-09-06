clear; close all; clc;

%% =========================
% 1. Define pairwise angles
% ==========================
% ang_SP: angle between saccade and perceptual axes
% ang_SS: angle between saccade and stimulus axes
% ang_PS: angle between perceptual and stimulus axes
savepath = 'F:\OneDrive\paperfig\mTDR\';
% savename = [savepath,'D_IPS_angle.pdf'];
ang_SP = 43;
ang_SS = 71;
ang_PS = 75;

% savename = [savepath,'F_IPS_angle.pdf'];
% ang_SP = 74;
% ang_SS = 74;
% ang_PS = 60;

% savename = [savepath,'D_STS_angle.pdf'];
% ang_SP = 85;
% ang_SS = 83;
% ang_PS = 53;

% savename = [savepath,'F_STS_angle.pdf'];
% ang_SP = 91;
% ang_SS = 89;
% ang_PS = 58;
%% ==========================================
% 2. Construct three 3D unit vectors
%    perceptual axis is the reference axis
% ==========================================

% Perceptual axis: fixed as x-axis in 3D
v_perc = [1 0 0];

% Saccade axis: placed in x-z plane
v_sacc = [cosd(ang_SP), 0, sind(ang_SP)];
v_sacc = v_sacc / norm(v_sacc);

% Stimulus axis: solve so that
% angle(v_stim, v_perc) = ang_PS
% angle(v_stim, v_sacc) = ang_SS

x = cosd(ang_PS);
z = (cosd(ang_SS) - x*v_sacc(1)) / v_sacc(3);

% Choose negative y so the stimulus axis points "downward"
% If you want the mirrored version, change the sign to +
y = -sqrt(max(0, 1 - x^2 - z^2));

v_stim = [x, y, z];
v_stim = v_stim / norm(v_stim);

%% =========================
% 3. Figure and style setup
% ==========================
fig = figure('Color','w','Position',[100 100 950 700]);
ax = axes(fig);
hold(ax,'on');
axis(ax,'equal');
axis(ax,'off');

% Display lengths (can be adjusted independently)
L_perc = 1.2;
L_sacc = 1.0;
L_stim = 0.8;   % slightly longer so stimulus axis does not look too short

% Plane appearance
planeAlpha1 = 0.14;
planeAlpha2 = 0.12;
planeColor1 = [0.82 0.82 0.82];
planeColor2 = [0.88 0.88 0.88];
edgeColor   = [0.72 0.72 0.72];

% Line/text appearance
lineColor   = [0 0 0];
fontSizeLab = 13;
fontSizeAng = 12;

%% ==========================================
% 4. Draw task-related planes
% ==========================================

% Plane 1: saccade-perceptual plane
% P1 = [0 0 0;
%       1.02*L_sacc*v_sacc;
%       1.02*(L_sacc*v_sacc + 0.72*L_perc*v_perc);
%       1.02*(0.72*L_perc*v_perc)];
% 
% patch('Vertices',P1, ...
%       'Faces',[1 2 3 4], ...
%       'FaceColor',planeColor1, ...
%       'FaceAlpha',planeAlpha1, ...
%       'EdgeColor',edgeColor, ...
%       'LineWidth',0.8);

% Plane 2: perceptual-stimulus plane
planeScale_perc = 0.90;
planeScale_stim = 0.85;

P2 = [0 0 0;
      planeScale_perc * L_perc * v_perc;
      planeScale_perc * L_perc * v_perc + planeScale_stim * L_stim * v_stim;
      planeScale_stim * L_stim * v_stim];

patch('Vertices',P2, ...
      'Faces',[1 2 3 4], ...
      'FaceColor',planeColor2, ...
      'FaceAlpha',planeAlpha2, ...
      'EdgeColor',edgeColor, ...
      'LineWidth',0.8);

%% ====================
% 5. Draw axes
% ====================
drawAxis(v_sacc, L_sacc, lineColor);
drawAxis(v_perc, L_perc, lineColor);
drawAxis(v_stim, L_stim, lineColor);

scatter3(0,0,0,20,'k','filled');

%% ====================
% 6. Add axis labels
% ====================
text(1.04*L_sacc*v_sacc(1)+0.03, ...
     1.04*L_sacc*v_sacc(2)+0.02, ...
     1.04*L_sacc*v_sacc(3)+0.02, ...
     'saccade choice axis', ...
     'FontSize',fontSizeLab, ...
     'Interpreter','none', ...
     'Rotation',62);

text(0.70*L_perc*v_perc(1)+0.02, ...
     0.70*L_perc*v_perc(2)+0.02, ...
     0.70*L_perc*v_perc(3)+0.01, ...
     'perceptual choice axis', ...
     'FontSize',fontSizeLab, ...
     'Interpreter','none');

text(0.83*L_stim*v_stim(1)+0.04, ...
     0.83*L_stim*v_stim(2)-0.03, ...
     0.83*L_stim*v_stim(3)-0.02, ...
     'stimulus axis', ...
     'FontSize',fontSizeLab, ...
     'Interpreter','none', ...
     'Rotation',-18);

%% ====================
% 7. Draw angle arcs
% ====================
drawArc3D(v_perc, v_stim, 0.24, sprintf('%d', ang_PS), [0.03 -0.01 -0.01], fontSizeAng);
drawArc3D(v_sacc, v_stim, 0.33, sprintf('%d', ang_SS), [-0.03 -0.02 0.00], fontSizeAng);
drawArc3D(v_sacc, v_perc, 0.42, sprintf('%d', ang_SP), [0.01 0.00 0.02], fontSizeAng);

%% ====================
% 8. Camera / view
% ====================
% You can change elevation here to make the figure more open.
% The helper function forceVectorHorizontal ensures that the
% perceptual axis remains horizontal on screen.
view(ax,[0 25]);
camproj(ax,'orthographic');
forceVectorHorizontal(ax, v_perc);
% 
% % % Make perceptual axis horizontal in the final rendered screen view
% % forceVectorHorizontal(ax, v_perc);
% % normal of the saccade-perceptual plane
% n_SP = cross(v_sacc, v_perc);
% n_SP = n_SP / norm(n_SP);
% 
% % in-plane tilt direction
% tiltDir = cross(n_SP, v_perc);
% tiltDir = tiltDir / norm(tiltDir);
% 
% % adjust alpha to balance angle fidelity and 3D appearance
% alpha = 0.18;
% 
% viewDir = n_SP + alpha * tiltDir;
% viewDir = viewDir / norm(viewDir);
% 
% camtarget(ax,[0 0 0]);
% campos(ax, -4*viewDir);     % if mirrored, change to +4*viewDir
% camproj(ax,'orthographic');
% 
% % keep perceptual axis horizontal on screen
% forceVectorHorizontal(ax, v_perc);
%% ====================
% 9. Final limits
% ====================
xlim([-0.20 1.65]);
ylim([-1.05 0.95]);
zlim([-0.15 1.15]);

lighting gouraud;
camlight headlight;

%% ====================
% 10. Export (optional)
% ====================
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [7 5.2]);          % PDF page size
set(gcf, 'PaperPosition', [0 0 7 5.2]);  % figure occupies the full page
set(gcf, 'PaperPositionMode', 'manual');

print(gcf, savename, '-dpdf', '-painters');

%% ====================
% Local functions
% ====================

function drawAxis(v, L, c)
    quiver3(0,0,0, ...
        L*v(1), L*v(2), L*v(3), 0, ...
        'Color',c, ...
        'LineWidth',1.5, ...
        'MaxHeadSize',0.12);
end

function drawArc3D(v1, v2, r, txt, txtOffset, fs)
    v1 = v1 / norm(v1);
    v2 = v2 / norm(v2);

    theta = acos(max(min(dot(v1,v2),1),-1));

    e1 = v1;
    temp = v2 - dot(v2,e1)*e1;
    e2 = temp / norm(temp);

    t = linspace(0, theta, 120)';
    pts = r*(cos(t)*e1 + sin(t)*e2);

    plot3(pts(:,1), pts(:,2), pts(:,3), ...
        'k', 'LineWidth',1.4);

    midIdx = round(numel(t)*0.58);
    pm = pts(midIdx,:) + txtOffset;

    text(pm(1), pm(2), pm(3), txt, ...
        'FontSize',fs, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
end

function forceVectorHorizontal(ax, v_horizontal)
    % Force the projection of v_horizontal to be horizontal on screen.

    v_horizontal = v_horizontal(:)' / norm(v_horizontal);

    camPos = campos(ax);
    camTar = camtarget(ax);

    % Viewing direction from camera to target
    viewDir = camTar - camPos;
    viewDir = viewDir / norm(viewDir);

    % Set camera-up vector so that the projected v_horizontal becomes horizontal
    newUp = cross(v_horizontal, viewDir);

    if norm(newUp) < 1e-8
        warning('Selected vector is nearly parallel to the viewing direction.');
        return;
    end

    newUp = newUp / norm(newUp);
    camup(ax, newUp);
end