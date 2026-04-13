function laser_coating_removal_dose_finder()
% LASER_COATING_REMOVAL_DOSE_FINDER
%   Dose-finding GUI for selective coating removal
%   (gray/metallic acrylic primer ~50 um on steel).
%
%   The user increments fluence by adjusting beam size (D4sigma)
%   and visually inspects results to find the lowest dose that
%   cleanly strips primer without visible steel damage or heat tinting.
%
%   Fixed laser: 2.5 mJ, 13 ns, 1000 Hz (Gaussian beam).
%   Primary control variable: beam diameter D4sigma [um].
%
%   MATLAB R2013a/b compatible:
%     - All graphics properties via get()/set() — no dot notation
%     - xline/yline replaced with explicit line() calls
%     - mean(...,'omitnan') replaced with manual NaN removal
%     - Colormap 'jet' (no 'turbo' / 'parula')
%     - 'SelectionChangeFcn'  (R2013 uibuttongroup property name)
%     - No FaceAlpha on patch/fill objects
%     - Bar chart via individual bar() calls — no isprop / CData
%     - round(x,N) replaced with round(x*10^N)/10^N
%     - string() replaced with num2str() / sprintf()
%
%   Usage:
%       laser_coating_removal_dose_finder()

% =========================================================================
% SECTION 1 — Constants & Defaults
% =========================================================================

LASER_E_MJ   = 2.5;      % Fixed pulse energy  [mJ]
LASER_TAU_NS = 13;        % Fixed pulse width   [ns]
LASER_RR_HZ  = 1000;      % Fixed rep rate      [Hz]

LASER_E_J    = LASER_E_MJ  * 1e-3;   % [J]
LASER_TAU_S  = LASER_TAU_NS * 1e-9;  % [s]  (informational)

PRIMER_T_UM  = 50;     % Nominal primer thickness   [um]

DEF_FTH      = 0.50;   % Default ablation threshold [J/cm^2]
DEF_FMIN     = 1.0;    % Default sweep start        [J/cm^2]
DEF_FMAX     = 8.0;    % Default sweep end          [J/cm^2]
DEF_FSTEP    = 1.0;    % Default sweep step         [J/cm^2]
DEF_OVERLAP  = 20;     % Default pulse overlap      [%]
DEF_CMAP     = 'jet';  % R2013-compatible colormap  (no turbo/parula)

% Recommended fluence bins  [F_low, F_high]  J/cm^2
REC_BINS = [1 2; 2 3; 3 4; 4 5; 5 6; 6 8];
N_BINS   = size(REC_BINS, 1);

% Spacing-mode identifiers
SM_STRIP = 1;   % spacing referenced to predicted strip / ablation width
SM_BEAM  = 2;   % spacing referenced to physical beam diameter D4sigma

% Scan-mode identifiers
SC_RASTER = 1;
SC_SERP   = 2;

% -------------------------------------------------------------------------
% Mutable GUI state (shared across nested functions via closure)
% -------------------------------------------------------------------------
st_fth      = DEF_FTH;
st_overlap  = DEF_OVERLAP;
st_spacingM = SM_STRIP;
st_scanM    = SC_RASTER;   % reserved for future path-planning integration
st_E_J      = LASER_E_J;   % may be overridden by LaserOverrides
st_RR_Hz    = LASER_RR_HZ; % may be overridden by LaserOverrides

% =========================================================================
% SECTION 2 — Laser Physics (nested helper functions)
% =========================================================================

    function d4s = fluence_to_d4s(F, E)
        % D4sigma [cm] for Gaussian peak fluence F [J/cm^2], energy E [J].
        % F0 = 8*E / (pi * D4s^2)  =>  D4s = sqrt(8*E / (pi*F))
        if F <= 0 || E <= 0
            d4s = NaN; return;
        end
        d4s = sqrt(8 * E / (pi * F));
    end

    function F0 = d4s_to_fluence(d4s, E)
        % Peak fluence [J/cm^2] from D4sigma [cm] and energy [J].
        w0 = d4s / 2;
        F0 = 2 * E / (pi * w0^2);
    end

    function ra = abl_r_cm(d4s_cm, F0, Fth)
        % Ablation radius [cm]: outermost radius where Gaussian >= Fth.
        w0 = d4s_cm / 2;
        if F0 <= Fth
            ra = 0; return;
        end
        ra = w0 * sqrt(log(F0 / Fth) / 2);
    end

    function dx = pulse_dx_um(d4s_um, F0, Fth, ovlp, smode)
        % Pulse-to-pulse spacing [um] for given spacing mode.
        d4s_cm = d4s_um * 1e-4;
        if smode == SM_STRIP
            ra_cm = abl_r_cm(d4s_cm, F0, Fth);
            ref_um = ra_cm * 2 * 1e4;           % strip diameter [um]
            if ref_um < 1
                ref_um = d4s_um * 0.5;          % fallback: half beam diam
            end
        else
            ref_um = d4s_um;                    % physical beam diameter
        end
        dx = ref_um * (1 - ovlp / 100);
        if dx < 0.5, dx = 0.5; end
    end

    function v = scan_speed_mms(dx_um, rr)
        % Scan speed [mm/s] = spacing [mm] * rep rate [Hz].
        v = dx_um * 1e-3 * rr;
    end

    function y = rnd1(x)
        % Round to 1 decimal place (replaces round(x,1) — introduced R2014b).
        y = round(x * 10) / 10;
    end

    function y = rnd2(x)
        % Round to 2 decimal places (replaces round(x,2)).
        y = round(x * 100) / 100;
    end

% =========================================================================
% SECTION 3 — Main Figure
% =========================================================================

hFig = figure( ...
    'Name',            'Laser Coating Removal — Dose Finder', ...
    'NumberTitle',     'off', ...
    'MenuBar',         'none', ...
    'ToolBar',         'figure', ...
    'Color',           [0.92 0.92 0.92], ...
    'Position',        [30 30 1400 840], ...
    'CloseRequestFcn', @cbClose, ...
    'Resize',          'on');

% =========================================================================
% SECTION 4 — Left Configuration Panel
% =========================================================================

pCfg = uipanel('Parent', hFig, ...
    'Title',           'Configuration', ...
    'FontWeight',      'bold', ...
    'Units',           'normalized', ...
    'Position',        [0.000 0.000 0.190 1.000], ...
    'BackgroundColor', [0.92 0.92 0.92]);

BG = [0.92 0.92 0.92];  % panel background colour shorthand
yy = 0.967;             % running y-cursor (normalised, top-down)

% --- Laser (Fixed) -------------------------------------------------------
uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'LASER  (Fixed)', ...
    'FontWeight',      'bold', ...
    'ForegroundColor', [0.00 0.28 0.60], ...
    'Units',           'normalized', ...
    'Position',        [0.02 yy 0.96 0.030], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', BG);
yy = yy - 0.033;

uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          sprintf('Energy :  %.1f mJ', LASER_E_MJ), ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.028], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BG);
yy = yy - 0.031;

uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          sprintf('Pulse width :  %d ns', LASER_TAU_NS), ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.028], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BG);
yy = yy - 0.031;

uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          sprintf('Rep rate :  %d Hz', LASER_RR_HZ), ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.028], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BG);
yy = yy - 0.042;

% --- Laser Overrides -----------------------------------------------------
uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'LASER OVERRIDES', ...
    'FontWeight',      'bold', ...
    'ForegroundColor', [0.60 0.10 0.10], ...
    'Units',           'normalized', ...
    'Position',        [0.02 yy 0.96 0.030], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', BG);
yy = yy - 0.033;

uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'Override Energy (mJ): [blank=fixed]', ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.028], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BG);
yy = yy - 0.037;

edEovr = uicontrol('Parent', pCfg, 'Style', 'edit', ...
    'String',          '', ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.56 0.038], ...
    'TooltipString',   'Leave blank to use fixed 2.5 mJ', ...
    'Callback',        @cbOverrideChange);
yy = yy - 0.040;

uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'Override Rep Rate (Hz): [blank=fixed]', ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.028], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BG);
yy = yy - 0.037;

edRRovr = uicontrol('Parent', pCfg, 'Style', 'edit', ...
    'String',          '', ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.56 0.038], ...
    'TooltipString',   'Leave blank to use fixed 1000 Hz', ...
    'Callback',        @cbOverrideChange);
yy = yy - 0.044;

% --- Material ------------------------------------------------------------
uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'MATERIAL', ...
    'FontWeight',      'bold', ...
    'ForegroundColor', [0.00 0.28 0.60], ...
    'Units',           'normalized', ...
    'Position',        [0.02 yy 0.96 0.030], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', BG);
yy = yy - 0.033;

uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'Ablation threshold F_th (J/cm^2):', ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.028], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BG);
yy = yy - 0.037;

edFth = uicontrol('Parent', pCfg, 'Style', 'edit', ...
    'String',          num2str(DEF_FTH), ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.56 0.038], ...
    'Callback',        @cbFthChange);
yy = yy - 0.040;

uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'Primer thickness (um):', ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.028], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BG);
yy = yy - 0.037;

% Primer thickness is informational — shown as a read-only display field.
uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',              sprintf('%.0f um (informational)', PRIMER_T_UM), ...
    'Units',               'normalized', ...
    'Position',            [0.04 yy 0.92 0.038], ...
    'HorizontalAlignment', 'left', ...
    'ForegroundColor',     [0.40 0.40 0.40], ...
    'BackgroundColor',     BG);
yy = yy - 0.048;

% --- Spacing Mode (uibuttongroup — SelectionChangeFcn for R2013) ----------
bgSpacing = uibuttongroup('Parent', pCfg, ...
    'Title',              'Spacing Reference', ...
    'Units',              'normalized', ...
    'Position',           [0.02 yy-0.132 0.96 0.140], ...
    'SelectionChangeFcn', @cbSpacingMode, ...
    'BackgroundColor',    BG);

rbStrip = uicontrol('Parent', bgSpacing, 'Style', 'radiobutton', ...
    'String',          'Strip / ablation width', ...
    'Tag',             'strip', ...
    'Units',           'normalized', ...
    'Position',        [0.04 0.55 0.92 0.36], ...
    'Value',           1, ...
    'BackgroundColor', BG);  %#ok<NASGU>

uicontrol('Parent', bgSpacing, 'Style', 'radiobutton', ...
    'String',          'Beam diameter D4sigma', ...
    'Tag',             'beam', ...
    'Units',           'normalized', ...
    'Position',        [0.04 0.10 0.92 0.36], ...
    'Value',           0, ...
    'BackgroundColor', BG);

yy = yy - 0.155;

% --- Scan Mode -----------------------------------------------------------
bgScan = uibuttongroup('Parent', pCfg, ...    %#ok<NASGU>
    'Title',              'Scan Mode', ...
    'Units',              'normalized', ...
    'Position',           [0.02 yy-0.110 0.96 0.118], ...
    'SelectionChangeFcn', @cbScanMode, ...
    'BackgroundColor',    BG);

uicontrol('Parent', bgScan, 'Style', 'radiobutton', ...
    'String',          'Raster', ...
    'Tag',             'raster', ...
    'Units',           'normalized', ...
    'Position',        [0.04 0.55 0.92 0.36], ...
    'Value',           1, ...
    'BackgroundColor', BG);

uicontrol('Parent', bgScan, 'Style', 'radiobutton', ...
    'String',          'Serpentine', ...
    'Tag',             'serpentine', ...
    'Units',           'normalized', ...
    'Position',        [0.04 0.10 0.92 0.36], ...
    'Value',           0, ...
    'BackgroundColor', BG);

yy = yy - 0.132;

% --- Dose Sweep Parameters -----------------------------------------------
uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'DOSE SWEEP', ...
    'FontWeight',      'bold', ...
    'ForegroundColor', [0.00 0.28 0.60], ...
    'Units',           'normalized', ...
    'Position',        [0.02 yy 0.96 0.030], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', BG);
yy = yy - 0.033;

uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'F min (J/cm^2):', ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.028], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BG);
yy = yy - 0.036;

edFmin = uicontrol('Parent', pCfg, 'Style', 'edit', ...
    'String',          num2str(DEF_FMIN), ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.56 0.038]);
yy = yy - 0.042;

uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'F max (J/cm^2):', ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.028], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BG);
yy = yy - 0.036;

edFmax = uicontrol('Parent', pCfg, 'Style', 'edit', ...
    'String',          num2str(DEF_FMAX), ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.56 0.038]);
yy = yy - 0.042;

uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'F step (J/cm^2):', ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.028], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BG);
yy = yy - 0.036;

edFstep = uicontrol('Parent', pCfg, 'Style', 'edit', ...
    'String',          num2str(DEF_FSTEP), ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.56 0.038]);
yy = yy - 0.042;

uicontrol('Parent', pCfg, 'Style', 'text', ...
    'String',          'Overlap (%):', ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.028], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BG);
yy = yy - 0.036;

edOverlap = uicontrol('Parent', pCfg, 'Style', 'edit', ...
    'String',          num2str(DEF_OVERLAP), ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.56 0.038]);
yy = yy - 0.050;

uicontrol('Parent', pCfg, 'Style', 'pushbutton', ...
    'String',          'Run Dose Sweep', ...
    'FontWeight',      'bold', ...
    'ForegroundColor', [1 1 1], ...
    'BackgroundColor', [0.18 0.55 0.18], ...
    'Units',           'normalized', ...
    'Position',        [0.04 yy 0.92 0.055], ...
    'Callback',        @cbRunSweep);  %#ok<NASGU>

% =========================================================================
% SECTION 5 — Centre-top: Dose Sweep Results Table
% =========================================================================

pSweep = uipanel('Parent', hFig, ...
    'Title',      'Dose Sweep Results', ...
    'FontWeight', 'bold', ...
    'Units',      'normalized', ...
    'Position',   [0.190 0.500 0.415 0.500]);

swColNames  = {'Fluence (J/cm2)', 'D4sigma (um)', 'w0 (um)', ...
               'Strip W (um)', 'dx (um)', 'Speed (mm/s)', 'F0/Fth', 'Status'};
swColFmts   = {'numeric','numeric','numeric','numeric','numeric','numeric','numeric','char'};
swColEdit   = [false false false false false false false false];
swColWidths = {100 90 75 90 70 90 70 100};

tblSweep = uitable('Parent', pSweep, ...
    'ColumnName',     swColNames, ...
    'ColumnFormat',   swColFmts, ...
    'ColumnEditable', swColEdit, ...
    'ColumnWidth',    swColWidths, ...
    'RowName',        {}, ...
    'Data',           {}, ...
    'Units',          'normalized', ...
    'Position',       [0.01 0.01 0.98 0.98], ...
    'CellSelectionCallback', @cbSweepSel);

% =========================================================================
% SECTION 6 — Right-top: Recommended Beam-Size Chart
% =========================================================================

pRec = uipanel('Parent', hFig, ...
    'Title',      'Recommended Beam Size Chart', ...
    'FontWeight', 'bold', ...
    'Units',      'normalized', ...
    'Position',   [0.605 0.500 0.395 0.500]);

recColNames  = {'Fluence Bin', 'D4sigma (um)', 'w0 (um)', 'Speed (mm/s)', 'Apply'};
recColFmts   = {'char', 'numeric', 'numeric', 'numeric', 'logical'};
recColEdit   = [false false false false true];
recColWidths = {120 90 75 90 50};

recData0 = buildRecData(REC_BINS, st_E_J, st_RR_Hz, st_overlap, st_fth);

tblRec = uitable('Parent', pRec, ...
    'ColumnName',     recColNames, ...
    'ColumnFormat',   recColFmts, ...
    'ColumnEditable', recColEdit, ...
    'ColumnWidth',    recColWidths, ...
    'RowName',        {}, ...
    'Data',           recData0, ...
    'Units',          'normalized', ...
    'Position',       [0.01 0.43 0.98 0.55], ...
    'CellEditCallback', @cbRecEdit);

uicontrol('Parent', pRec, 'Style', 'pushbutton', ...
    'String',          'Apply Checked Dose', ...
    'FontWeight',      'bold', ...
    'ForegroundColor', [1 1 1], ...
    'BackgroundColor', [0.10 0.46 0.78], ...
    'Units',           'normalized', ...
    'Position',        [0.05 0.31 0.90 0.10], ...
    'Callback',        @cbApplyDose);  %#ok<NASGU>

lblStatus = uicontrol('Parent', pRec, 'Style', 'text', ...
    'String',              'No dose applied.', ...
    'Units',               'normalized', ...
    'Position',            [0.04 0.19 0.92 0.10], ...
    'HorizontalAlignment', 'center', ...
    'ForegroundColor',     [0.40 0.40 0.40]);

% Bar chart axes (inside pRec, below status label)
axBar = axes('Parent', pRec, ...
    'Units',    'normalized', ...
    'Position', [0.12 0.02 0.86 0.16]);
set(axBar, 'Box', 'on');
xlabel(axBar, 'Fluence Bin');
ylabel(axBar, 'D4sigma (um)');
title(axBar,  'D4sigma per Fluence Bin');

% =========================================================================
% SECTION 7 — Centre-bottom: Fluence Profile Plot
% =========================================================================

pPlot = uipanel('Parent', hFig, ...
    'Title',      'Fluence Profile & Ablation Preview', ...
    'FontWeight', 'bold', ...
    'Units',      'normalized', ...
    'Position',   [0.190 0.000 0.415 0.500]);

axFlux = axes('Parent', pPlot, ...
    'Units',    'normalized', ...
    'Position', [0.12 0.13 0.84 0.81]);
xlabel(axFlux, 'Radial position (um)');
ylabel(axFlux, 'Fluence (J/cm^2)');
title(axFlux,  'Select a dose to preview profile');
grid(axFlux, 'on');
box(axFlux, 'on');

% =========================================================================
% SECTION 8 — Right-bottom: Metrics Panel
% =========================================================================

pMet = uipanel('Parent', hFig, ...
    'Title',      'Current Dose Metrics', ...
    'FontWeight', 'bold', ...
    'Units',      'normalized', ...
    'Position',   [0.605 0.000 0.395 0.500]);

metRowDefs = { ...
    'D4sigma (beam diam.)',    'um';     ...
    'w0 (beam radius)',        'um';     ...
    'Peak fluence F0',         'J/cm^2'; ...
    'Ablation threshold Fth',  'J/cm^2'; ...
    'Strip width (2*r_abl)',   'um';     ...
    'Pulse spacing dx',        'um';     ...
    'Scan speed',              'mm/s';   ...
    'Fluence ratio F0/Fth',    '—';      ...
    'Risk assessment',         ''};

nMet    = size(metRowDefs, 1);
metHnd  = cell(nMet, 1);
yMet    = 0.940;
dMet    = 0.088;

for km = 1:nMet
    uicontrol('Parent', pMet, 'Style', 'text', ...
        'String',          [metRowDefs{km,1} ':'], ...
        'FontWeight',      'bold', ...
        'Units',           'normalized', ...
        'Position',        [0.02 yMet-0.058 0.54 0.055], ...
        'HorizontalAlignment', 'right');
    metHnd{km} = uicontrol('Parent', pMet, 'Style', 'text', ...
        'String',          '—', ...
        'Units',           'normalized', ...
        'Position',        [0.57 yMet-0.058 0.28 0.055], ...
        'HorizontalAlignment', 'left');
    uicontrol('Parent', pMet, 'Style', 'text', ...
        'String',          metRowDefs{km,2}, ...
        'Units',           'normalized', ...
        'Position',        [0.86 yMet-0.058 0.12 0.055], ...
        'HorizontalAlignment', 'left', ...
        'ForegroundColor', [0.50 0.50 0.50]);
    yMet = yMet - dMet;
end

uicontrol('Parent', pMet, 'Style', 'pushbutton', ...
    'String',   'Export Sweep Results (CSV)', ...
    'Units',    'normalized', ...
    'Position', [0.05 0.01 0.90 0.07], ...
    'Callback', @cbExport);  %#ok<NASGU>

% =========================================================================
% SECTION 9 — Initialise plots
% =========================================================================

updateBarChart();

axes(axFlux);  %#ok<LAXES>
cla;
text(0.5, 0.5, 'Run a sweep or apply a dose to preview', ...
    'Units',               'normalized', ...
    'HorizontalAlignment', 'center', ...
    'Color',               [0.55 0.55 0.55], ...
    'FontSize',            11);

% =========================================================================
% SECTION 10 — Callbacks
% =========================================================================

    % --- Parameter change: ablation threshold ----------------------------
    function cbFthChange(~, ~)
        v = str2double(get(edFth, 'String'));
        if isnan(v) || v <= 0
            set(edFth, 'String', num2str(st_fth));
            return;
        end
        st_fth = v;
        refreshRecTable();
        updateBarChart();
    end

    % --- Laser override values -------------------------------------------
    function cbOverrideChange(~, ~)
        eStr = strtrim(get(edEovr,  'String'));
        rStr = strtrim(get(edRRovr, 'String'));
        if isempty(eStr)
            st_E_J  = LASER_E_J;
        else
            v = str2double(eStr);
            if ~isnan(v) && v > 0
                st_E_J = v * 1e-3;  % mJ to J
            end
        end
        if isempty(rStr)
            st_RR_Hz = LASER_RR_HZ;
        else
            v = str2double(rStr);
            if ~isnan(v) && v > 0
                st_RR_Hz = v;
            end
        end
        refreshRecTable();
        updateBarChart();
    end

    % --- Spacing mode (uibuttongroup SelectionChangeFcn — R2013 name) ----
    function cbSpacingMode(~, event)
        % In R2013 the event arg is an EventData object — use get(), not dot.
        newRb = get(event, 'NewValue');
        tag   = get(newRb, 'Tag');
        if strcmp(tag, 'strip')
            st_spacingM = SM_STRIP;
        else
            st_spacingM = SM_BEAM;
        end
        refreshRecTable();
    end

    % --- Scan mode -------------------------------------------------------
    function cbScanMode(~, event)
        newRb = get(event, 'NewValue');
        tag   = get(newRb, 'Tag');
        if strcmp(tag, 'raster')
            st_scanM = SC_RASTER;
        else
            st_scanM = SC_SERP;
        end
    end

    % --- Run dose sweep --------------------------------------------------
    function cbRunSweep(~, ~)
        fmin_v  = str2double(get(edFmin,    'String'));
        fmax_v  = str2double(get(edFmax,    'String'));
        fstep_v = str2double(get(edFstep,   'String'));
        ovlp_v  = str2double(get(edOverlap, 'String'));
        fth_v   = str2double(get(edFth,     'String'));

        if any(isnan([fmin_v fmax_v fstep_v ovlp_v fth_v]))
            errordlg('One or more sweep parameters are invalid.', 'Sweep Error');
            return;
        end
        if fstep_v <= 0 || fmin_v >= fmax_v || fmin_v <= 0
            errordlg('Require: F_min > 0, F_min < F_max, F_step > 0.', 'Sweep Error');
            return;
        end

        st_overlap = ovlp_v;
        st_fth     = fth_v;

        fluences = fmin_v : fstep_v : fmax_v;
        nF = numel(fluences);
        tbl = cell(nF, 8);

        for fi = 1:nF
            Ft = fluences(fi);
            d4c  = fluence_to_d4s(Ft, st_E_J);   % cm
            d4u  = d4c * 1e4;                     % um
            w0u  = d4u / 2;
            F0i  = d4s_to_fluence(d4c, st_E_J);
            rac  = abl_r_cm(d4c, F0i, fth_v);
            swu  = rac * 2 * 1e4;                 % strip width um
            dxi  = pulse_dx_um(d4u, F0i, fth_v, ovlp_v, st_spacingM);
            vi   = scan_speed_mms(dxi, st_RR_Hz);
            ri   = F0i / fth_v;

            if ri < 1.0
                st_i = 'Below threshold';
            elseif ri < 1.5
                st_i = 'Near threshold';
            elseif ri > 15
                st_i = 'Risk: steel damage';
            elseif ri > 8
                st_i = 'High — check damage';
            else
                st_i = 'OK';
            end

            tbl{fi,1} = rnd2(Ft);
            tbl{fi,2} = rnd1(d4u);
            tbl{fi,3} = rnd1(w0u);
            tbl{fi,4} = rnd1(swu);
            tbl{fi,5} = rnd1(dxi);
            tbl{fi,6} = rnd1(vi);
            tbl{fi,7} = rnd2(ri);
            tbl{fi,8} = st_i;
        end

        set(tblSweep, 'Data', tbl);
        refreshRecTable();
        updateBarChart();
    end

    % --- Sweep table row selected ----------------------------------------
    function cbSweepSel(~, event)
        % event is a struct in R2013 — .Indices is a struct field, not a
        % graphics handle property, so dot access is correct here.
        idx = event.Indices;
        if isempty(idx), return; end
        row = idx(1, 1);
        data = get(tblSweep, 'Data');
        if row < 1 || row > size(data, 1), return; end

        d4v  = data{row, 2};
        dxv  = data{row, 5};
        spv  = data{row, 6};

        if isnumeric(d4v) && ~isnan(d4v) && d4v > 0
            drawFluencePlot(d4v);
            refreshMetrics(d4v, dxv, spv);
        end
    end

    % --- Recommended table checkbox toggled ------------------------------
    function cbRecEdit(~, event)
        idx = event.Indices;
        if isempty(idx) || idx(1,2) ~= 5, return; end
        % Ensure only one row checked at a time
        data = get(tblRec, 'Data');
        for ri = 1:size(data, 1)
            data{ri, 5} = false;
        end
        data{idx(1,1), 5} = true;
        set(tblRec, 'Data', data);
    end

    % --- Apply checked dose ----------------------------------------------
    function cbApplyDose(~, ~)
        data    = get(tblRec, 'Data');
        chosen  = 0;
        for ri = 1:size(data, 1)
            chk = data{ri, 5};
            if isequal(chk, true) || isequal(chk, 1)
                chosen = ri;
                break;
            end
        end
        if chosen == 0
            set(lblStatus, ...
                'String',          'Check one bin row then click Apply.', ...
                'ForegroundColor', [0.70 0.20 0.10]);
            return;
        end

        d4v    = data{chosen, 2};
        spv    = data{chosen, 4};
        binStr = data{chosen, 1};

        d4c  = d4v * 1e-4;
        F0a  = d4s_to_fluence(d4c, st_E_J);
        dxa  = pulse_dx_um(d4v, F0a, st_fth, st_overlap, st_spacingM);

        drawFluencePlot(d4v);
        refreshMetrics(d4v, dxa, spv);

        set(lblStatus, ...
            'String',          sprintf('Applied: %s | D4s=%.0f um | v=%.1f mm/s', ...
                                       binStr, d4v, spv), ...
            'ForegroundColor', [0.00 0.48 0.00]);
    end

    % --- Export sweep to CSV ---------------------------------------------
    function cbExport(~, ~)
        data = get(tblSweep, 'Data');
        if isempty(data)
            msgbox('Run a dose sweep first.', 'Export');
            return;
        end
        [fname, fpath] = uiputfile('*.csv', 'Save Sweep Results as CSV');
        if isequal(fname, 0), return; end
        fid = fopen(fullfile(fpath, fname), 'w');
        if fid < 0
            errordlg('Could not open file for writing.', 'Export Error');
            return;
        end
        fprintf(fid, ...
            'Fluence(J/cm2),D4sigma(um),w0(um),StripW(um),dx(um),Speed(mm/s),F0/Fth,Status\n');
        for ri = 1:size(data, 1)
            r = data(ri, :);
            fprintf(fid, '%.2f,%.1f,%.1f,%.1f,%.1f,%.1f,%.2f,%s\n', ...
                r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, r{7}, r{8});
        end
        fclose(fid);
        msgbox(['Saved: ' fullfile(fpath, fname)], 'Export OK');
    end

    % --- Close figure ----------------------------------------------------
    function cbClose(src, ~)
        delete(src);
    end

% =========================================================================
% SECTION 11 — Update / Refresh Functions
% =========================================================================

    % --- Refresh recommended table with current state --------------------
    function refreshRecTable()
        ov = str2double(get(edOverlap, 'String'));
        if isnan(ov) || ov < 0 || ov >= 100
            ov = st_overlap;
        end
        newData = buildRecData(REC_BINS, st_E_J, st_RR_Hz, ov, st_fth);
        set(tblRec, 'Data', newData);
    end

    % --- Draw Gaussian fluence profile -----------------------------------
    function drawFluencePlot(d4s_um)
        % Use axes(h) + working on gca — safe pattern for R2013.
        axes(axFlux);  %#ok<LAXES>
        cla;

        if isnan(d4s_um) || d4s_um <= 0
            text(0.5, 0.5, 'Invalid D4sigma value', ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'center', ...
                'Color', [0.6 0.6 0.6]);
            return;
        end

        fth_p   = st_fth;
        E_p     = st_E_J;
        d4c_p   = d4s_um * 1e-4;
        w0c_p   = d4c_p / 2;
        w0u_p   = d4s_um / 2;
        F0_p    = d4s_to_fluence(d4c_p, E_p);

        % Radial coordinate
        r_max   = d4s_um * 1.9;
        r_vec   = linspace(-r_max, r_max, 600);

        % Gaussian profile  F(r) = F0 * exp(-2*r^2/w0^2)
        F_vec   = F0_p .* exp(-2 .* (r_vec ./ w0u_p).^2);

        % Plot profile
        plot(r_vec, F_vec, 'b-', 'LineWidth', 2, 'DisplayName', 'Fluence F(r)');
        hold on;

        % Set axis limits before drawing reference lines
        ymax_p = F0_p * 1.18;
        set(gca, 'XLim', [-r_max r_max], 'YLim', [0 ymax_p]);

        % --- Ablation threshold (horizontal line — replaces yline) -------
        line([-r_max r_max], [fth_p fth_p], ...
            'LineStyle', '--', ...
            'Color',     [0.85 0.10 0.10], ...
            'LineWidth',  1.5, ...
            'DisplayName', sprintf('F_{th} = %.2f J/cm^2', fth_p));

        % --- 1/e^2 radius w0 markers (vertical — replaces xline) ---------
        line([-w0u_p -w0u_p], [0 ymax_p], ...
            'LineStyle',        ':', ...
            'Color',            [0.00 0.50 0.00], ...
            'LineWidth',         1.2, ...
            'HandleVisibility', 'off');
        line([w0u_p w0u_p], [0 ymax_p], ...
            'LineStyle',   ':', ...
            'Color',       [0.00 0.50 0.00], ...
            'LineWidth',    1.2, ...
            'DisplayName', sprintf('w_0 = %.0f um', w0u_p));

        % --- Ablation boundary (vertical — replaces xline) ---------------
        ra_p = abl_r_cm(d4c_p, F0_p, fth_p);
        if ra_p > 0
            rau_p = ra_p * 1e4;
            line([-rau_p -rau_p], [0 ymax_p], ...
                'LineStyle',        '-.', ...
                'Color',            [0.85 0.45 0.00], ...
                'LineWidth',         1.3, ...
                'HandleVisibility', 'off');
            line([rau_p rau_p], [0 ymax_p], ...
                'LineStyle',   '-.', ...
                'Color',       [0.85 0.45 0.00], ...
                'LineWidth',    1.3, ...
                'DisplayName', sprintf('Strip edge: +-%.0f um', rau_p));
        end

        hold off;
        xlabel('Radial position (um)');
        ylabel('Fluence (J/cm^2)');
        title(sprintf('Gaussian Profile  |  D4sigma = %.0f um  |  F_0 = %.3f J/cm^2', ...
                       d4s_um, F0_p));
        legend('show', 'Location', 'NorthEast');
        grid on;
        box on;
    end

    % --- Refresh metrics panel -------------------------------------------
    function refreshMetrics(d4s_um, dx_um, spd)
        if isnan(d4s_um) || d4s_um <= 0
            for km2 = 1:nMet
                set(metHnd{km2}, 'String', '—', 'ForegroundColor', [0 0 0]);
            end
            return;
        end

        fth_m    = st_fth;
        d4c_m    = d4s_um * 1e-4;
        w0u_m    = d4s_um / 2;
        F0_m     = d4s_to_fluence(d4c_m, st_E_J);
        ra_m     = abl_r_cm(d4c_m, F0_m, fth_m);
        sw_m     = ra_m * 2 * 1e4;
        ratio_m  = F0_m / fth_m;

        vals = { ...
            sprintf('%.1f', d4s_um), ...
            sprintf('%.1f', w0u_m),  ...
            sprintf('%.4f', F0_m),   ...
            sprintf('%.4f', fth_m),  ...
            sprintf('%.1f', sw_m),   ...
            sprintf('%.1f', dx_um),  ...
            sprintf('%.1f', spd),    ...
            sprintf('%.2f', ratio_m), ...
            ''};

        if ratio_m < 1.0
            riskStr   = 'BELOW THRESHOLD — no ablation';
            riskColor = [0.80 0.00 0.00];
        elseif ratio_m < 1.5
            riskStr   = 'Near threshold — marginal removal';
            riskColor = [0.85 0.50 0.00];
        elseif ratio_m > 15
            riskStr   = 'VERY HIGH — risk of steel damage';
            riskColor = [0.90 0.00 0.00];
        elseif ratio_m > 8
            riskStr   = 'High — check for heat tinting';
            riskColor = [0.85 0.38 0.00];
        else
            riskStr   = 'OK — safe operating range';
            riskColor = [0.00 0.50 0.00];
        end
        vals{9} = riskStr;

        for km3 = 1:(nMet - 1)
            set(metHnd{km3}, 'String', vals{km3}, 'ForegroundColor', [0 0 0]);
        end
        set(metHnd{nMet}, 'String', riskStr, 'ForegroundColor', riskColor);
    end

    % --- Redraw D4sigma bar chart ----------------------------------------
    function updateBarChart()
        % Use axes(h) pattern — bar() does not reliably accept axes handle
        % as first argument in MATLAB R2013 (HG1).
        axes(axBar);  %#ok<LAXES>
        cla;

        cmap  = colormap(DEF_CMAP);
        nC    = size(cmap, 1);
        step  = max(1, floor(nC / N_BINS));
        idxC  = 1 : step : nC;
        if numel(idxC) < N_BINS
            idxC = round(linspace(1, nC, N_BINS));
        end
        idxC = idxC(1:N_BINS);

        lbls = cell(N_BINS, 1);
        hold on;
        for bi = 1:N_BINS
            Fmid_b  = (REC_BINS(bi,1) + REC_BINS(bi,2)) / 2;
            d4c_b   = fluence_to_d4s(Fmid_b, st_E_J);
            d4u_b   = d4c_b * 1e4;
            lbls{bi} = sprintf('%d-%d', REC_BINS(bi,1), REC_BINS(bi,2));
            bar(bi, d4u_b, 0.65, 'FaceColor', cmap(idxC(bi),:), 'EdgeColor', 'none');
        end
        hold off;

        set(gca, 'XTick', 1:N_BINS, 'XTickLabel', lbls);
        ylabel('D4sigma (um)');
        title('D4sigma vs Fluence Bin');
        grid on;
    end

% =========================================================================
% SECTION 12 — Helper: Build Recommended Table Data
% =========================================================================

    function data_out = buildRecData(bins, Ej, rr, ovlp, fth)
        nB       = size(bins, 1);
        data_out = cell(nB, 5);
        for bi2 = 1:nB
            Fmid2  = (bins(bi2,1) + bins(bi2,2)) / 2;
            bstr   = sprintf('%d-%d J/cm^2', bins(bi2,1), bins(bi2,2));
            d4c2   = fluence_to_d4s(Fmid2, Ej);
            d4u2   = d4c2 * 1e4;
            w0u2   = d4u2 / 2;
            F0_2   = d4s_to_fluence(d4c2, Ej);
            dx2    = pulse_dx_um(d4u2, F0_2, fth, ovlp, st_spacingM);
            spd2   = scan_speed_mms(dx2, rr);
            data_out{bi2, 1} = bstr;
            data_out{bi2, 2} = rnd1(d4u2);
            data_out{bi2, 3} = rnd1(w0u2);
            data_out{bi2, 4} = rnd1(spd2);
            data_out{bi2, 5} = false;
        end
    end

end  % laser_coating_removal_dose_finder
