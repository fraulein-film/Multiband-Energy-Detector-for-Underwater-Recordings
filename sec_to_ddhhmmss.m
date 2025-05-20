function time_str = sec_to_ddhhmmss(seconds)
%
% Converts time in seconds to 'dd:hh:mm:ss' format.
%
% --------------------------- Input Parameters ----------------------------
%
% seconds: Time in seconds
%
% --------------------------- Output Parameters ---------------------------
%
% time_str: Time in format 'dd:hh:mm:ss'
% d - days
% h - hours
% m - minutes
% s - seconds
%
% -------------------------------------------------------------------------

% Calculate d, h, m, s
days = floor(seconds / 86400);
hours = mod( floor(seconds / 3600), 24);
minutes = mod( floor(seconds / 60), 60);
remaining_seconds = mod( seconds, 60);

% Create the time string
time_str = sprintf('%02d:%02d:%02d:%02d', days, hours, minutes, round(remaining_seconds));

end