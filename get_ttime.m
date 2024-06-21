function [t]=get_ttime(dist,depth,phase)

%dist=4 ; depth = 10; phase = 'P'

% usage 
% this script get the travel time (sec) from the taup API given the ::
% dist = source distance in degree
% depth = source depth in km
% phase = phase name ( e.g. Pn, Pg, Sn)


dist=num2str(dist); depth=num2str(depth);
%phase=num2str(phase);
mod='ak135'; % velocity model used

try
    full_url = ['https://service.iris.edu/irisws/traveltime/1/query?model=' ...
    mod '&phases=' phase '&evdepth=' depth '&distdeg=' dist ...
    '&noheader=true&traveltimeonly=true'];
    tout = weboptions('Timeout',120); % set the timout value
    str = webread(full_url,tout);
    %t=str2double(str(293:299));
    %t=str2double(str);
    t1 = strsplit(str, ' ');
    t1 = str2double(t1);
    
%     if dist > 2.4 & strcmp(phase, 'S')
%         t = max(t1);
%     else
%         t = min(t1);
%     end
    t = min(t1);
catch ME
    % Handle the error
    disp(['Error occurred: ' ME.message]);
    t = NaN;
    
end


end

