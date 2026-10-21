function result = gfz_http_get(url, timeoutSeconds)
%GFZ_HTTP_GET Only HTTP work; safe to call on a PROCESS worker, without app.
% Do not call this function on backgroundPool/ThreadPool on older MATLAB.
opts = matlab.net.http.HTTPOptions('ConnectTimeout', 15, ...
    'ResponseTimeout', timeoutSeconds, 'ConvertResponse', false);
% DataTimeout is not available in every MATLAB release.
if isprop(opts, 'DataTimeout')
    opts.DataTimeout = timeoutSeconds;
end
request = matlab.net.http.RequestMessage('GET');
response = request.send(matlab.net.URI(url), opts);
bytes = response.Body.Data;
% Some releases still return textual Content-Type bodies as string/char
% with ConvertResponse=false; binary NetCDF remains uint8.
if isstring(bytes) || ischar(bytes)
    bytes = unicode2native(char(bytes), 'UTF-8');
else
    bytes = uint8(bytes);
end
result = struct('Status', double(response.StatusCode), ...
    'Bytes', bytes, 'RetryAfter', 0);
fields = response.getFields('Retry-After');
if ~isempty(fields)
    value = char(fields(1).Value);
    delay = str2double(value);
    if isnan(delay)
        try
            t = datetime(value, 'InputFormat', 'eee, dd MMM yyyy HH:mm:ss ''GMT''', ...
                'Locale', 'en_US', 'TimeZone', 'UTC');
            delay = seconds(t - datetime('now', 'TimeZone', 'UTC'));
        catch
            delay = 0;
        end
    end
    result.RetryAfter = max(0, delay);
end
end
