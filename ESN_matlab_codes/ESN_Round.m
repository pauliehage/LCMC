function out_ = ESN_Round(in_, accuracy_, type_)

if nargin == 1
    accuracy_ = 1;
    type_ = 'round';
end

if nargin == 2
    type_ = 'round';
end

in_ = double(in_);

if (strcmp(type_, 'round'))
    out_ = round(in_ ./ accuracy_) .* accuracy_;
elseif (strcmp(type_, 'ceil'))
    out_ = ceil(in_ ./ accuracy_) .* accuracy_;
elseif (strcmp(type_, 'floor'))
    out_ = floor(in_ ./ accuracy_) .* accuracy_;
else
    out_ = round(in_ ./ accuracy_) .* accuracy_;
end
