function set_fonts(axis)

if nargin < 1
    axis = gca;
end

set(get(axis, 'XLabel'), 'FontSize', 6);%6);
set(get(axis, 'YLabel'), 'FontSize', 6);%6);
set(axis, 'FontSize', 6);%6);
set(axis, 'FontName', 'Arial');

set(get(axis, 'Title'), 'FontSize', 8);%8);


end