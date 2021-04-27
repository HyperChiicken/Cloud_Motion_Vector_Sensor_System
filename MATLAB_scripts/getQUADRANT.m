function QUADRANT = getQUADRANT(theta)
%% This function receives an input angle in radiance and finds the quadrant it belongs to
QUADRANT = cell(numel(theta),1);

for idx = 1:numel(theta)
    x_val = cos(theta(idx));
    y_val = sin(theta(idx));
    
    if y_val > 0 && x_val > 0
        QUADRANT{idx,1} = 'Q1';
    elseif y_val > 0 && x_val < 0
        QUADRANT{idx,1} = 'Q2';
    elseif y_val < 0 && x_val < 0
        QUADRANT{idx,1} = 'Q3';
    elseif y_val < 0 && x_val > 0
        QUADRANT{idx,1} = 'Q4';
    end
end