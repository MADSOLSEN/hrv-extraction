function x_out = change_EKG_direction(x)
% this functions uses the gradient to determine the direction of the EKG
% signal. The direction is changed if the peak is pointed down. 

x_above_zero = x((x > 0) & abs(gradient(x) > 0));
x_below_zero = x((x < 0) & abs(gradient(x) > 0));

x_diff_above = median(abs(diff(x_above_zero)));
x_diff_below = median(abs(diff(x_below_zero)));

if x_diff_below > x_diff_above
    x = -x;
end

x_out = x; 

end