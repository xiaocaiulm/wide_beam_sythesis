function x_projected = projection_con(x)
    theta = angle(x);
    x_projected = exp(1i*theta );
end
