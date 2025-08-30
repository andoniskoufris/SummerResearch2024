% this function will do rejection sampling to create n points that are distributed approximately according to a Lorentzian distribution
function accepted_points = rejection_sampling(n,A)

accepted_points = zeros(n,1);

j = 0;

while j<n
    x = rand();
    y = rand();
    
    if y < lorentzian_function(x,A)
        j = j+1;
        accepted_points(j) = x;
    else
        continue
    end

end

end