function [wpts] = Radial_Loiter_Plan(theta,r_array,alt_array,direction)
    if direction == -1
        r_array = fliplr(r_array);
    end
    [X,Y] = pol2cart(theta,r_array);
    Z = alt_array;
    wpts = [X;Y;Z];
end

