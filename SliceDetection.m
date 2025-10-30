function i = SliceDetection(r, location, N)

    i = ceil((location-r.r_o2)/r.t_ctn*N);

end