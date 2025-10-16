function r = setRadius(r_rf, t_si, t_tox1, t_tn1, t_tox2, t_ctn, t_box, t_top_channel)
r = struct();
r_0 = r_rf + t_si;
r.r_rf = r_rf;
r.t_si = t_si;
r.r_si = r_rf + t_si;
r.r_0 = r_0;
r.t_tox1 = t_tox1;
r.t_tn1 = t_tn1;
r.t_tox2 = t_tox2;
r.r_tox1 = r_0 + t_tox1;
r.r_tn1 = r.r_tox1 + t_tn1;
r.r_tox2 = r.r_tn1 + t_tox2;
r.r_tox = r.r_tox2;          % top = box, blo
r.t_ctn = t_ctn;
r.r_ctn = r.r_tox + t_ctn;
r.t_box = t_box;
r.r_box = r.r_ctn + t_box;
r.t_top_channel = t_top_channel;
r.r_top_channel = r.r_box + t_top_channel;
R = r.r_top_channel;
n_num = 1e4;
r.n_num = n_num;
r.dr = (R-r_0)/n_num; dr = r.dr;
r.r = ((r_0-1e-7):dr:R)';
r.l_WL = 22e-7;
end