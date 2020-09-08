function []=global_init();
global sin_table;
global cos_table;
sin_table = [0,6,11,15,16,15,11, 6, 0,-6,-11,-15,-16,-15,-11,-6];
cos_table = [16,15,11, 6, 0,-6,-11,-15,-16,-15,-11,-6, 0, 6,11,15];
end