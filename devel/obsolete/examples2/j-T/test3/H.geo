algebraic3d

solid left = orthobrick (-500,-25,-10;-25,25,10) -maxh=5;
solid middle = orthobrick (-25,-10,-10;25,10,10) -maxh=5;
solid right = orthobrick (25,-25,-10;500,25,10) -maxh=5;



solid all = left or middle or right;
	
tlo all; 
#tlo tr;
#tlo tl;
#tlo left;
#tlo right;

 