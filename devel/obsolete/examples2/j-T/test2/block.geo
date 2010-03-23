algebraic3d

solid Py = orthobrick (-500,-25,-10;500,25,10) -maxh=20;

solid substrate = orthobrick (-500,-25,-110;500,25,-10) -maxh=20;

solid all = Py or substrate;
	
#tlo all; 

tlo substrate;
tlo Py;
#tlo tr;
#tlo tl;
#tlo left;
#tlo right;

 