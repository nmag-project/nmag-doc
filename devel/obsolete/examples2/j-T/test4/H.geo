algebraic3d

solid left = orthobrick (-500,-25,-10;-25,25,10) -maxh=15;
solid middle = orthobrick (-25,-10,-10;25,10,10) -maxh=15;
solid right = orthobrick (25,-25,-10;500,25,10) -maxh=15;

solid substrate = orthobrick (-700,-50,-110;700,50,-10) -maxh=15;

solid conductor = left or middle or right;

solid contact_left = orthobrick (-500,-25,10;-400,25,40) -maxh=15;


solid contact_left_wedge = orthobrick (-530,-25,-10;-500,25,40) and plane (-505,0,40; -1,0,1);

solid contact_right = orthobrick (400,-25,10;1000,25,40) -maxh=15;

solid contact_left_lower = orthobrick (-1000,-25,-10;-520,25,20) -maxh=15;

solid contact_right = orthobrick (400,-25,10;500,25,40) -maxh=15;

solid contact_right_wedge = orthobrick (500,-25,-10;530,25,40) and 
plane (505,0,40; 1,0,1);

solid contact_right_lower = orthobrick (520,-25,-10;1000,25,20) -maxh=15;

solid contacts_left = contact_left or  contact_left_wedge or contact_left_lower;

solid contacts_right = contact_right or contact_right_wedge or contact_right_lower;

solid contacts = contacts_left or contacts_right;
	
tlo contacts;
tlo conductor; 
tlo substrate;

 