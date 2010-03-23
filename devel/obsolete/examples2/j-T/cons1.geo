algebraic3d

solid left = orthobrick (-200,-50,-5;-50,50,5);
solid right = orthobrick (50,-50,-5;200,50,5);

#triangle left (in constrction): tl
solid tl_base = orthobrick (-50,-50,-5;0,50,5);
#cut plane lower
solid tl_cutlow = plane (-50,-50,0; 0.8,-1,0);
solid tl_cuthigh = plane (-50,50,0; 0.8,1,0);

solid tl = tl_base and tl_cutlow and tl_cuthigh;



#triangle right (in constrction): tr
solid tr_base = orthobrick (0,-50,-5;50,50,5);
#cut plane lower
solid tr_cutlow = plane (50,-50,0; -0.8,-1,0);
solid tr_cuthigh = plane (50,50,0; -0.8,1,0);

solid tr = tr_base and tr_cutlow and tr_cuthigh;

solid all = left or tl or tr or right -maxh=10;
	
tlo all; 
#tlo tr;
#tlo tl;
#tlo left;
#tlo right;

 