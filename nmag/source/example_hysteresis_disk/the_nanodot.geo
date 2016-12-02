algebraic3d

solid fincyl = cylinder ( 0, 0, 1; 0, 0, -1; 1.0 )
	and plane (0, 0, -0.1; 0, 0, -1)
	and plane (0, 0, 0.1; 0, 0, 1) -maxh = 0.05;

tlo fincyl;
