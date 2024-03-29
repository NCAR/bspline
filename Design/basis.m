
global dx = 1;

function y = basis (x)
	y = 0;
	global dx;
	x = x / dx;
	if (abs(x) < 2)
		if (x < 0)
			y = 0.25 * (2 + x)^3;
		else
			y = 0.25 * (2 - x)^3;
		endif
		if (abs(x) < 1)
			if (x < 0)
				y = y - (1 + x)^3;
			else
				y = y - (1 - x)^3;
			endif
		endif
	endif
endfunction

function y = dbasis (x)
	y = 0;
	global dx;
	x = x / dx;
	if (abs(x) < 2)
		if (x < 0)
			y = 3* 0.25 * (2 + x)^2;
		else
			y = -3 * 0.25 * (2 - x)^2;
		endif
		if (abs(x) < 1)
			if (x < 0)
				y = y - 3 * (1 + x)^2;
			else
				y = y + 3 * (1 - x)^2;
			endif
		endif
	endif

#	# If x == 0, sign(x) is zero, which is ok
#	if (abs(x) <= 2) 
#		y = - sign(x) * 3 * 0.25 * (2 - abs(x))^2;
#	endif
#	if (abs(x) <= 1)
#		y = y + sign(x) * 3 * (1 - abs(x))^2;
#	endif
endfunction

function y = ddbasis (x)
	y = 0;
	global dx;
	x = x / dx;
	if (abs(x) < 2)
		if (x < 0)
			y = 6 * 0.25 * (2 + x);
		else
			y = 6 * 0.25 * (2 - x);
		endif
		if (abs(x) < 1)
			if (x < 0)
				y = y - 6 * (1 + x);
			else
				y = y - 6 * (1 - x);
			endif
		endif
	endif
endfunction

function y = dddbasis (x)
	y = 0;
	global dx;
	x = x / dx;
	if (abs(x) < 2)
		if (x < 0)
			y = 6 * 0.25;
		else
			y = -6 * 0.25;
		endif
		if (abs(x) < 1)
			if (x < 0)
				y = y - 6;
			else
				y = y + 6;
			endif
		endif
	endif
endfunction


# m is the distance between the two dbasis functions
function y = qm (x, m)
	y = dddbasis(x) * dddbasis(x - m);
endfunction

function y = q0 (x)
	y = qm (x, 0);
endfunction

function y = q1 (x)
	y = qm (x, 1);
endfunction

function y = q2 (x)
	y = qm (x, 2);
endfunction

function y = q3 (x)
	y = qm (x, 3);
endfunction

x0 = -3 * dx;
xm = 3 * dx;
x = ( linspace (x0,xm,((xm-x0)/0.05)) );
#x = (-3:0.05:3);

for i = [ 1 : columns(x) ]
y(i) = basis(x(i));
dy(i) = dbasis(x(i));
ddy(i) = ddbasis(x(i));
dddy(i) = dddbasis(x(i));
dy2(1,i) = qm(x(i), 0);
dy2(2,i) = qm(x(i), 1);
dy2(3,i) = qm(x(i), 2);
dy2(4,i) = qm(x(i), 3);
endfor

plot (x,y,x,dy,x,ddy,x,dddy);
#(1,:),x,dy2(2,:),x,dy2(3,:),x,dy2(4,:));

printf ("dx = %f\n", dx);

printf ("For derivative constraint, where the row is the distance,\n");
printf ("and column is the X range, -2:-1, -1:0, 0:1, 1:2.\n");

#part(4,4);
for m = [ 1 : 4 ]
	x = m - 3;
	part(1,m) = quad("q0", x, x+1);
	part(2,m) = quad("q1", x, x+1);
	part(3,m) = quad("q2", x, x+1);
	part(4,m) = quad("q3", x, x+1);
endfor

part

#printf ("Integral of y(x): %f\n", quad("basis", x0, xm));
#printf ("Integral of dy(x): %f\n", quad("dbasis", x0, xm));
#printf ("Integral of q0(x): %f\n", quad("q0", x0, xm));
#printf ("Integral of q1(x): %f\n", quad("q1", x0, xm));
#printf ("Integral of q2(x): %f\n", quad("q2", x0, xm));
#printf ("Integral of q3(x): %f\n", quad("q3", x0, xm));


