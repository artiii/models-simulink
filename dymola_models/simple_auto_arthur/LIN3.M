Tmax = 10;

[x0,n_X0]=initload; [p,n_par]=parload; p1=p; x0z=x0;

[s1,n1] = dymosim([0,Tmax,0,5000,1.e-6,8], x0, p);
Nmax = size(x0,1); x0=(s1(size(s1,1),2:Nmax+1)).';

dyn_in_create([0,10,0,500,1.e-6,8], x0, p);
eval( ['! C:\wp\Dymola\Mfiles\traj\alist.exe -a dsin.mat dsin.txt'] );

lin9

