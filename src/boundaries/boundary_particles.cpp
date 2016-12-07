void part_periodic(double *x, double *v, const double x0, const double L){
	if(*x > L)
		*x-=L;
	else if(*x < x0)
		*x+=L;


}

void part_reflecting(double *x, double *v, const double x0, const double L){
	if(*x > L){
		*x = 2.0*L - *x;
		*v = -*v;
	}
	else if(*x < x0)
	{
		*x = 2.0*x0 - *x;
		*v = -*v;
	}
}

void part_absorbing(double *x, double *v, const double x0, const double L){}
