class param_set
{
public:
    double r; //growth speed
    double la,lb; // apic and basal zones length
    double ln; // space between branches
    double N; // maximal number of braches
    double angle; // angle between root and predecessor (rad)
    double max_len; // maximal root length
    double w1; // weight for angle change using gradient field
    param_set() {}
    param_set(double _r,double _la,double _lb,double _ln,double _N,double _a,double _ml,double _w1)
    {
	r=_r;la=_la;lb=_lb;ln=_ln;N=_N;angle=_a;max_len=_ml;w1=_w1;
    }
    void write(FILE *fi)
    {
	fprintf(fi,"%g %g %g %g %g %g %g %g\n",r,la,lb,ln,N,angle,max_len,w1);
    }
    void read(FILE *fi)
    {
	fscanf(fi,"%lg %lg %lg %lg %lg %lg %lg %lg\n",&r,&la,&lb,&ln,&N,&angle,&max_len,&w1);
    }
};
std::vector<param_set> param_sets;
class segment
{
public:
    double x[2][2];
    segment() {}
    segment(double x00,double x01,double x10,double x11)
    {
	x[0][0]=x00;x[0][1]=x01;x[1][0]=x10;x[1][1]=x11;
    }
    double length()
    {
	return sqrt((x[1][0]-x[0][0])*(x[1][0]-x[0][0])+(x[1][1]-x[0][1])*(x[1][1]-x[0][1]));
    }
    void dir(double &d0,double &d1)
    {
	d0=(x[1][0]-x[0][0])/length();
	d1=(x[1][1]-x[0][1])/length();
    }
    void print()
    {
	printf("%g %g %g %g\n",x[0][0],x[0][1],x[1][0],x[1][1]);
    }
    void write(FILE *fi)
    {
	fprintf(fi,"%g %g %g %g\n",x[0][0],x[0][1],x[1][0],x[1][1]);
    }
    void read(FILE *fi)
    {
	fscanf(fi,"%lg %lg %lg %lg\n",&x[0][0],&x[0][1],&x[1][0],&x[1][1]);
    }
};
void grad_func_0(double x,double y,double &d0,double &d1,void *a)
{
    d0=d1=0;
}
void grad_func_1(double x,double y,double &d0,double &d1,void *a)
{
    d0=sin(x/180.0)*cos(y/180.0)+0.25*(((rand()%1000)/10000.0)-0.05);
    d1=cos(x/180.0)*cos(y/180.0)+0.25*(((rand()%1000)/10000.0)-0.05);
}
class root;
class root
{
    param_set p;
    double level;
    std::vector<segment> segments;
    void (*grad_func)(double x,double y,double &d0,double &d1,void *);
    std::vector<root> branches;
    void *fp;
public:
    root(param_set _p,segment init,void (*_gf)(double,double,double &,double &,void *),double _l,void * _fp=NULL)
    {
	p=_p;level=_l;fp=_fp;
	grad_func=_gf;
	segments.push_back(init);
    }
    root()
    {
    }
    ~root()
    {
    }
    void print()
    {
	for (int i=0;i<segments.size();i++)
		segments[i].print();
	for (int i=0;i<branches.size();i++)
		branches[i].print();
    }
    void write(FILE *fi)
    {
	p.write(fi);
	fprintf(fi,"%g %d %d\n",level,segments.size(),branches.size());
	for (int i=0;i<segments.size();i++)
		segments[i].write(fi);
	for (int i=0;i<branches.size();i++)
		branches[i].write(fi);
    }
    void read(FILE *fi,void (*_gf)(double,double,double &,double &,void *),void * _fp=NULL)
    {
	fp=_fp;
	grad_func=_gf;
	p.read(fi);
	int ns,nb;
	fscanf(fi,"%lg %d %d\n",&level,&ns,&nb);
	for (int i=0;i<ns;i++)
	{
		segments.push_back(segment());
		segments[i].read(fi);
	}
	for (int i=0;i<nb;i++)
	{
		branches.push_back(root());
		branches[i].read(fi,_gf,_fp);
	}
    }
    double length()
    {
	double ret=0;
	for (int i=0;i<segments.size();i++)
	    ret+=segments[i].length();
	return ret;
    }
    int size()
    {
	int ret=segments.size();
	for (int i=0;i<branches.size();i++)
	    ret+=branches[i].size();
	return ret;
    }
    double bounds(int c,double _cw=0)
    {
	double w=0;
	double cw=((level==0)?segments[0].x[0][c]:_cw);
	for (int i=0;i<segments.size();i++)
	{
	    double lw=fabs(segments[i].x[1][c]-cw);
	    if (lw>w) w=lw;
	}
	for (int i=0;i<branches.size();i++)
	{
	    double bw=branches[i].bounds(c,cw);
	    if (bw>w) w=bw;
	}
	return w;
    }
    double *density(int nx,int ny,double *od=NULL,double _cwx=0,double _cwy=0,double _b0=0,double _b1=0)
    {
	double *da;
	if (od==NULL)
	{
	     da=new double[nx*ny];
	     memset(da,0, nx*ny*sizeof(double));
	}
	else 
	    da=od;
	double cwx=((level==0)?segments[0].x[0][0]:_cwx);
	double cwy=((level==0)?segments[0].x[0][1]:_cwy);
	double b0=((level==0)?bounds(0):_b0);
	double b1=((level==0)?bounds(1):_b1);
	if (b0==0.0) b0=p.r;
	if (b1==0.0) b1=p.r;
	double x0=cwx-b0;
	double y0=cwy;
	double lx=2*b0/nx;
	double ly=b1/ny;
	double step=lx;
	if (ly<lx) step=ly;
	for (int i=0;i<segments.size();i++)
	{
	    int nsteps=segments[i].length()/step;
	    int oi,oj,curi,curj;
	    double d0,d1;
	    if (nsteps>(nx+ny))
	    {
		nsteps=nx+ny;
		step=segments[i].length()/nsteps;
	    }
	    segments[i].dir(d0,d1);
	    oi=(segments[i].x[0][0]-x0)/lx;
	    oj=(segments[i].x[0][1]-y0)/ly;
	    if ((oi>=0)&&(oi<nx)&&(oj>=0)&&(oj<ny))
		da[oi*ny+oj]+=1.0;
	    for (int j=1;j<=nsteps;j++)
	    {
		curi=(segments[i].x[0][0]+d0*step*j-x0)/lx;
		curj=(segments[i].x[0][1]+d1*step*j-y0)/ly;
		if ((curi!=oi)&&(curj!=oj))
		{
		    if ((curi>=0)&&(curi<nx)&&(curj>=0)&&(curj<ny))
			da[curi*ny+curj]+=1.0;
		    oi=curi;
		    oj=curj;
		}
	    }
	    curi=(segments[i].x[1][0]-x0)/lx;
	    curj=(segments[i].x[1][1]-y0)/ly;
	    if ((curi>=0)&&(curi<nx)&&(curj>=0)&&(curj<ny))
	        if ((curi!=oi)&&(curj!=oj))
			da[curi*ny+curj]+=1.0;
	}
	for (int i=0;i<branches.size();i++)
	    branches[i].density(nx,ny,da,cwx,cwy,b0,b1);
	if (level==0)
	{
	    double sum=0;
	    for (int i=0;i<nx;i++)
		for (int j=0;j<ny;j++)
		    sum+=da[i*ny+j];
	    if (sum!=0)
	    for (int i=0;i<nx;i++)
		for (int j=0;j<ny;j++)
		    da[i*ny+j]/=sum;
	}
	return da;
    } 
    void growth(double dt) // root growth
    {
	if (length()+p.r*dt<p.max_len)
	{
	    double d0,d1;
	    double g0,g1;
	    double x,y;
	    x=segments[segments.size()-1].x[1][0];
	    y=segments[segments.size()-1].x[1][1];
	    segments[segments.size()-1].dir(d0,d1); //direction of the last segment
	    grad_func(x,y,g0,g1,fp); //gradient in the last point
	    d0+=p.w1*g0;
	    d1+=p.w1*g1;
	    segments.push_back(segment(x,y,x+d0*p.r*dt,y+d1*p.r*dt));
	}
	for (int i=0;i<branches.size();i++)
	    branches[i].growth(dt);
    }
    void branching(double dt) 
    {
	if (level<param_sets.size()-1)
	if (branches.size()<p.N)
	if (length()>=(p.la+p.lb+(branches.size()+1)*p.ln)) // branching zone must be already created and there is a space for the next branch
	{
	    double new_angle=p.angle*2*((rand()%10000)-5000)/10000.0;
	    double xpos=p.lb+(length()-p.la-p.lb)*((rand()%10000)/10000.0);
	    double d0,d1;
	    double x,y,cur_l=0;
	    int found=0;
	    // initial point
	    int segm_i=0;
	    for (int i=0;i<segments.size();i++)
	    {
		if ((xpos>=cur_l)&&(xpos<cur_l+segments[i].length())) // initial point is within the segment
		{
		    double k=(xpos-cur_l)/segments[i].length();
		    x=segments[i].x[0][0]+(segments[i].x[1][0]-segments[i].x[0][0])*k;
		    y=segments[i].x[0][1]+(segments[i].x[1][1]-segments[i].x[0][1])*k;
		    segm_i=i;
		    found=1;
		    break;
		}
		else
		    cur_l+=segments[i].length();
	    }
	    if (found)
	    {
		// direction 
		segments[segm_i].dir(d0,d1); //direction of the segment
		d0+=sin(new_angle*M_PI/180.0);
		d1+=cos(new_angle*M_PI/180.0);
		double l=sqrt(d0*d0+d1*d1);
		if (l!=0)
		{
		    d0/=l;
		    d1/=l;
		}
		branches.push_back(root(param_sets[level+1],segment(x,y,x+d0*p.r*dt,y+d1*p.r*dt),grad_func,level+1,fp));
	    }
	}
	for (int i=0;i<branches.size();i++)
	    branches[i].branching(dt);
    }
};
void init_paramsets()
{
//    param_sets.push_back(param_set(2,15.7,0.7,0.7,50,63,30,0.5));
//    param_sets.push_back(param_set(6.4,2.7,0.7,0.7,7,63,15,0.5));
//    param_sets.push_back(param_set(6.4,2.7,0.7,0.7,3,63,10,0.5));
    param_sets.push_back(param_set(0.75/(30*86400),0.05,0.01,0.01,50,63,1.0,0.5));
    param_sets.push_back(param_set(7.5/(30*86400),0.01,0.01,0.01,7,63,0.15,0.5));
    param_sets.push_back(param_set(7.5/(30*86400),0.01,0.01,0.01,3,63,0.10,0.5));
}
