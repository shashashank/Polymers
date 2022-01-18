#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
#include<fstream>
#include<iostream>
#include<complex>

clock_t start = clock();

using namespace  std;

int N;
int T=500;

int alpha=6;
//int M=66;
int t;

double pi=4*atan(1.0);
double dt=1;
double rcut=pow(2,1.0/alpha);
//double BOXL=M*rcut;
double BOXL=46.0209;
double psi6;
double rcutsq=rcut*rcut;

vector<double> x;
vector<double> y;

vector<vector<double>> nbrs;

typedef complex<double> dcomp;
dcomp a,b,iota,psi;
vector<dcomp> q6;
vector<dcomp> qq;


void read_data(ifstream& list,ofstream& read)
{
    string s;
    char c;
    double data1,data2,data3;

    list>>N;
    list>>s;
    x.resize(N);
    y.resize(N);
    q6.resize(N);
    qq.resize(N);
    read<<N<<endl;
    read<<s<<endl;
    for(int i=0; i<N; i++)
    {
        list>>c>>data1>>data2>>data3;
        x[i]=data1;
        y[i]=data2;
        read<<i<<"   "<<x[i]<<"   "<<data1<<"   "<<y[i]<<"  "<<data2<<"  "<<endl;
    }
}
void neighbor(double rcutsq,int N)
{
    nbrs.resize(N);
    for(int i=0; i<N; i++)
    {
        nbrs[i].resize(0);
    }
    for(int i=0; i<N-1; i++)
    {
        for(int j=i+1; j<N; j++)
        {
            double  xij  = x[i] - x[j];
            double  yij  = y[i] - y[j];
            xij=xij-BOXL*rint(xij/BOXL);
            yij=yij-BOXL*rint(yij/BOXL);

            double rijsq = xij * xij + yij * yij ;
            double  rijq = yij * yij ;
            if ( rijsq < rcutsq )
            {
                nbrs[i].push_back(j);
                nbrs[j].push_back(i);
            }
        }
    }
}

void q6_(int N)
{
    /****complex no. defined****/
    b=-1;
    iota=sqrt(b);
double nb;
    for(int i=0; i<N; i++)
    {
        q6[i]=(0.0,0.0);
    }
    double th_ij,th_i0;
    for(int i=0; i<N; i++)
    {
	 nb=nbrs[i].size();
        double er=0.0,ep=0.0;
        for(int j=0; j<nb; j++)
        {
            th_ij=0.0;
            int m=nbrs[i][j];
            double  xij  = x[i] - x[m];
            double  yij  = y[i] - y[m];

            double xi0=x[i]-x[nbrs[i][0]];
            double yi0=y[i]-y[nbrs[i][0]];
            xij=xij-BOXL*rint(xij/BOXL);
            yij=yij-BOXL*rint(yij/BOXL);

            if(xij==0.0)
                th_ij=2*acos(0.0);
		else
                th_ij=atan(yij/xij);
            if(xi0==0.0)
                th_i0=2*acos(0.0);
		else
                th_i0=atan(yi0/xi0);

                th_ij=th_ij -th_i0;


            //    q6[i]+=exp(6*th_ij*iota);
                q6[i]+=(1.0/nb)*exp(6*th_ij*iota);
        }
    }
}

void psi6_(int N)
{
    psi=(0,0);
    int m=0;
    for(int i=0; i<N; i++)
    {
        //if(nbrs[i].size()==6)
        {
            psi+=q6[i];
            m++;
        }
    }
    psi=(1.0/m)*psi;
    psi6=sqrt(psi.real()*psi.real()+psi.imag()*psi.imag());
}
void write_global_parameter(ostream& out)
{

    out.setf(ios::fixed,ios::floatfield);
    out.precision(5);
    out<<t*dt<<" "<<psi6<<" "<<endl;
    out.flush();

}
void write_local_parameter(ostream& out)
{

    out.setf(ios::fixed,ios::floatfield);
    out.precision(5);
    out<<"T:"<<t*dt<<"\n"<<endl;
    for(int i=0; i<N; i++)
    {
        out<<x[i]<<" "<<y[i]<<" "<<q6[i].real()<<" "<<q6[i].imag()<<"  "<<" "<<nbrs[i].size()<<" "<<i<<" "<<endl;
    }
    out.flush();

}

int main()
{
    ifstream list;
    ofstream prt("reading.txt");
    ofstream og("global_order_parameter.txt");
    ofstream ol("local_order_parameter.txt");

    list.open("vmd_data.xyz");
     //list.open("../cluster.xyz");
    for(t=0; t<T; t++)
    {
        read_data(list,prt);
        cout<<t*dt<<endl;
        neighbor(rcutsq,N);
        q6_(N);
        psi6_(N);
        write_global_parameter(og);
        //write_local_parameter(ol);
    }
    list.close();
    cout<<"Time elapsed: "<<( (double(clock()) - start) / CLOCKS_PER_SEC)<<endl;

}
