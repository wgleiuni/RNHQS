#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <string>
#include <fstream>
#include <time.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define Nblock 16
#define NThreadPerBlock 64
#define NTotal Nblock*NThreadPerBlock
#define N_ 200
class RNHQS_ZB
{
    public:
        __device__ RNHQS_ZB();
        __device__ void initial(int,double,double,double,double,double,int,double*,double*,double*,double*,double*,double*,double*);
        __device__ void go(double *);
        __device__ void record(double *);
    private:
        __device__ void onestep();
        __device__ void getMean();
        __device__ void da(double,double*,double*,double*,double*);

        int tN_,numdt_,i_;
//        int N_;
        double h_,sigma_,ri_,omega_,phi_,t_;
        double *a1_,*a2_,*k_[4],*l_[4],*ta1_,*ta2_,*X_;
        double outMo_,outX_;
};

__device__
RNHQS_ZB::RNHQS_ZB()
{
}

__device__
void RNHQS_ZB::initial(int tN,double h,double sigma,double ri,double w,double phi,int numdt,double* X,double* a1,double* a2,double* ta1,double* ta2,double* k,double* l)
{
    tN_=tN;
    h_=h;
    sigma_=sigma;
    ri_=ri;
    omega_=w;
    phi_=phi;
    numdt_=numdt;
//    N_=200;
    i_=0;

    t_=0.0;

//    a1_=(double *)malloc(N_*sizeof(double));
//    a2_=(double *)malloc(N_*sizeof(double));
//    X_=(double *)malloc(N_*sizeof(double));
//    ta1_=(double *)malloc(N_*sizeof(double));
//    ta2_=(double *)malloc(N_*sizeof(double));
    a1_=a1;
    a2_=a2;
    X_=X;
    ta1_=ta1;
    ta2_=ta2;
    int i;
//    for (i=0;i<4;i++)
//    {
//        *(k_+i)=new double [N_];
//        *(l_+i)=new double [N_];
//    }

    for (i=0;i<4;i++)
    {
        *(k_+i)=k+N_*i;
        *(l_+i)=l+N_*i;
    }
    double w_d=6.0;
    for (i=0;i<N_;i++)
    {
        *(X_+i)=-N_/2.0+1.0/2.0+i;
        *(a1_+i)=exp(-pow(*(X_+i)/w_d,2.0))*cos(M_PI**(X_+i)/2.0);
        *(a2_+i)=exp(-pow(*(X_+i)/w_d,2.0))*sin(M_PI**(X_+i)/2.0);
    }
}

__device__ 
void RNHQS_ZB::go(double *_outMo)
{
    int ntime;

    for (ntime=0;ntime<numdt_;ntime++)
    {
        t_+=h_;
        RNHQS_ZB::onestep();
        if (ntime%1==0)
        {
            RNHQS_ZB::record(_outMo);
        }
    }
}

__device__ 
void RNHQS_ZB::onestep()
{
    int i,j;
    for (i=0;i<4;i++)
    {
        if (i==0)
        {
            RNHQS_ZB::da(t_,a1_,a2_,*(k_+i),*(l_+i));
        }
        else if (i==1 || i==2)
        {
            for (j=0;j<N_;j++)
            {
                *(ta1_+j)=*(a1_+j)+h_/2.0**(*(k_+i-1)+j);
                *(ta2_+j)=*(a2_+j)+h_/2.0**(*(l_+i-1)+j);
            }
            RNHQS_ZB::da(t_+h_/2.0,ta1_,ta2_,*(k_+i),*(l_+i));
        }
        else if (i==3)
        {
            for (j=0;j<N_;j++)
            {
                *(ta1_+j)=*(a1_+j)+h_**(*(k_+i-1)+j);
                *(ta2_+j)=*(a2_+j)+h_**(*(l_+i-1)+j);
            }
            RNHQS_ZB::da(t_+h_,ta1_,ta2_,*(k_+i),*(l_+i));
        }
    }
    for (j=0;j<N_;j++)
    {
        *(a1_+j)+=h_/6.0*(*(*k_+j)+2.0**(*(k_+1)+j)+2.0**(*(k_+2)+j)+*(*(k_+3)+j));
        *(a2_+j)+=h_/6.0*(*(*l_+j)+2.0**(*(l_+1)+j)+2.0**(*(l_+2)+j)+*(*(l_+3)+j));
    }
}

__device__ 
void RNHQS_ZB::da(double t,double *a1,double *a2,double *k,double *l)
{
    int i;
    double sigmai;

    sigmai=sigma_*ri_*sin(omega_*t+M_PI*phi_);

    for (i=1;i<N_-1;i++)
    {
        *(k+i)=-(*(a2+i+1)+*(a2+i-1))+pow(-1.0,i)*(sigma_**(a2+i)+sigmai**(a1+i));
        *(l+i)=(*(a1+i+1)+*(a1+i-1))+pow(-1.0,i+1)*(sigma_**(a1+i)-sigmai**(a2+i));
    }

    *(k+0)=-(*(a2+1)+*(a2+N_-1))+pow(-1.0,0)*(sigma_**(a2+0)+sigmai**(a1+0));
    *(l+0)=(*(a1+1)+*(a1+N_-1))+pow(-1.0,0+1)*(sigma_**(a1+0)-sigmai**(a2+0));

    *(k+N_-1)=-(*(a2+0)+*(a2+N_-2))+pow(-1.0,N_-1)*(sigma_**(a2+N_-1)+sigmai**(a1+N_-1));
    *(l+N_-1)=(*(a1+0)+*(a1+N_-2))+pow(-1.0,N_)*(sigma_**(a1+N_-1)-sigmai**(a2+N_-1));
}

__device__ 
void RNHQS_ZB::getMean()
{
    int i;
    outMo_=0.0;
    outX_=0.0;
    for (i=0;i<N_;i++)
    {
        outMo_+=pow(*(a1_+i),2.0)+pow(*(a2_+i),2.0);
        outX_+=*(X_+i)*(pow(*(a1_+i),2.0)+pow(*(a2_+i),2.0));
    }

//    if (wasExecuted_)
//   {
//        return;
//    }
//    else
//    {
//        double rest;
//        for (i=N_-10;i<N_-1;i++)
//        {
//            rest+=pow(*(a1_+i),2.0)+pow(*(a2_+i),2.0);
//        }

//        if (rest/outMo_>0.001)
//        {
 //           std::cout << tN_ << ": waveguide not long enough" << std::endl;
//        }
//        wasExecuted_ = true;
//    }
}

__device__ 
void RNHQS_ZB::record(double *_outMo)
{
    RNHQS_ZB::getMean();
    *(_outMo+i_)=outMo_;
    i_++;
}

__global__ void sRNHQS_ZB(int tN,double h,double sigma,double ri,double* omega,double phi,int numdt,double* _outMo,size_t pitch)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;

    RNHQS_ZB one;
//    __shared__ double a1[N_];
//    __shared__ double a2[N_];
//    __shared__ double X[N_];
//    __shared__ double ta1[N_];
//    __shared__ double ta2[N_];
//    __shared__ double k[4*N_];
//    __shared__ double l[4*N_];
    double a1[N_];
    double a2[N_];
    double X[N_];
    double ta1[N_];
    double ta2[N_];
    double k[4*N_];
    double l[4*N_];
    one.initial(tN,h,sigma,ri,*(omega+i),phi,numdt,X,a1,a2,ta1,ta2,k,l);
    one.go(_outMo+pitch/sizeof(double)*i);
}

int main (int argc, char *argv[])
{
    int tN,numdt,i,j;
    double h,phi,sigma,ri;
    double *omega;

    omega=new double [NTotal];
    for (i=0;i<NTotal;i++)
    {
        *(omega+i)=10.0/(NTotal-1.0)*i;
    }

//    size_t s1,s2,s3;
//    cudaDeviceSetLimit(cudaLimitStackSize,4*1024*sizeof(double));
    cudaError cE;
    cE=cudaDeviceSetLimit(cudaLimitMallocHeapSize,16*1024*1024*sizeof(double));
    if (cE!=0) std::cout << cE << std::endl;
//    cudaDeviceSetLimit(cudaLimitPrintfFifoSize,1024*1024*sizeof(double));
//    cudaDeviceGetLimit(&s1,cudaLimitStackSize);
//    cudaDeviceGetLimit(&s2,cudaLimitPrintfFifoSize);
//    cudaDeviceGetLimit(&s3,cudaLimitMallocHeapSize);
//    std::cout << s1/sizeof(double)/1024/1024 << "M " << s2/sizeof(double)/1024/1024 << " " << s3/sizeof(double)/1024/1024 << std::endl;
//    std::cout << "Stack size per thread = " << s1/1024 << "Kb" << std::endl;
//    std::cout << "IO in total = " << s2/1024 << "Kb" << std::endl;
//    std::cout << "IO per thread = " << s2/1024/1024 << "Kb" << std::endl;
//    std::cout << "Heap size in total = " << s3/1024 << "Kb" << std::endl;
//    std::cout << "Heap size in total = " << s3/1024/1024 << "Mb" << std::endl;
//    std::cout << "Heap size per thread = " << s3/1024/1024 << "Kb" << std::endl;
//    std::cout << "Total double per thread = " << s3/sizeof(double)/1024  << std::endl;

    if (argc>1)
    {
        
        if (argc==7)
        {
            tN=(int)atof(argv[1]);
            h=(double)atof(argv[2]);
            sigma=(double)atof(argv[3]);
            ri=(double)atof(argv[4]);
            phi=(double)atof(argv[5]);
            numdt=(int)atof(argv[6]);

            double* _omega;
            cE=cudaMalloc(&_omega,1024*sizeof(double));
            if (cE!=0) std::cout << cE << std::endl;
            cE=cudaMemcpy(_omega,omega,1024*sizeof(double),cudaMemcpyHostToDevice);
            if (cE!=0) std::cout << cE << std::endl;

            double *_outMo,*outMo;
            size_t pitch;
            cE=cudaMallocPitch(&_outMo,&pitch,sizeof(double)*numdt,NTotal);
            if (cE!=0) std::cout << cE << std::endl;
            outMo=new double [numdt*NTotal];
            sRNHQS_ZB<<<Nblock,NThreadPerBlock>>>(tN,h,sigma,ri,_omega,phi,numdt,_outMo,pitch);
            cE=cudaDeviceSynchronize();
            if (cE!=0) std::cout << cE << std::endl;
            std::cout << "CUDA computation complete. Start to copy data into host." << std::endl;
            cE=cudaMemcpy2D(outMo,numdt*sizeof(double),_outMo,pitch,numdt*sizeof(double),NTotal,cudaMemcpyDeviceToHost);
            if (cE!=0) std::cout << cE << std::endl;
            std::cout << "Copy finished. Start to write to file." << std::endl;

            char* filename[NTotal];
            std::ofstream out[NTotal];
            for (i=0;i<NTotal;i++)
            {
                filename[i]=new char [20];
                sprintf(filename[i],"ZB%d.txt",i+1);
                out[i].open(filename[i],std::ostream::out);
                for (j=0;j<numdt;j++)
                {
                    out[i] << *((double *)outMo+i*numdt+j) << std::endl;
                }
            }
        }
    }

    std::cout << "time used = " << clock()/CLOCKS_PER_SEC << std::endl;

    return 0;
}
