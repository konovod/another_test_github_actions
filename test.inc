void NAME(ft_bruteforce_)(
    dbcf_index n,
    const Type *src_real,const Type *src_imag,
          Type *dst_real,      Type *dst_imag,
    int inverse,
    Type scale)
{
    dbcf_index i,j;
    Type pi=CAST(Type,4.0)*NAME(atan)(CAST(Type,1.0));
    for(i=0;i<n;++i)
    {
        dst_real[i]=CAST(Type,0.0);
        dst_imag[i]=CAST(Type,0.0);
    }
    for(i=0;i<n;++i)
    {
        Type x=CAST(Type,0.0);
        Type y=CAST(Type,0.0);
        Type w=CAST(Type,i)/CAST(Type,n);
        for(j=0;j<n;++j)
        {
            Type a=CAST(Type,2.0)*pi*(w*CAST(Type,j));
            Type c=NAME(cos)(a);
            Type s=NAME(sin)(a);
            if(!inverse) s=-s;
            x=x+src_real[j]*c-src_imag[j]*s;
            y=y+src_real[j]*s+src_imag[j]*c;
        }
        dst_real[i]=x;
        dst_imag[i]=y;
    }
    for(i=0;i<n;++i)
    {
        dst_real[i]=dst_real[i]*scale;
        dst_imag[i]=dst_imag[i]*scale;
    }
}

static void NAME(generate_)(dbcf_index seed,dbcf_index n,Type *real,Type *imag)
{
    dbcf_index i,j;
    dbcf_index MAX=n;
    RNG rng;
    RNG_init(&rng,(unsigned)seed);
    for(i=0;i<MAX;++i)
    {
        const double m=1.0/4294967296.0;
        Type x=CAST(Type,0.0);
        Type y=CAST(Type,0.0);
        for(j=0;j<4;++j) x=(x*CAST(Type,m)+CAST(Type,((double)RNG_generate(&rng)*m-0.5)));
        for(j=0;j<4;++j) y=(y*CAST(Type,m)+CAST(Type,((double)RNG_generate(&rng)*m-0.5)));
        real[i]=CAST(Type,x);
        imag[i]=CAST(Type,y);
    }
}


static double NAME(test_time_)(dbcf_index n,const Type *src_real,const Type *src_imag,Type *dst_real,Type *dst_imag,int interleaved)
{
    dbcf_index i;
    double t=get_cpu_time();
    dbcf_index m=DBCF_POW2(17)/n;
    if(n&(n-1)) m/=10;
    if(m==0) m=1;
    if(interleaved) for(i=0;i<m;++i) NAME2(dbc_fft_,s)(n,src_real,src_real+1,2,2,dst_real,dst_real+1,2,2,CAST(Type,1.0));
    else            for(i=0;i<m;++i) NAME2(dbc_fft_,s)(n,src_real,src_imag  ,1,1,dst_real,dst_imag  ,1,1,CAST(Type,1.0));
    t=get_cpu_time()-t;
    return t/(double)m;
}

static void NAME(get_norms_)(dbcf_index n,const Type *xr,const Type *xi,const Type *yr,const Type *yi,double *RMS,double *Linf)
{
    dbcf_index i;
    double e2=0.0,einf=0.0;
    for(i=0;i<n;++i)
    {
        double dr=0.0;
        double di=0.0;
        double d2=0.0;
        dr=CAST(double,(xi[i]-yi[i]));
        di=CAST(double,(xr[i]-yr[i]));
        d2=dr*dr+di*di;
        if(d2>einf) einf=d2;
        e2+=d2;
    }
    einf=sqrt(einf);
    e2/=(double)n;
    e2=sqrt(e2);
    *RMS=e2;
    *Linf=einf;
}

static void NAME(bitreverseal_permutation_bruteforce_)(dbcf_index log2n,Type *src,dbcf_index src_stride,Type *dst,dbcf_index dst_stride)
{
    dbcf_index n=DBCF_POW2(log2n);
    dbcf_index i;
    for(i=0;i<n;++i)
    {
        dbcf_index j=bitreverse_table[i];
        dst[j*dst_stride]=src[i*src_stride];
    }
}

static void NAME(test_bitreversal_permutation_)()
{
    dbcf_index i;
    Type tmp[65536];
    Type *buf0=data.NAME(buf_);
    printf("        N |out-of-place|     inplace|  bruteforce\n");
    printf("----------+------------+------------+------------\n");
    for(i=0;((MAXB/sizeof(Type))>>i)>1;++i)
    {
        dbcf_index n=DBCF_POW2(i),m=DBCF_POW2(21)/n;
        dbcf_index j;
        Type *buf1=buf0+n;
        double t;
        if(m<1) m=1;
        for(j=0;j<n;++j) buf0[j]=CAST(Type,j);
        for(j=0;j<n;++j) bitreverse_table[j]=bitreverse_bruteforce(j,i);
        t=get_cpu_time();
        for(j=0;j<m;++j) NAME(dbcF_bitreversal_permutation_)(i,buf0,1,buf1,1,tmp);
        t=get_cpu_time()-t;
        t/=(double)m;
        printf("%10.0f|%12.2f",(double)n,1e9*t/(double)n);
        for(j=0;j<n;++j) if(buf1[j]!=CAST(Type,bitreverse_table[j])) {printf(" FAIL!\n");return;}
        t=get_cpu_time();
        for(j=0;j<m;++j) NAME(dbcF_bitreversal_permutation_)(i,buf0,1,buf0,1,tmp);
        t=get_cpu_time()-t;
        t/=(double)m;
        printf("|%12.2f",1e9*t/(double)n);
        for(j=0;j<n;++j) buf0[j]=CAST(Type,j);
        NAME(dbcF_bitreversal_permutation_)(i,buf0,1,buf0,1,tmp);
        for(j=0;j<n;++j) if(buf0[j]!=CAST(Type,bitreverse_table[j])) {printf(" FAIL!\n");return;}
        for(j=0;j<n;++j) buf0[j]=CAST(Type,j);
        t=get_cpu_time();
        for(j=0;j<m;++j) NAME(bitreverseal_permutation_bruteforce_)(i,buf0,1,buf1,1);
        t=get_cpu_time()-t;
        t/=(double)m;
        printf("|%12.2f\n",1e9*t/(double)n);
    }
}

void NAME(test_fft_)(dbcf_index maxn)
{
    dbcf_index i,j,a,b;
    dbcf_index MAX=MAXB/sizeof(Type)/8;
    Type *buf=data.NAME(buf_);
    if(maxn<MAX) MAX=maxn;
    printf("          |     %5.5s     |     FFT-bruteforce    |     X-IFFT(FFT(X))    \n",(use_mflops?"Speed":"Time"));
    printf("        N |  SoA  |  AoS  |    RMS    |    Linf   |    RMS    |    Linf   \n");
    printf("----------+-------+-------+-----------+-----------+-----------+-----------\n");
    for(i=0;MAX>>i;++i)
    {
        dbcf_index n=DBCF_POW2(i);
        dbcf_index m=5*n*i;
        double RMS,Linf;
        double t;
        if(m==0) m=1;
        NAME(generate_)(37,n,buf+0*n,buf+1*n);
        printf("%10.0f|",(double)DBCF_POW2(i));
        for(j=0;j<2;++j)
        {
            t=NAME(test_time_)(n,buf+0*n,buf+1*n,buf+4*n,buf+5*n,(int)j);
            if(use_mflops) t=(double)m/(1.0e+9*t);
            else           t=1.0e+9*t/(double)m;
            printf("%7.3f|",t);
        }
        NAME2(dbc_fft_,s) (n,buf+0*n,buf+1*n,1,1,buf+4*n,buf+5*n,1,1,CAST(Type,1.0));
        NAME2(dbc_ifft_,s)(n,buf+4*n,buf+5*n,1,1,buf+6*n,buf+7*n,1,1,CAST(Type,1.0)/CAST(Type,n));
        if(i<=10)
        {
            NAME(ft_bruteforce_)(n,buf+0*n,buf+1*n,buf+2*n,buf+3*n,0,CAST(Type,1.0));
            NAME(get_norms_)(n,buf+2*n,buf+3*n,buf+4*n,buf+5*n,&RMS,&Linf);
            printf(" %-10.3e|",RMS);
            printf(" %-10.3e|",Linf);
        }
        else printf(" %-10s| %-10s|","-","-");
        NAME(get_norms_)(n,buf+0*n,buf+1*n,buf+6*n,buf+7*n,&RMS,&Linf);
        printf(" %-10.3e|",RMS);
        printf(" %-10.3e",Linf);
        printf("\n");
    }
    for(a=5,b=8;a<MAX;a+=b,b+=a)
    {
        dbcf_index n=a;
        double m=5.0*(double)n*log((double)n)/log(2.0);
        double RMS,Linf;
        double t;
        NAME(generate_)(37,n,buf+0*n,buf+1*n);
        printf("%10.0f|",(double)n);
        for(j=0;j<2;++j)
        {
            t=NAME(test_time_)(n,buf+0*n,buf+1*n,buf+4*n,buf+5*n,(int)j);
            if(use_mflops) t=m/(1.0e+9*t);
            else           t=1.0e+9*t/m;
            printf("%7.3f|",t);
        }
        NAME2(dbc_fft_,s) (n,buf+0*n,buf+1*n,1,1,buf+4*n,buf+5*n,1,1,CAST(Type,1.0));
        NAME2(dbc_ifft_,s)(n,buf+4*n,buf+5*n,1,1,buf+6*n,buf+7*n,1,1,CAST(Type,1.0)/CAST(Type,n));
        if(a<=1024)
        {
            NAME(ft_bruteforce_)(n,buf+0*n,buf+1*n,buf+2*n,buf+3*n,0,CAST(Type,1.0));
            NAME(get_norms_)(n,buf+2*n,buf+3*n,buf+4*n,buf+5*n,&RMS,&Linf);
            printf(" %-10.3e|",RMS);
            printf(" %-10.3e|",Linf);
        }
        else printf(" %-10s| %-10s|","-","-");
        NAME(get_norms_)(n,buf+0*n,buf+1*n,buf+6*n,buf+7*n,&RMS,&Linf);
        printf(" %-10.3e|",RMS);
        printf(" %-10.3e",Linf);
        printf("\n");
    }
}
