#include "master.h"

class mglEigenVec : public mglDataA
{
public:
  long nx;
  long ny;
  long nz;
  Eigen::SparseVector<float>* a;

  inline mglEigenVec(){}

  inline mglEigenVec(Plate* plate, Eigen::SparseVector<float>* d)
  { nx=plate->x; ny=plate->y; nz=plate->z; a = d; }

  virtual ~mglEigenVec()  { if(a)  delete a; }

  inline long GetNx() const { return nx; }
  inline long GetNy() const { return ny; }
  inline long GetNz() const { return nz; }

  inline mreal Maximal() const  { return (a->coeffs()).maxCoeff(); }
  inline mreal Minimal() const  { return /*(a->coeffs()).minCoeff()*/0; }

protected:
  inline mreal v(long i,long j=0,long k=0) const
  { return a->coeff(i+nx*(j+ny*k)); }

  inline mreal vthr(long i) const
  { return a->coeff(i); }

  inline mreal dvx(long i,long j=0,long k=0) const
  { long i0=i+nx*(j+ny*k);
    float res=i>0? (i<nx-1? (a->coeff(i0+1)-a->coeff(i0-1))/2.:a->coeff(i0)-a->coeff(i0-1)) : a->coeff(i0+1)-a->coeff(i0);
    return res;  }

  inline mreal dvy(long i,long j=0,long k=0) const
  { long i0=i+nx*(j+ny*k);
    float res=j>0? (j<ny-1? (a->coeff(i0+nx)-a->coeff(i0-nx))/2.:a->coeff(i0)-a->coeff(i0-nx)) : a->coeff(i0+nx)-a->coeff(i0);
    return res;  }

  inline mreal dvz(long i,long j=0,long k=0) const
  { long i0=i+nx*(j+ny*k), n=nx*ny;
    float res=k>0? (k<nz-1? (a->coeff(i0+n)-a->coeff(i0-n))/2.:a->coeff(i0)-a->coeff(i0-n)) : a->coeff(i0+n)-a->coeff(i0);
    return res;  }

  inline mreal valueD(mreal x,mreal y=0,mreal z=0,mreal *dx=0,mreal *dy=0,mreal *dz=0) const
  { return 0;  }

  inline mreal value(mreal x,mreal y=0,mreal z=0) const
  { return 0;  }

};


class Foo: public mglDraw{
public:
    int Draw(mglGraph *gr);
    mglEigenVec data;
    Foo(Plate** p, Eigen::SparseVector<float>* d)
        : data(*p,d){}
};
int Foo::Draw(mglGraph *gr)
{
  gr->Title("MathGL Demo");
  gr->Rotate(60,40);
  gr->Alpha(true);
  //gr->Box();
  gr->Cloud(data,".wyrRk");
  return 0;
}


int run_mathgl_test(Plate** p, Eigen::SparseVector<float> *d) {
Foo foo(p,d);
mglFLTK gr(&foo, "MathGL run test");
return gr.Run();
}
