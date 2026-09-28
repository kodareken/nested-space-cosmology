// Fourth-order unitary transport for independent Dirac angular channels.
// This accelerates the same two-exponential method used by the Python owner.
#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

extern "C" int angular_transport(int nm, int nk, int nt, const double* masses,
    const double* momentum, int ratio_mode, const double* dt,
    const double* w1, const double* wp1, const double* w2, const double* wp2,
    double* positive, double* negative, int workers) {
  if (nm < 1 || nk < 1 || nt < 0 || workers < 1) return 1;
  const int total = nm * nk;
  const double a = (3.0-2.0*std::sqrt(3.0))/12.0;
  const double b = (3.0+2.0*std::sqrt(3.0))/12.0;
  auto evolve = [&](int lo, int hi) {
    for (int index=lo; index<hi; ++index) {
      const double m=masses[index/nk];
      const double k=momentum[index%nk]*(ratio_mode ? m : 1.0);
      const double m2=m*m, k2=k*k;
      double pr=positive[2*index], pi=positive[2*index+1];
      double nr=negative[2*index], ni=negative[2*index+1];
      auto step = [&](double y, double z, double delta) {
        const double norm=std::sqrt(y*y+z*z), phase=delta*norm;
        const double c=std::cos(phase), sn=std::sin(phase)/norm;
        const double ys=y*sn, zs=z*sn;
        const double npr=c*pr+zs*pi-ys*nr;
        const double npi=c*pi-zs*pr-ys*ni;
        const double nnr=ys*pr+c*nr-zs*ni;
        const double nni=ys*pi+c*ni+zs*nr;
        pr=npr; pi=npi; nr=nnr; ni=nni;
      };
      for (int t=0; t<nt; ++t) {
        const double e1=m2*w1[t]+k2, e2=m2*w2[t]+k2;
        const double y1=-m*k*wp1[t]/(4.0*std::sqrt(w1[t])*e1);
        const double y2=-m*k*wp2[t]/(4.0*std::sqrt(w2[t])*e2);
        const double z1=-std::sqrt(e1)/w1[t], z2=-std::sqrt(e2)/w2[t];
        step(b*y1+a*y2,b*z1+a*z2,dt[t]);
        step(a*y1+b*y2,a*z1+b*z2,dt[t]);
      }
      positive[2*index]=pr;positive[2*index+1]=pi;
      negative[2*index]=nr;negative[2*index+1]=ni;
    }
  };
  workers=std::min(workers,total);
  std::vector<std::thread> threads;
  for (int i=0;i<workers;++i) threads.emplace_back(evolve,total*i/workers,total*(i+1)/workers);
  for (auto& worker:threads) worker.join();
  return 0;
}
