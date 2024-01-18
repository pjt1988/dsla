#include "bcsrmatrix.h"
#include "densematrix.h"
#include "imatrix.h"
#include "settings.h"

#include "miscfunctions.h"
#include "timing.h"

#include <iostream>
#include <numeric>

#include <set>
#include <vector>

using namespace DSLA;

int main(){

  auto timer = Timer::Instance();
  timer->start();
  BCSRMatrix t(10,10);
  BCSRMatrix t2(10,10,4);

  std::cout << "mat type " << t.getTypeAsString() << std::endl;

  auto settings = Settings::Instance();
  settings->setSparsityThresh(1e-6);

  std::cout << "default sparsity " << settings->getSparsityThresh() << std::endl;

#if 0
  //DenseMatrix d(10,10);
  auto* d = new DenseMatrix(4,4);
  generateRandomMatrix(4,4,-1.0,10.0,d->getBuffer(),true);

  d->print();
  timer->stop(Timings::BCSR_COPY);


  timer->printTimes();

  auto timer2 = Timer::Instance();
  timer2->start();

  t.generateBlocking(4,d->getBuffer());
  std::cout << "first " << std::endl;
  t.print();
  t.scale(2.0);
  std::cout << "scale 2 " << std::endl;
  t.print();
  t.scale(0.1);
  std::cout << "scale 0.1 " << std::endl;
  t.print();
  std::cout << "scale 0.00001 " << std::endl;
  t.scale(0.000001);

  t.print();

  std::cout << "bcsr norm " << t.norm() << std::endl;

  //t2.copy(*d);

  t.write("test.dat");
  t.clear();
  t.print();
  t.read("test.dat");
  t.print();

  timer2->stop(Timings::BCSR_ADD);
#endif  

  auto* d2 = new DenseMatrix(15000,15000);
  generateRandomMatrix(15000,15000,-100,100,d2->getBuffer(),true);
  //auto* d3 = new DenseMatrix(20000,20000);
  //generateRandomMatrix(20000,20000,0.0,10,d3->getBuffer(),true);

#if 0
  timer->resetAll();
  timer->start();
  for(int i=0;i<100;++i){
    std::vector<double> vec;
    vec.reserve(200*200);
    for(int j=0;j<200*200;++j){
      vec.push_back(d2->getBuffer()[j]);
    }
    std::sort(vec.begin(), vec.end());
  }
  timer->stop(Timings::MISC);
  std::cout << "vec insert + sort \n";
  timer->printTimes();
  timer->resetAll();
  timer->start();
  for(int i=0;i<100;++i){
    std::set<double> s;
    for(int j=0;j<200*200;++j){
      s.insert(d2->getBuffer()[j]);
    }
  }
  timer->stop(Timings::MISC);
  std::cout << "with set \n";
  timer->printTimes();
#endif


  BCSRMatrix c1(15000,15000);
  c1.generateBlocking(250,d2->getBuffer());
  c1.occupancy();

#if 0
  BCSRMatrix c2(20000,20000);
  c2.generateBlocking(250,d2->getBuffer());


  c1.occupancy();
  c2.occupancy();

  BCSRMatrix s1(20000,20000);
  s1.generateBlocking(250,d3->getBuffer());
  s1.occupancy();
  timer->resetAll();
  timer->printTimes();
  timer->start();
  c1.copyFromBCSR(s1);
  timer->stop(Timings::BCSR_COPY);
  c1.occupancy();
  timer->start();
  c2.copyFromBCSR2(s1);
  timer2->stop(Timings::BCSR_ADD);
  c2.occupancy();

  auto* addTest = new DenseMatrix(200,200);
  generateRandomMatrix(200,200,-10,10,addTest->getBuffer(),false);
  auto* addTest2 = new DenseMatrix(200,200);
  generateRandomMatrix(200,200,-10,10,addTest2->getBuffer(),false);
  DenseMatrix addTest3(*addTest);
  DenseMatrix addTest4(*addTest);
  DenseMatrix addTest5(*addTest);

  const auto num = 15000;
  std::vector<double> res(num);

  //separate add + scale
  timer->resetAll();
  timer->start();
  for(auto i=0;i<num;++i){
    _dadd(addTest->getBuffer(), addTest2->getBuffer(), 200*200);
    res[i] = fNorm(addTest->getBuffer(), 200*200);
  }
  timer->stop(Timings::MISC);
  std::cout << "Final summed up res " << std::reduce(res.begin(),res.end()) << std::endl; 
  std::fill(res.begin(),res.end(), 0.0);

  //combined test
  timer->start();
  for(auto i=0;i<num;++i){
    res[i] = _daddS(addTest3.getBuffer(), addTest2->getBuffer(), 200*200);
  }
  timer->stop(Timings::BCSR_ADD);
  std::cout << "Combo Final summed up res " << std::reduce(res.begin(),res.end()) << std::endl;
  timer->printTimes();
  timer->resetAll();

  timer->start();
  for(auto i=0;i<num;++i){
    res[i] = _daddS(addTest4.getBuffer(), addTest2->getBuffer(), 200*200);
  }
  timer->stop(Timings::BCSR_ADD);
  std::cout << "2 Combo Final summed up res " << std::reduce(res.begin(),res.end()) << std::endl;

  timer->start();
  for(auto i=0;i<num;++i){
    auto dadds = [&](double* __restrict__ a, double* __restrict__ b, const size_t dim){
      double sum =0;
      #pragma omp parallel for reduction(+:sum)
      for(auto i=0ul;i<dim;++i){
        a[i] += b[i];
        sum += a[i]*a[i];
      }
      return sum;
    };
    res[i] = dadds(addTest5.getBuffer(), addTest2->getBuffer(),200*200);

  }
  timer->stop(Timings::MISC);
  std::cout << "lambda Final summed up res " << std::reduce(res.begin(),res.end()) << std::endl;

#endif

double* buf = new double[81];
generateRandomMatrix(9,9,1.0,5.0,buf);
for(int i=0;i<3;++i){
  for(int j=6;j<9;++j){
    buf[i*9+j] = 0.0;
  }
}
BCSRMatrix tmat(9,9);
tmat.generateBlocking(3,buf);
tmat.print();

printf("new transpose\n\n");
tmat.transpose_in_place();
tmat.print();


#if 0
double sum = 0.0;
timer->start();
for(auto i=0;i<50;++i){
  c1.transpose_in_place();
  sum += c1.norm();
}
timer->stop(Timings::BCSR_TRANSPOSE);
std::cout << "summed norm is " << sum << std::endl;
timer->printTimes();

timer->reset(Timings::BCSR_TRANSPOSE);

sum = 0.0;
timer->start();
for(auto i=0;i<50;++i){
  c1.transpose_in_place_old();
  sum += c1.norm();
}
timer->stop(Timings::BCSR_TRANSPOSE);
std::cout << "summed norm is " << sum << std::endl;
timer->printTimes();
#endif


double sum = 0.0;
timer->start();
for(auto i=0;i<10000;++i){
  sum += c1.trace();
}
timer->stop(Timings::BCSR_TRACE);
std::cout << "NEW summed trace is " << sum << std::endl;
timer->printTimes();



}
