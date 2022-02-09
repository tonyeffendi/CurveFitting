// to test gaussian model curve fitting with 5 parameters.
// This example below is taken from the link below and modified by Tony Effendi  
// https://github.com/RainerKuemmerle/g2o/blob/master/g2o/examples/data_fitting/curve_fit.cpp.

// g2o - General Graph Optimization
// Copyright (C) 2012 R. Kümmerle
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//#include <Eigen/Core>
#include <iostream>
#include <sstream>
#include <string>

#include "g2o/stuff/sampler.h"
#include "g2o/stuff/command_args.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/base_vertex.h"
#include "g2o/core/base_unary_edge.h"
#include "g2o/solvers/dense/linear_solver_dense.h"

using namespace std;

//typedef Eigen::Matrix<double,5,1> MyParam;

/**
 * Params: A, x0, s_x, y0, s_y
 * \brief the params: z = A * exp(-(0.5*(x-x0).^2/sx^2)-(0.5*(y-y0).^2/sy^2));
 */
class VertexParams : public g2o::BaseVertex<5, Eigen::Matrix<double,5,1>>
//class VertexParams : public g2o::BaseVertex<3, Eigen::vector5>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    VertexParams()
    {
    }

    virtual bool read(std::istream& /*is*/)
    {
      cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
      return false;
    }

    virtual bool write(std::ostream& /*os*/) const
    {
      cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
      return false;
    }

    virtual void setToOriginImpl()
    {
      cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
      _estimate  << 0,0,0,0,0;
    }

    virtual void oplusImpl(const double* update)
    {
      //Eigen::Vector3d::ConstMapType v(update);
      _estimate += Eigen::Matrix<double,5,1>(update);
    }
};

/**
 * \brief measurement for a point on the curve
 *
 * Here the measurement is the point which is lies on the curve.
 * The error function computes the difference between the curve
 * and the point.
 */
class EdgePointOnCurve : public g2o::BaseUnaryEdge<1, Eigen::Vector3d, VertexParams>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgePointOnCurve()
    {
    }
    virtual bool read(std::istream& /*is*/)
    {
      cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
      return false;
    }
    virtual bool write(std::ostream& /*os*/) const
    {
      cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
      return false;
    }

    void computeError()
    {
      const VertexParams* params = static_cast<const VertexParams*>(vertex(0));
      const double& A = params->estimate()(0);
      const double& x0 = params->estimate()(1);
      const double& sx = params->estimate()(2);
      const double& y0 = params->estimate()(3);
      const double& sy = params->estimate()(4);

      //%z = A * exp(-(0.5*(x-x0).^2/sx^2)-(0.5*(y-y0).^2/sy^2));
      const double x = measurement()(0);
      const double y = measurement()(1);
      const double z = measurement()(2);

      double fval = A * exp(-(0.5*pow((x-x0),2.0)/(sx*sx))-(0.5*pow((y-y0),2.0)/(sy*sy)));
      _error(0) = fval - z;
    }
};

void ReadFile(const string& filename, vector<vector<double>>& res)
{
  string strdata;
  ifstream myfile(filename);
  while(getline(myfile,strdata)) {
    vector<double> vectmp;
    istringstream is;
    is.str(strdata);
    double tmp;
    while(is>>tmp)    
      vectmp.push_back(tmp);      
    res.push_back(vectmp);
  }
  myfile.close();
} 


int main(int argc, char** argv)
{
  std::cout<<" ...."<<std::endl;

  int numPoints;
  int maxIterations = 130;
  bool verbose = true;
  std::vector<int> gaugeList;
  string dumpFilename = "/home/emulti/workspace/MyCurveFitting02/gauge.txt";
  g2o::CommandArgs arg;
  

  cout<<"read file "<<endl;


  // load txt data file
  string file_loc = "/home/emulti/workspace/MyCurveFitting02/zn.txt";
  vector<vector<double>> res;
  ReadFile(file_loc,res);


  numPoints = res.size();

  // construct input matrix
  Eigen::Vector3d* points = new Eigen::Vector3d[res.size()];
  for (int i = 0; i < res.size(); ++i) {
    points[i].x() = res[i][0];
    points[i].y() = res[i][1];
    points[i].z() = res[i][2];
  }

  if (dumpFilename.size() > 0) {
    ofstream fout(dumpFilename.c_str());
    for (int i = 0; i < numPoints; ++i)
      fout << points[i].transpose() << endl;
  }

  // some handy typedefs
  typedef g2o::BlockSolver< g2o::BlockSolverTraits<Eigen::Dynamic, Eigen::Dynamic> >  MyBlockSolver;
  typedef g2o::LinearSolverDense<MyBlockSolver::PoseMatrixType> MyLinearSolver;

  // setup the solver
  g2o::SparseOptimizer optimizer;
  optimizer.setVerbose(true);
  MyLinearSolver* linearSolver = new MyLinearSolver();
  MyBlockSolver* solver_ptr = new MyBlockSolver(linearSolver);
  g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
  optimizer.setAlgorithm(solver);

  // build the optimization problem given the points
  // 1. add the parameter vertex
  VertexParams* params = new VertexParams();
  params->setId(0);
  //params->setEstimate(Eigen::Matrix<double,5,1>(1,1,1,1,1)); // some initial value for the params
  Eigen::Matrix<double,5,1> tmp;
  tmp << 1,1,1,3,1;
  params->setEstimate(tmp); // some initial value for the params

  optimizer.addVertex(params);
  // 2. add the points we measured to be on the curve
  for (int i = 0; i < numPoints; ++i) {
    EdgePointOnCurve* e = new EdgePointOnCurve;
    e->setInformation(Eigen::Matrix<double, 1, 1>::Identity());
    e->setVertex(0, params);
    e->setMeasurement(points[i]);
    optimizer.addEdge(e);
  }

  // perform the optimization
  optimizer.initializeOptimization();
  optimizer.setVerbose(verbose);
  optimizer.optimize(maxIterations);

  if (verbose)
    cout << endl;

  // print out the result
  cout << "Target curve" << endl;
  cout << "%z = A * exp(-(0.5*(x-x0).^2/sx^2)-(0.5*(y-y0).^2/sy^2));" << endl;
  cout << "Iterative least squares solution" << endl;
  cout << "A      = " << params->estimate()(0) << endl;
  cout << "x0     = " << params->estimate()(1) << endl;
  cout << "s_x    = " << params->estimate()(2) << endl;
  cout << "y0     = " << params->estimate()(3) << endl;
  cout << "s_y    = " << params->estimate()(4) << endl;
  cout << endl;

  // clean up
  delete[] points;

  return 0;
}
