#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;
using Eigen::Lower;

// [[Rcpp::export]]
MatrixXd GetGInt(const Map<MatrixXi> & X, const Map<MatrixXd> & muhat, const int n, const int p){
    int         j       = 0;
    double      scale   = 0.0;
    //double      phat    = 0.0;
    MatrixXd    Z(n, p);
    MatrixXd    One     = MatrixXd::Ones(p, 1);
    MatrixXd    X_j;
    for(j = 0; j < p; j++){
        X_j         = X.col(j).cast<double>();
        //phat        = X_j.sum();
        Z.col(j)    = X_j - One*muhat(j)*2.0;//X_j - One*phat;
        scale       = scale + muhat(j)*(1.0 - muhat(j));
    }
    MatrixXd AAt(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(Z));
    return(AAt/(2.0*scale));
}

// [[Rcpp::export]]
MatrixXd GetGInd(const Map<MatrixXd> & X, const Map<MatrixXd> & muhat, const int n, const int p){
    int         j       = 0;
    double      scale   = 0.0;
    //MatrixXd    Z(n, p);
    //MatrixXd    One     = MatrixXd::Ones(p, 1);
    //MatrixXd    muhat   = X.colwise().sum()/(2.0*n);
    for(j = 0; j < p; j++){
        //Z.col(j)    = X.col(j) - One*muhat(j);
        scale       = scale + muhat(j)*(1.0 - muhat(j));
    }
    MatrixXd AAt(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(X));
    return(AAt/(2.0*scale));
}

// [[Rcpp::export]]
MatrixXd GetGPool(const Map<MatrixXd> & X, const Map<MatrixXd> & muhat, const int n, const int p){
    int         j       = 0;
    double      scale   = 0.0;
    MatrixXd    Z(n, p);
    MatrixXd    One     = MatrixXd::Ones(p, 1);
    //MatrixXd    muhat   = X.colwise().sum()/(2.0*n);
    for(j = 0; j < p; j++){
        Z.col(j)    = X.col(j) - One*muhat(j);
        scale       = scale + muhat(j)*(1.0 - muhat(j));
    }
    MatrixXd AAt(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(Z));
    return(2.0*AAt/scale);
}