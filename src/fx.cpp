/**
* @file fx.cpp
* @author Auralius Manurung
* @date 28 Jan 2017
* @brief Implement a user defined mathematical function.
*/

#include "fx.h"

FX::FX(colvec(*f) (colvec &x, colvec &some_constants))
{
    // Epsilon for computing the Jacobian numerically
    Epsilon_ = 1e-5;

    // Remember the callback function
    Fxk_ = f;
}

FX::FX(colvec(*f) (colvec &x))
{
    // Epsilon for computing the Jacobian numerically
    Epsilon_ = 1e-5;

    // Remember the callback function
    Fx_ = f;
}

FX::~FX()
{
}

colvec FX::SolveAt(colvec &x, colvec &some_constants)
{
    colvec ret = Fxk_(x, some_constants);
    return ret;
}

colvec FX::SolveAt(colvec &x)
{
    colvec ret = Fx_(x);
    return ret;
}

mat FX::JacobianAt(colvec &x, colvec &some_constants)
{
    colvec f0 = Fxk_(x, some_constants);
    colvec fn;

    int c = x.size();
    int r = f0.size();
    mat jac;
    jac.resize(r, c);

    for (int j = 0; j < c; j++) {
        colvec x_eps = x;
        x_eps(j) = x_eps(j) + Epsilon_;
        fn = Fxk_(x_eps, some_constants);

        for (int i = 0; i < r; i++) {
            jac(i, j) = (fn(i) - f0(i)) / Epsilon_;
        }
    }

    //jac.print("Jac=");
    return jac;
}

mat FX::JacobianAt(colvec &x)
{
    colvec f0 = Fx_(x);
    colvec fn;

    int c = x.size();
    int r = f0.size();
    mat jac;
    jac.resize(r, c);

    for (int j = 0; j < c; j++) {
        colvec x_eps = x;
        x_eps(j) = x_eps(j) + Epsilon_;
        fn = Fx_(x_eps);

        for (int i = 0; i < r; i++) {
            jac(i, j) = (fn(i) - f0(i)) / Epsilon_;
        }
    }

    //jac.print("Jac=");
    return jac;
}

mat FX::HessianAt(int i, colvec &x, colvec &some_constants)
{
    int n = x.size();
    mat hess;
    hess.resize(n, n);

    for (int j = 0; j < n; j++) {
        colvec x0 = x;

        x0(j) = x0(j) + Epsilon_;
        mat hess_plus = JacobianAt(x0, some_constants);

        x0 = x;
        x0(j) = x0(j) - Epsilon_;
        mat hess_min = JacobianAt(x0, some_constants);

        rowvec delta = (hess_plus.row(i) - hess_min.row(i)) / (2 * Epsilon_);
        hess.row(j) = delta;
    }

    //hess.print("Hessian=");
    return hess;
}

mat FX::HessianAt(int i, colvec &x)
{
    int n = x.size();
    mat hess;
    hess.resize(n, n);
    
    for (int j = 0; j < n; j++) {
        colvec x0 = x;

        x0(j) = x0(j) + Epsilon_;
        mat hess_plus = JacobianAt(x0);

        x0 = x;
        x0(j) = x0(j) - Epsilon_;
        mat hess_min = JacobianAt(x0);

        rowvec delta = (hess_plus.row(i) - hess_min.row(i)) / (2 * Epsilon_);
        hess.row(j) = delta;
    }

    //hess.print("Hessian=");
    return hess;
}

void FX::SetEpsilon(double epsilon)
{
    Epsilon_ = epsilon;
}
