/**
 * @file main1.cpp
 * @author Auralius Manurung
 * @date 27 Jan 2017
 * 
 * @brief This is to test the class FX.
 *
 */

#include "fx.h"


colvec foo(colvec &x, colvec &a)
{
    colvec ret(2);
    double y0 = a.at(0)*x.at(0)*x.at(1)*x.at(2);
    double y1 = a.at(1)*x.at(0)*x.at(1)*x.at(2);
    ret << y0 << y1;
    return ret;
}

colvec foo2(colvec &x)
{
    colvec ret(2);
    double y0 = x.at(0)*x.at(1)*x.at(2);
    double y1 = x.at(0)*x.at(1)*x.at(2);
    ret << y0 << y1;
    return ret;
}

/////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    FX fxk(&foo);
    colvec x(1);
    colvec a(1);
    x << 1.0 << 2.0 << 3.0;
    a << 2.0 << 3.0;

    colvec y = fxk.SolveAt(x, a);
    mat jac = fxk.JacobianAt(x, a);
    mat hess1 = fxk.HessianAt(0, x, a);
    mat hess2 = fxk.HessianAt(1, x, a);

    y.print("y");
    jac.print("jac");
    hess1.print("hess1");
    hess2.print("hess2");
    
    std:cout << "-----------------------------------\n";
    
    FX fx(&foo2);

    y = fx.SolveAt(x);
    jac = fx.JacobianAt(x);
    hess1 = fx.HessianAt(0, x);
    hess2 = fx.HessianAt(1, x);

    y.print("y");
    jac.print("jac");
    hess1.print("hess1");
    hess2.print("hess2");

    return 0;
}