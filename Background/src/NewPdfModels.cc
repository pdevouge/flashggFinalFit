#include "../interface/NewPdfModels.h"
#include <cmath>
#include <math.h>
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"


// InvPow Implementation
InvPow::InvPow(const char *name, const char *title, RooAbsReal& _mgg,
               RooAbsReal& _p1, RooAbsReal& _p2) :
    RooAbsPdf(name, title),
    mgg("mgg", "mgg", this, _mgg),
    p1("p1", "p1", this, _p1),
    p2("p2", "p2", this, _p2) {}

Double_t InvPow::evaluate() const {
    return pow(1 + mgg*p1, p2);
}

// InvPowLin Implementation
InvPowLin::InvPowLin(const char *name, const char *title, RooAbsReal& _mgg,
                    RooAbsReal& _p1, RooAbsReal& _p2, RooAbsReal& _p3) :
    RooAbsPdf(name, title),
    mgg("mgg", "mgg", this, _mgg),
    p1("p1", "p1", this, _p1),
    p2("p2", "p2", this, _p2),
    p3("p3", "p3", this, _p3) {}

Double_t InvPowLin::evaluate() const {
    return pow(1 + mgg*p1, p2 + p3*mgg);
}

// Expow Implementation
Expow::Expow(const char *name, const char *title, RooAbsReal& _mgg,
            RooAbsReal& _p1, RooAbsReal& _p2) :
    RooAbsPdf(name, title),
    mgg("mgg", "mgg", this, _mgg),
    p1("p1", "p1", this, _p1),
    p2("p2", "p2", this, _p2) {}

Double_t Expow::evaluate() const {
    return exp(p1*mgg) * pow(mgg, p2);
}

// Dijet Implementation
Dijet::Dijet(const char *name, const char *title, RooAbsReal& _mgg,
            RooAbsReal& _p1, RooAbsReal& _p2) :
    RooAbsPdf(name, title),
    mgg("mgg", "mgg", this, _mgg),
    p1("p1", "p1", this, _p1),
    p2("p2", "p2", this, _p2) {}

Double_t Dijet::evaluate() const {
    return pow(mgg, p1 + p2*log(mgg));
}

// InvPow copy constructor
InvPow::InvPow(const InvPow& other, const char* name) :
    RooAbsPdf(other, name),
    mgg("mgg", this, other.mgg),
    p1("p1", this, other.p1),
    p2("p2", this, other.p2) {}

// InvPowLin copy constructor
InvPowLin::InvPowLin(const InvPowLin& other, const char* name) :
    RooAbsPdf(other, name),
    mgg("mgg", this, other.mgg),
    p1("p1", this, other.p1),
    p2("p2", this, other.p2),
    p3("p3", this, other.p3) {}

// Expow copy constructor
Expow::Expow(const Expow& other, const char* name) :
    RooAbsPdf(other, name),
    mgg("mgg", this, other.mgg),
    p1("p1", this, other.p1),
    p2("p2", this, other.p2) {}

// Dijet copy constructor
Dijet::Dijet(const Dijet& other, const char* name) :
    RooAbsPdf(other, name),
    mgg("mgg", this, other.mgg),
    p1("p1", this, other.p1),
    p2("p2", this, other.p2) {}
