#ifndef NEW_PDF_MODELS
#define NEW_PDF_MODELS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class InvPow : public RooAbsPdf {
public:
    InvPow(const char *name, const char *title, RooAbsReal& _mgg,
           RooAbsReal& _p0, RooAbsReal& _p1, RooAbsReal& _p2);
    InvPow(const InvPow& other, const char* name=0);
    TObject* clone(const char* newname) const override { return new InvPow(*this,newname); }
    
protected:
    RooRealProxy mgg, p0, p1, p2;
    Double_t evaluate() const override;
};

class InvPowLin : public RooAbsPdf {
public:
    InvPowLin(const char *name, const char *title, RooAbsReal& _mgg,
             RooAbsReal& _p0, RooAbsReal& _p1, RooAbsReal& _p2, RooAbsReal& _p3);
    InvPowLin(const InvPowLin& other, const char* name=0);
    TObject* clone(const char* newname) const override { return new InvPowLin(*this,newname); }
    
protected:
    RooRealProxy mgg, p0, p1, p2, p3;
    Double_t evaluate() const override;
};

class Expow : public RooAbsPdf {
public:
    Expow(const char *name, const char *title, RooAbsReal& _mgg,
          RooAbsReal& _p0, RooAbsReal& _p1, RooAbsReal& _p2);
    Expow(const Expow& other, const char* name=0);
    TObject* clone(const char* newname) const override { return new Expow(*this,newname); }
    
protected:
    RooRealProxy mgg, p0, p1, p2;
    Double_t evaluate() const override;
};

class Dijet : public RooAbsPdf {
public:
    Dijet(const char *name, const char *title, RooAbsReal& _mgg,
          RooAbsReal& _p1, RooAbsReal& _p2);
    Dijet(const Dijet& other, const char* name=0);
    TObject* clone(const char* newname) const override { return new Dijet(*this,newname); }
    
protected:
    RooRealProxy mgg, p1, p2;
    Double_t evaluate() const override;
};

#endif
