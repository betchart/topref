#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "Math/BoostZ.h"
#include "Math/RotationZ.h"
#ifdef __CINT__
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;
#pragma link C++ function ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >::operator+(LV);
#pragma link C++ function ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >::operator-(LV);
#pragma link C++ function ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >::Dot(LV);
#pragma link C++ function ROOT::Math::BoostZ::operator()(LV);
#pragma link C++ function ROOT::Math::Boost::operator()(LV);
#pragma link C++ function ROOT::Math::RotationZ::operator()(LV);
#endif
