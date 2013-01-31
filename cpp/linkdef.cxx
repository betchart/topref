#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "Math/BoostZ.h"
#include "Math/RotationZ.h"
#include "Rtypes.h"
#include "extendVectorUtil.h"
#ifdef __CINT__
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Double32_t> > LV;
#pragma link C++ function ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Double32_t> >::operator+(LV);
#pragma link C++ function ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Double32_t> >::operator-(LV);
#pragma link C++ function ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Double32_t> >::Dot(LV);
#pragma link C++ function ROOT::Math::BoostZ::operator()(LV);
#pragma link C++ function ROOT::Math::Boost::operator()(LV);
#pragma link C++ function ROOT::Math::RotationZ::operator()(LV);
#pragma link C++ namespace ROOT::Math::VectorUtil+;
#endif
