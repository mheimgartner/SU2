#pragma once

#include <vector>
#include <memory>

#include "CFluidModel.hpp"

class CFluidScalar : public CFluidModel {

protected:
  unsigned short n_scalars;
  su2double viscosityMixture = 0.0;
  // su2double Density;
  // su2double Temperature;
  su2double Gas_Constant = 0.0; 
  su2double Gamma = 0.0; 
  
  bool wilke;
  bool davidson; 

  std::vector<su2double> massFractions; 
  std::vector<su2double> moleFractions;  
  std::vector<su2double> molarMasses;
  std::vector<su2double> laminarViscosity = {1.118e-05, 1.795e-05}; 

  unique_ptr<CViscosityModel> LaminarViscosityPointers[10];  
  // using LaminarViscosityPointers = std::unique_ptr<std::unique_ptr<CViscosityModel>[100]>; 
  // unique_ptr<CViscosityModel[]> LaminarViscosityPointers(new CViscosityModel[100]); 
  // auto LaminarViscosityPointers = std::make_unique<CViscosityModel[]>(100);

  CFluidModel *fluidModel = nullptr;

 public:
  // CFluidScalar(CConfig *config, CFluidModel *fluidModel);
  CFluidScalar(CConfig *config, su2double value_pressure_operating);

  ~CFluidScalar() {delete fluidModel;};

  void SetLaminarViscosityModel(const CConfig* config);

  unsigned long SetTDState_T(su2double val_temperature, su2double *val_scalars);

  inline su2double GetTemperature() {return fluidModel->GetTemperature(); }

  inline su2double GetDensity() {return fluidModel->GetDensity(); }

  // std::vector<su2double> massToMoleFractions(su2double* val_scalars);
  std::vector<su2double> massToMoleFractions(std::vector<su2double> massFractions);

  su2double wilkeViscosity(su2double* val_scalars);

  // su2double davidsonViscosity(su2double* val_scalars);
  su2double davidsonViscosity(std::vector<su2double> massFractions);

  inline su2double GetLaminarViscosity() {return Mu; }

  inline unsigned short GetNScalars() { return n_scalars; }
};