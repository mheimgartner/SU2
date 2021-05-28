#include <cmath>
#include <math.h>
#include <vector>
#include <numeric>

#include "../../include/fluid/CFluidScalar.hpp"
#include "../../include/fluid/CConstantViscosity.hpp"
#include "../../include/fluid/CSutherland.hpp"
#include "../../include/fluid/CPolynomialViscosity.hpp"

#include "../../include/fluid/CIncIdealGas.hpp"

// CFluidScalar::CFluidScalar(CConfig *config, CFluidModel *FluidModel, su2double value_pressure_operating) : CFluidModel(), FluidModel(FluidModel){
CFluidScalar::CFluidScalar(CConfig *config, su2double value_pressure_operating) : CFluidModel() {

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  n_scalars = config->GetNScalars(); 
  // LaminarViscosityPointers[n_scalars]; 
  molarMasses.resize(n_scalars); 
  moleFractions.resize(n_scalars);
  laminarViscosity.resize(n_scalars);

  for (int iVar = 0; iVar < n_scalars; iVar++) {
    molarMasses.at(iVar) = config->GetMolecular_Weight(iVar); 
  }

  wilke = false;
  davidson = true;

  Pressure = value_pressure_operating;
  fluidModel = new CIncIdealGas(config->GetSpecific_Heat_Cp(), config->GetGas_Constant(), Pressure);

  // SetLaminarViscosityModel(config);
}

std::vector<su2double> CFluidScalar::massToMoleFractions(su2double* val_scalars){
  su2double mixtureMolarMass = 0.0; 
  su2double mixtureFractionDenumerator = 0.0; 
  // moleFractions.clear(); // is vector moleFractions automatically cleared after one SetTDState_T cycle? --> not needed because cleared in destructor of vector

  for(int iVar = 0; iVar < n_scalars; iVar++){
    mixtureMolarMass += val_scalars[iVar] / molarMasses[iVar]; 
  }

  for(int iVar = 0; iVar < n_scalars; iVar++){
    moleFractions.at(iVar) = (val_scalars[iVar] / molarMasses[iVar]) / mixtureMolarMass;
  }
  return moleFractions; 
}

void CFluidScalar::SetLaminarViscosityModel(const CConfig* config) {
  switch (config->GetKind_ViscosityModel()) {
    case VISCOSITYMODEL::CONSTANT:
    // Build a list of LaminarViscosity pointers to be used in wilkeViscosity to get the species viscosities. 
      for (int iVar = 0; iVar < n_scalars; iVar++){
        LaminarViscosityPointers[iVar] = unique_ptr<CConstantViscosity>(new CConstantViscosity(config->GetMu_Constant(iVar)));
      }
      break;
    case VISCOSITYMODEL::SUTHERLAND:
      for (int iVar = 0; iVar < n_scalars; iVar++){
        LaminarViscosityPointers[iVar] = unique_ptr<CSutherland>(new CSutherland(config->GetMu_Ref(iVar), config->GetMu_Temperature_Ref(iVar), config->GetMu_S(iVar)));
      }
      break;
    case VISCOSITYMODEL::POLYNOMIAL:
      for (int iVar = 0; iVar < n_scalars; iVar++){
        LaminarViscosityPointers[iVar] = unique_ptr<CPolynomialViscosity<N_POLY_COEFFS>>(
          new CPolynomialViscosity<N_POLY_COEFFS>(config->GetMu_PolyCoeffND()));
      }
      break;
    case VISCOSITYMODEL::FLAMELET:
      /* do nothing. Viscosity is obtained from the table and set in setTDState_T */
      break;
    default:
      SU2_MPI::Error("Viscosity model not available.", CURRENT_FUNCTION);
      break;
  }
}

su2double CFluidScalar::wilkeViscosity(su2double* val_scalars){
  std::vector<su2double> phi; 
  std::vector<su2double> wilkeNumerator;
  std::vector<su2double> wilkeDenumeratorSum;
  su2double wilkeDenumerator = 0.0; 
  wilkeDenumeratorSum.clear();
  wilkeNumerator.clear(); 
   
  massToMoleFractions(val_scalars); 

  Density = GetDensity();
  Temperature = GetTemperature();
/*  
  for (int iVar = 0; iVar < n_scalars; iVar++){
    LaminarViscosityPointers[iVar]->SetViscosity(Temperature, Density);
    laminarViscosity.at(iVar) = LaminarViscosityPointers[iVar]->GetViscosity();
  }
*/  
  for(int i = 0; i < n_scalars; i++){
    for(int j = 0; j < n_scalars; j++){
      phi.push_back(pow((1 + pow(laminarViscosity[i] / laminarViscosity[j],0.5) * pow(molarMasses[j] / molarMasses[i],0.25)),2) / pow(8 * (1 + molarMasses[i] / molarMasses[j]),0.5)); 
      wilkeDenumerator += moleFractions[j] * phi[j];   
    }
    wilkeDenumeratorSum.push_back(wilkeDenumerator); 
    wilkeDenumerator = 0.0; 
    phi.clear();
    wilkeNumerator.push_back(moleFractions[i] * laminarViscosity[i]); 
    viscosityMixture += wilkeNumerator[i] / wilkeDenumeratorSum[i];  
  }
  return viscosityMixture; 
}

su2double CFluidScalar::davidsonViscosity(su2double* val_scalars){
  su2double fluidity = 0.0;
  su2double E = 0.0; 
  su2double mixtureFractionDenumerator = 0.0;
  const su2double A = 0.375;   
  std::vector<su2double> mixtureFractions; 
  mixtureFractions.clear(); 
    
  massToMoleFractions(val_scalars);
  Density = GetDensity();
  Temperature = GetTemperature();

/*  
  for (int iVar = 0; iVar < n_scalars; iVar++){
    LaminarViscosityPointers[iVar]->SetViscosity(Temperature, Density);
    laminarViscosity.at(iVar) = LaminarViscosityPointers[iVar]->GetViscosity();
  }
*/
  for(int i = 0; i < n_scalars; i++){
    mixtureFractionDenumerator += moleFractions[i] * sqrt(molarMasses[i]);
  }

  for(int j = 0; j < n_scalars; j++){
    mixtureFractions.push_back((moleFractions[j] * sqrt(molarMasses[j])) / mixtureFractionDenumerator); 
  }
    
  for(int i = 0; i < n_scalars; i++){
    for(int j = 0; j < n_scalars; j++){
      E = (2*sqrt(molarMasses[i]) * sqrt(molarMasses[j])) / (molarMasses[i] + molarMasses[j]);
      fluidity += ((mixtureFractions[i] * mixtureFractions[j]) / (sqrt(laminarViscosity[i]) * sqrt(laminarViscosity[j]))) * pow(E, A);  
    }
  }
  return viscosityMixture = 1 / fluidity;
}

unsigned long CFluidScalar::SetTDState_T(su2double val_temperature, su2double *val_scalars){
  fluidModel->SetTDState_T(val_temperature, val_scalars);  

  if(wilke){
    Mu  = wilkeViscosity(val_scalars);
  }
  else if(davidson){
    Mu = davidsonViscosity(val_scalars); 
  }

  delete [] val_scalars;
  return 0;
}

