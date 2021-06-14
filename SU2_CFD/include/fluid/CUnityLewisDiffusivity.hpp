/*!
 * \file CUnityLewisDiffusivity.hpp
 * \brief Defines unity Lewis mass diffusivity.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna, T. Economon
 * \version 7.0.6 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "CDiffusivityModel.hpp"

/*!
 * \class CUnityLewisDiffusivity
 * \brief Defines a unity lewis mass diffusivity model for species equations.
 * \author T. Economon
 */
class CUnityLewisDiffusivity final : public CDiffusivityModel {
public:
  /*!
   * \brief Constructor of the class.
   */
  // CUnityLewisDiffusivity() {}
  CUnityLewisDiffusivity(su2double kt_const) : kt_(kt_const) {}
  // CUnityLewisDiffusivity(su2double diff_const) :  diff_(diff_const){}
  
  su2double GetDiffusivity() const override {return diff_;}

  /*!
   * \brief Set diffusivity.
   */
  void SetDiffusivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) override { 
    diff_ = kt_ / (Lewis * rho * cp); 
  }

  private:
    su2double diff_{0.0}; 
    su2double kt_{0.0}; 
    su2double Lewis = 1;  
};
