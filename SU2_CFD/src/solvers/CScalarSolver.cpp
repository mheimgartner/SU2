/*!
 * \file CScalarSolver.cpp
 * \brief Main subroutines for the  transported scalar model.
 * \author D. Mayer, T. Economon
 * \version 7.1.1 "Blackbird"
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

#include "../../include/solvers/CScalarSolver.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"


CScalarSolver::CScalarSolver(void) : CSolver() {

  lowerlimit       = nullptr;
  upperlimit       = nullptr;
  //Scalar_Inf       = NULL;  
}

CScalarSolver::CScalarSolver(CGeometry* geometry, CConfig *config) : CSolver() {

 
  lowerlimit       = nullptr;
  upperlimit       = nullptr;
  //Scalar_Inf       = NULL;
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  nMarker = config->GetnMarker_All();

  /*--- Store the number of vertices on each marker for deallocation later ---*/
  nVertex.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

#ifdef HAVE_OMP
  /*--- Get the edge coloring, see notes in CEulerSolver's constructor. ---*/
  su2double parallelEff = 1.0;
  const auto& coloring = geometry->GetEdgeColoring(&parallelEff);

  ReducerStrategy = parallelEff < COLORING_EFF_THRESH;

  if (ReducerStrategy && (coloring.getOuterSize()>1))
    geometry->SetNaturalEdgeColoring();

  if (!coloring.empty()) {
    auto groupSize = ReducerStrategy? 1ul : geometry->GetEdgeColorGroupSize();
    auto nColor = coloring.getOuterSize();
    EdgeColoring.reserve(nColor);

    for(auto iColor = 0ul; iColor < nColor; ++iColor)
      EdgeColoring.emplace_back(coloring.innerIdx(iColor), coloring.getNumNonZeros(iColor), groupSize);
  }

  nPoint = geometry->GetnPoint();
  omp_chunk_size = computeStaticChunkSize(nPoint, omp_get_max_threads(), OMP_MAX_SIZE);
#else
  EdgeColoring[0] = DummyGridColor<>(geometry->GetnEdge());
#endif
}

CScalarSolver::~CScalarSolver(void) {

  if (lowerlimit != nullptr)    delete [] lowerlimit;
  if (upperlimit != nullptr)    delete [] upperlimit;

  //if (Scalar_Inf != NULL)    delete [] Scalar_Inf;

  for (auto& mat : SlidingState) {
    for (auto ptr : mat) delete [] ptr;
  }

  delete nodes;
}

void CScalarSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                                  CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool muscl = config->GetMUSCL_Scalar();
  const bool limiter = (config->GetKind_SlopeLimit_Scalar() != NO_LIMITER);

  /*--- Only reconstruct flow variables if MUSCL is on for flow (requires upwind) and scalar. ---*/
  const bool musclFlow = config->GetMUSCL_Flow() && muscl &&
                        (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND);
  /*--- Only consider flow limiters for cell-based limiters, edge-based would need to be recomputed. ---*/
  const bool limiterFlow = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) &&
                           (config->GetKind_SlopeLimit_Flow() != VAN_ALBADA_EDGE);

  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Static arrays of MUSCL-reconstructed flow primitives and turbulence variables (thread safety). ---*/
  su2double solution_i[MAXNVAR] = {0.0}, flowPrimVar_i[MAXNVARFLOW] = {0.0};
  su2double solution_j[MAXNVAR] = {0.0}, flowPrimVar_j[MAXNVARFLOW] = {0.0};

  /*--- Loop over edge colors. ---*/
  for (auto color : EdgeColoring)
  {
  /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
  SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
  for(auto k = 0ul; k < color.size; ++k) {

    auto iEdge = color.indices[k];

    unsigned short iDim, iVar;

    /*--- Points in edge and normal vectors ---*/

    auto iPoint = geometry->edges->GetNode(iEdge,0);
    auto jPoint = geometry->edges->GetNode(iEdge,1);

    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    /*--- Primitive variables w/o reconstruction ---*/

    const auto V_i = flowNodes->GetPrimitive(iPoint);
    const auto V_j = flowNodes->GetPrimitive(jPoint);
    numerics->SetPrimitive(V_i, V_j);

    /*--- Scalar variables w/o reconstruction ---*/

    const auto Scalar_i = nodes->GetSolution(iPoint);
    const auto Scalar_j = nodes->GetSolution(jPoint);
    numerics->SetScalarVar(Scalar_i, Scalar_j);

    /*--- Grid Movement ---*/

    if (dynamic_grid)
      numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                           geometry->nodes->GetGridVel(jPoint));

    if (muscl || musclFlow) {
      const su2double *Limiter_i = nullptr, *Limiter_j = nullptr;

      const auto Coord_i = geometry->nodes->GetCoord(iPoint);
      const auto Coord_j = geometry->nodes->GetCoord(jPoint);

      su2double Vector_ij[MAXNDIM] = {0.0};
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_ij[iDim] = 0.5*(Coord_j[iDim] - Coord_i[iDim]);
      }

      if (musclFlow) {
        /*--- Reconstruct mean flow primitive variables. ---*/

        auto Gradient_i = flowNodes->GetGradient_Reconstruction(iPoint);
        auto Gradient_j = flowNodes->GetGradient_Reconstruction(jPoint);

        if (limiterFlow) {
          Limiter_i = flowNodes->GetLimiter_Primitive(iPoint);
          Limiter_j = flowNodes->GetLimiter_Primitive(jPoint);
        }

        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnPrimVarGrad(); iVar++) {
          su2double Project_Grad_i = 0.0, Project_Grad_j = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Project_Grad_i += Vector_ij[iDim]*Gradient_i[iVar][iDim];
            Project_Grad_j -= Vector_ij[iDim]*Gradient_j[iVar][iDim];
          }
          if (limiterFlow) {
            Project_Grad_i *= Limiter_i[iVar];
            Project_Grad_j *= Limiter_j[iVar];
          }
          flowPrimVar_i[iVar] = V_i[iVar] + Project_Grad_i;
          flowPrimVar_j[iVar] = V_j[iVar] + Project_Grad_j;
        }

        numerics->SetPrimitive(flowPrimVar_i, flowPrimVar_j);
      }

      if (muscl) {
        /*--- Reconstruct turbulence variables. ---*/

        auto Gradient_i = nodes->GetGradient_Reconstruction(iPoint);
        auto Gradient_j = nodes->GetGradient_Reconstruction(jPoint);

        if (limiter) {
          Limiter_i = nodes->GetLimiter(iPoint);
          Limiter_j = nodes->GetLimiter(jPoint);
        }

        for (iVar = 0; iVar < nVar; iVar++) {
          su2double Project_Grad_i = 0.0, Project_Grad_j = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Project_Grad_i += Vector_ij[iDim]*Gradient_i[iVar][iDim];
            Project_Grad_j -= Vector_ij[iDim]*Gradient_j[iVar][iDim];
          }
          if (limiter) {
            Project_Grad_i *= Limiter_i[iVar];
            Project_Grad_j *= Limiter_j[iVar];
          }
          solution_i[iVar] = Scalar_i[iVar] + Project_Grad_i;
          solution_j[iVar] = Scalar_j[iVar] + Project_Grad_j;
        }

        numerics->SetScalarVar(solution_i, solution_j);
      }
    }

    /*--- Update convective residual value ---*/

    auto residual = numerics->ComputeResidual(config);

    if (ReducerStrategy) {
      EdgeFluxes.SetBlock(iEdge, residual);
      if (implicit) Jacobian.SetBlocks(iEdge, residual.jacobian_i, residual.jacobian_j);
    }
    else {
      LinSysRes.AddBlock(iPoint, residual);
      LinSysRes.SubtractBlock(jPoint, residual);
      if (implicit) Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
    }

    /*--- Viscous contribution. ---*/

    Viscous_Residual(iEdge, geometry, solver_container,
                     numerics_container[VISC_TERM + omp_get_thread_num()*MAX_TERMS], config);
  }
  } // end color loop

  if (ReducerStrategy) {
    SumEdgeFluxes(geometry);
    if (implicit) Jacobian.SetDiagonalAsColumnSum();
  }
}

void CScalarSolver::Viscous_Residual(unsigned long iEdge, CGeometry *geometry, CSolver **solver_container,
                                   CNumerics *numerics, CConfig *config) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  /*--- Points in edge ---*/

  auto iPoint = geometry->edges->GetNode(iEdge,0);
  auto jPoint = geometry->edges->GetNode(iEdge,1);

  /*--- Points coordinates, and normal vector ---*/

  numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                     geometry->nodes->GetCoord(jPoint));
  numerics->SetNormal(geometry->edges->GetNormal(iEdge));

  /*--- Conservative variables w/o reconstruction ---*/

  numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint),
                         flowNodes->GetPrimitive(jPoint));
    
  /*--- Scalar variables w/o reconstruction, and its gradients ---*/  

  numerics->SetScalarVar(nodes->GetSolution(iPoint),
                         nodes->GetSolution(jPoint));
  numerics->SetScalarVarGradient(nodes->GetGradient(iPoint),
                                 nodes->GetGradient(jPoint));

  /*--- Mass diffusivity coefficients. ---*/

  numerics->SetDiffusionCoeff(nodes->GetDiffusivity(iPoint),
                              nodes->GetDiffusivity(jPoint));

  /*--- Compute residual, and Jacobians ---*/

  auto residual = numerics->ComputeResidual(config);

  if (ReducerStrategy) {
    EdgeFluxes.SubtractBlock(iEdge, residual);
    if (implicit) Jacobian.UpdateBlocksSub(iEdge, residual.jacobian_i, residual.jacobian_j);
  }
  else {
    LinSysRes.SubtractBlock(iPoint, residual);
    LinSysRes.AddBlock(jPoint, residual);
    if (implicit) Jacobian.UpdateBlocksSub(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
  }
}

void CScalarSolver::SumEdgeFluxes(CGeometry* geometry) {

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {

    LinSysRes.SetBlock_Zero(iPoint);

    for (auto iEdge : geometry->nodes->GetEdges(iPoint)) {
      if (iPoint == geometry->edges->GetNode(iEdge,0))
        LinSysRes.AddBlock(iPoint, EdgeFluxes.GetBlock(iEdge));
      else
        LinSysRes.SubtractBlock(iPoint, EdgeFluxes.GetBlock(iEdge));
    }
  }

}

void CScalarSolver::BC_Sym_Plane(CGeometry      *geometry,
                               CSolver        **solver_container,
                               CNumerics      *conv_numerics,
                               CNumerics      *visc_numerics,
                               CConfig        *config,
                               unsigned short val_marker) {

  /*--- Convective and viscous fluxes across symmetry plane are equal to zero. ---*/

}

void CScalarSolver::BC_Euler_Wall(CGeometry      *geometry,
                                CSolver        **solver_container,
                                CNumerics      *conv_numerics,
                                CNumerics      *visc_numerics,
                                CConfig        *config,
                                unsigned short val_marker) {

  /*--- Convective fluxes across euler wall are equal to zero. ---*/

}

void CScalarSolver::BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                                  CNumerics *numerics, CConfig *config) {

  /*--- Complete residuals for periodic boundary conditions. We loop over
   the periodic BCs in matching pairs so that, in the event that there are
   adjacent periodic markers, the repeated points will have their residuals
   accumulated corectly during the communications. For implicit calculations
   the Jacobians and linear system are also correctly adjusted here. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
  }

}

void CScalarSolver::PrepareImplicitIteration(CGeometry *geometry, CSolver** solver_container, CConfig *config) {

 const auto flowNodes = solver_container[FLOW_SOL]->GetNodes();


/*--- use preconditioner for scalar transport equations ---*/
 const bool incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
// if (incompressible) SetPreconditioner(geometry, solver_container, config);


  /*--- Set shared residual variables to 0 and declare
   *    local ones for current thread to work on. ---*/

  SetResToZero();

  su2double resMax[MAXNVAR] = {0.0}, resRMS[MAXNVAR] = {0.0};
  const su2double* coordMax[MAXNVAR] = {nullptr};
  unsigned long idxMax[MAXNVAR] = {0};

  /*--- Build implicit system ---*/

  SU2_OMP(for schedule(static,omp_chunk_size) nowait)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /// TODO: This could be the SetTime_Step of this solver.
     // nijso: we use a global cfl scaling factor
    //Delta = Vol / (config->GetCFLRedCoeff_Scalar()*solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint));
    //su2double dt = nodes->GetLocalCFL(iPoint) / flowNodes->GetLocalCFL(iPoint) * flowNodes->GetDelta_Time(iPoint);
    su2double dt = config->GetCFLRedCoeff_Scalar() * flowNodes->GetDelta_Time(iPoint);
    nodes->SetDelta_Time(iPoint, dt);

    /*--- Modify matrix diagonal to improve diagonal dominance. ---*/

    if (dt != 0.0) {
      su2double Vol = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
      Jacobian.AddVal2Diag(iPoint, Vol / dt);
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      LinSysRes.SetBlock_Zero(iPoint);
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      unsigned long total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = -LinSysRes[total_index];
      LinSysSol[total_index] = 0.0;

      su2double Res = fabs(LinSysRes[total_index]);
      resRMS[iVar] += Res*Res;
      if (Res > resMax[iVar]) {
        resMax[iVar] = Res;
        idxMax[iVar] = iPoint;
        coordMax[iVar] = geometry->nodes->GetCoord(iPoint);
      }
    }
  }

  // SetPreconditioner();


  SU2_OMP_CRITICAL
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Residual_RMS[iVar] += resRMS[iVar];
    AddRes_Max(iVar, resMax[iVar], geometry->nodes->GetGlobalIndex(idxMax[iVar]), coordMax[iVar]);
  }
  SU2_OMP_BARRIER

  /*--- Compute the root mean square residual ---*/
  SU2_OMP_MASTER
  SetResidual_RMS(geometry, config);
  SU2_OMP_BARRIER
}

void CScalarSolver::CompleteImplicitIteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  const bool compressible = (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  const auto flowNodes = solver_container[FLOW_SOL]->GetNodes();


  // nijso: TODO: we also still have an underrelaxation factor as config option
  ComputeUnderRelaxationFactor(config);

  /*--- Update solution (system written in terms of increments) ---*/

  if (!adjoint) {

  /*--- Update the scalar solution. ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

      su2double density = flowNodes->GetDensity(iPoint);
      su2double density_old = density;

      if (compressible)
        density_old = flowNodes->GetSolution_Old(iPoint,0);

      // nijso: check difference between conservative and regular solution
      // nijso: conservative solution is for transport of rho*Y, 
      // nijso: check underrelaxation
      for (unsigned short iVar = 0u; iVar < nVar; iVar++) {
        //nodes->AddSolution(iPoint,iVar, nodes->GetUnderRelaxation(iPoint)*LinSysSol[iPoint*nVar+iVar]);
        // FIXME dan: Underrelaxation is turned off here beacuse problems.
          //nodes->GetUnderRelaxation(iPoint)*LinSysSol(iPoint,iVar),
        nodes->AddClippedSolution(iPoint, iVar, LinSysSol(iPoint,iVar), lowerlimit[iVar], upperlimit[iVar]);
      }
    }                
  }


  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  InitiateComms(geometry, config, SOLUTION_EDDY);
  CompleteComms(geometry, config, SOLUTION_EDDY);

}

void CScalarSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  PrepareImplicitIteration(geometry, solver_container, config);

  /*--- Solve or smooth the linear system. ---*/

  SU2_OMP(for schedule(static,OMP_MIN_SIZE) nowait)
  for (unsigned long iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    LinSysRes.SetBlock_Zero(iPoint);
    LinSysSol.SetBlock_Zero(iPoint);
  }

  auto iter = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  SU2_OMP_MASTER {
    SetIterLinSolver(iter);
    SetResLinSolver(System.GetResidual());
  }
  SU2_OMP_BARRIER

  CompleteImplicitIteration(geometry, solver_container, config);

}

void CScalarSolver::ComputeUnderRelaxationFactor(const CConfig *config) {


  /* Loop over the solution update given by relaxing the linear
   system for this nonlinear iteration. */

  su2double localUnderRelaxation =  1.00;
  
  /*--- We could use the relaxationfactor from the config file and use it, as an allowableratio  ---*/
  //const su2double allowableRatio =  1.0 - config->GetRelaxation_Factor_Scalar();

  const su2double allowableRatio = 0.90;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    localUnderRelaxation = 1.0;
    //if (sa_model) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {

        /* We impose a limit on the maximum percentage that the
         scalar variables can change over a nonlinear iteration. */

        const unsigned long index = iPoint * nVar + iVar;
        su2double ratio = fabs(LinSysSol[index]) / (fabs(nodes->GetSolution(iPoint, iVar)) + EPS);
        if (ratio > allowableRatio) {
          localUnderRelaxation = min(allowableRatio / ratio, localUnderRelaxation);
        }

      }
    //}

    /* Threshold the relaxation factor in the event that there is
     a very small value. This helps avoid catastrophic crashes due
     to non-realizable states by canceling the update. */

    if (localUnderRelaxation < 1e-10) localUnderRelaxation = 0.0;

    /* Store the under-relaxation factor for this point. */

    nodes->SetUnderRelaxation(iPoint, localUnderRelaxation);

  }

}

void CScalarSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                       unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {

  const bool sst_model = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool first_order = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST);
  const bool second_order = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);

  /*--- Flow solution, needed to get density. ---*/

  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  /*--- Store the physical time step ---*/

  const su2double TimeStep = config->GetDelta_UnstTimeND();

  /*--- Local variables ---*/

  unsigned short iVar, iMarker, iDim, iNeigh;
  unsigned long iPoint, jPoint, iVertex, iEdge;

  const su2double *U_time_nM1 = nullptr, *U_time_n = nullptr, *U_time_nP1 = nullptr;
  su2double Volume_nM1, Volume_nP1;
  su2double Density_nM1, Density_n, Density_nP1;
  const su2double *Normal = nullptr, *GridVel_i = nullptr, *GridVel_j = nullptr;
  su2double Residual_GCL;

  /*--- Compute the dual time-stepping source term for static meshes ---*/

  if (!dynamic_grid) {

    /*--- Loop over all nodes (excluding halos) ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/

      U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      U_time_n   = nodes->GetSolution_time_n(iPoint);
      U_time_nP1 = nodes->GetSolution(iPoint);

      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/

      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/
 
       /*--- If this is the SST model, we need to multiply by the density
         in order to get the conservative variables. 
         This is also the case in general for our transported scalars: 
         we transport rho*Y                                            ---*/

        if (incompressible){
          /*--- This is temporary and only valid for constant-density problems:
          density could also be temperature dependent, but as it is not a part
          of the solution vector it's neither stored for previous time steps
          nor updated with the solution at the end of each iteration. */
          Density_nM1 = flowNodes->GetDensity(iPoint);
          Density_n   = flowNodes->GetDensity(iPoint);
          Density_nP1 = flowNodes->GetDensity(iPoint);
        }
        else{
          Density_nM1 = flowNodes->GetSolution_time_n1(iPoint)[0];
          Density_n   = flowNodes->GetSolution_time_n(iPoint,0);
          Density_nP1 = flowNodes->GetSolution(iPoint,0);
        }

        for (iVar = 0; iVar < nVar; iVar++) {
          if (first_order)
            LinSysRes(iPoint,iVar) += ( Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*Volume_nP1 / TimeStep;
          if (second_order)
            LinSysRes(iPoint,iVar) += ( 3.0*Density_nP1*U_time_nP1[iVar] - 4.0*Density_n*U_time_n[iVar]
                                       +1.0*Density_nM1*U_time_nM1[iVar] ) * Volume_nP1/(2.0*TimeStep);
        }


      /*--- Compute the Jacobian contribution due to the dual time source term. ---*/
      if (implicit) {
        if (first_order) Jacobian.AddVal2Diag(iPoint, Volume_nP1/TimeStep);
        if (second_order) Jacobian.AddVal2Diag(iPoint, (Volume_nP1*3.0)/(2.0*TimeStep));
      }

    }

  } else {

    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
     dynamically deforming), the Geometric Conservation Law (GCL) should be
     satisfied in conjunction with the ALE formulation of the governing
     equations. The GCL prevents accuracy issues caused by grid motion, i.e.
     a uniform free-stream should be preserved through a moving grid. First,
     we will loop over the edges and boundaries to compute the GCL component
     of the dual time source term that depends on grid velocities. ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

      GridVel_i = geometry->nodes->GetGridVel(iPoint);
      U_time_n  = nodes->GetSolution_time_n(iPoint);
      Density_n = 1.0;

      //if (sst_model) {
        if (incompressible)
          Density_n = flowNodes->GetDensity(iPoint); // Temporary fix
      //  else
      //    Density_n = flowNodes->GetSolution_time_n(iPoint,0);
      //}

      for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {

        iEdge = geometry->nodes->GetEdge(iPoint, iNeigh);
        Normal = geometry->edges->GetNormal(iEdge);

        jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
        GridVel_j = geometry->nodes->GetGridVel(jPoint);

        /*--- Determine whether to consider the normal outward or inward. ---*/
        su2double dir = (iPoint < jPoint)? 0.5 : -0.5;

        Residual_GCL = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual_GCL += dir*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

        Residual_GCL *= Density_n;

        for (iVar = 0; iVar < nVar; iVar++)
          LinSysRes(iPoint,iVar) += U_time_n[iVar]*Residual_GCL;
      }
    }

    /*--- Loop over the boundary edges ---*/

    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {

        SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

          /*--- Get the index for node i plus the boundary face normal ---*/

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          /*--- Grid velocities stored at boundary node i ---*/

          GridVel_i = geometry->nodes->GetGridVel(iPoint);

          /*--- Compute the GCL term by dotting the grid velocity with the face
           normal. The normal is negated to match the boundary convention. ---*/

          Residual_GCL = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Residual_GCL -= 0.5*(GridVel_i[iDim]+GridVel_i[iDim])*Normal[iDim];

          /*--- Compute the GCL component of the source term for node i ---*/

          U_time_n = nodes->GetSolution_time_n(iPoint);

          /*--- Multiply by density at node i for the SST model ---*/

          //if (sst_model) {
            if (incompressible)
              Density_n = flowNodes->GetDensity(iPoint); // Temporary fix
            else
              Density_n = flowNodes->GetSolution_time_n(iPoint,0);

            for (iVar = 0; iVar < nVar; iVar++)
              LinSysRes(iPoint,iVar) += Density_n*U_time_n[iVar]*Residual_GCL;
          //}
          //else {
          //  for (iVar = 0; iVar < nVar; iVar++)
          //    LinSysRes(iPoint,iVar) += U_time_n[iVar]*Residual_GCL;
          //}

        }
      }
    }

    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/

      U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      U_time_n   = nodes->GetSolution_time_n(iPoint);
      U_time_nP1 = nodes->GetSolution(iPoint);

      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/

      Volume_nM1 = geometry->nodes->GetVolume_nM1(iPoint);
      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/

      //if (sst_model) {

        /*--- If this is the SST model, we need to multiply by the density
         in order to get the conservative variables ---*/
        if (incompressible) {
          /*--- This is temporary and only valid for constant-density problems:
          density could also be temperature dependent, but as it is not a part
          of the solution vector it's neither stored for previous time steps
          nor updated with the solution at the end of each iteration. */
          Density_nM1 = flowNodes->GetDensity(iPoint);
          Density_n   = flowNodes->GetDensity(iPoint);
          Density_nP1 = flowNodes->GetDensity(iPoint);
        }
        else {
          Density_nM1 = flowNodes->GetSolution_time_n1(iPoint)[0];
          Density_n   = flowNodes->GetSolution_time_n(iPoint,0);
          Density_nP1 = flowNodes->GetSolution(iPoint,0);
        }

        for (iVar = 0; iVar < nVar; iVar++) {
          if (first_order)
            LinSysRes(iPoint,iVar) += (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(Volume_nP1/TimeStep);
          if (second_order)
            LinSysRes(iPoint,iVar) += (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
                                      + (Density_nM1*U_time_nM1[iVar] - Density_n*U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
        }

      //} else {
      //
      //  for (iVar = 0; iVar < nVar; iVar++) {
      //    if (first_order)
      //      LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
      //    if (second_order)
      //      LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
      //                                 + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      //  }
      //}

      /*--- Compute the Jacobian contribution due to the dual time source term. ---*/
      if (implicit) {
        if (first_order) Jacobian.AddVal2Diag(iPoint, Volume_nP1/TimeStep);
        if (second_order) Jacobian.AddVal2Diag(iPoint, (Volume_nP1*3.0)/(2.0*TimeStep));
      }
    }

  } // end dynamic grid

}


void CScalarSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/

  unsigned short iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;
  su2double Area_Children, Area_Parent;
  const su2double *Solution_Fine = nullptr;
  //bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
  //                  (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));
  //bool time_stepping = (config->GetTime_Marching() == TIME_STEPPING);
  //unsigned short iZone = config->GetiZone();
  //unsigned short nZone = config->GetnZone();
  
  //string UnstExt, text_line;
  //ifstream restart_file;

  string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  bool turbulent = ((config->GetKind_Solver() == RANS) ||
                    (config->GetKind_Solver() == INC_RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_INC_RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_RANS));
  
  unsigned short turbSkip = 0;
  if (turbulent) turbSkip = solver[MESH_0][TURB_SOL]->GetnVar();

  /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
  SU2_OMP_MASTER
  {

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Skip flow variables ---*/

  unsigned short skipVars = 0;

  if (nDim == 2) skipVars += 6;
  if (nDim == 3) skipVars += 8;
  
  /*--- Skip turbulent variables if necessary variables ---*/

  if (turbulent) skipVars += turbSkip;

  /*--- Adjust the number of solution variables in the incompressible
   restart. We always carry a space in nVar for the energy equation in the
   mean flow solver, but we only write it to the restart if it is active.
   Therefore, we must reduce skipVars here if energy is inactive so that
   the turbulent variables are read correctly. ---*/

  bool incompressible       = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
  bool energy               = config->GetEnergy_Equation();
  bool weakly_coupled_heat  = config->GetWeakly_Coupled_Heat();

  if (incompressible && ((!energy) && (!weakly_coupled_heat))) skipVars--;

  /*--- Load data from the restart into correct containers. ---*/

  unsigned long counter = 0, iPoint_Global = 0;
  for (; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {


    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; ++iVar)
        nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index+iVar]);

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (counter != nPointDomain) {
    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  } // end SU2_OMP_MASTER, pre and postprocessing are thread-safe.
  SU2_OMP_BARRIER

  /*--- MPI solution and compute the scalars ---*/

  solver[MESH_0][SCALAR_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][SCALAR_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);

  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  solver[MESH_0][SCALAR_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->nodes->GetVolume(iPoint);
      su2double Solution_Coarse[MAXNVAR] = {0.0};
      for (iChildren = 0; iChildren < geometry[iMesh]->nodes->GetnChildren_CV(iPoint); iChildren++) {
        Point_Fine = geometry[iMesh]->nodes->GetChildren_CV(iPoint, iChildren);
        Area_Children = geometry[iMesh-1]->nodes->GetVolume(Point_Fine);
        Solution_Fine = solver[iMesh-1][SCALAR_SOL]->GetNodes()->GetSolution(Point_Fine);
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution_Coarse[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][SCALAR_SOL]->GetNodes()->SetSolution(iPoint,Solution_Coarse);
    }
    solver[iMesh][SCALAR_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][SCALAR_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    solver[iMesh][SCALAR_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
  }

  /*--- Go back to single threaded execution. ---*/
  SU2_OMP_MASTER
  {
  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars; Restart_Vars = nullptr;
  delete [] Restart_Data; Restart_Data = nullptr;

  } // end SU2_OMP_MASTER
  SU2_OMP_BARRIER
  
}

void CScalarSolver::BC_HeatFlux_Wall(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config,
                                     unsigned short val_marker) {
  /*--- Convective fluxes across viscous walls are equal to zero. ---*/
  
}

void CScalarSolver::BC_Far_Field(CGeometry *geometry,
                                 CSolver **solver_container,
                                 CNumerics *conv_numerics,
                                 CNumerics *visc_numerics, CConfig *config,
                                 unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *Normal, *V_infty, *V_domain;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->nodes->GetDomain(iPoint)) {
      
      /*--- Allocate the value at the infinity ---*/
      
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Set scalar variable at the wall, and at infinity ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);
        Solution_j[iVar] = Scalar_Inf[iVar];
      }
      conv_numerics->SetScalarVar(Solution_i, Solution_j);
      
      /*--- Set Normal (it is necessary to change the sign) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute residuals and Jacobians ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  delete [] Normal;
  
}



void CScalarSolver::SetInletAtVertex(const su2double *val_inlet,
                                     unsigned short iMarker,
                                     unsigned long iVertex) {
  unsigned short iVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Inlet_ScalarVars[iMarker][iVertex][iVar] = val_inlet[Inlet_Position+iVar];
  }
  
}

su2double CScalarSolver::GetInletAtVertex(su2double       *val_inlet,
                                          unsigned long    val_inlet_point,
                                          unsigned short   val_kind_marker,
                                          string           val_marker,
                                          const CGeometry *geometry,
                                          const CConfig   *config) const {
  
  /*--- Local variables ---*/
  
  unsigned short iMarker, iDim, iVar;
  unsigned long iPoint, iVertex;
  su2double Area = 0.0;
  su2double Normal[3] = {0.0,0.0,0.0};
  
  /*--- Alias positions within inlet file for readability ---*/
  
  if (val_kind_marker == INLET_FLOW) {
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) &&
          (config->GetMarker_All_TagBound(iMarker) == val_marker)) {
        
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++){
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (iPoint == val_inlet_point) {
            
            /*-- Compute boundary face area for this vertex. ---*/
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            Area = sqrt(Area);
            
            /*--- Access and store the inlet variables for this vertex. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++)
              val_inlet[Inlet_Position+iVar] = Inlet_ScalarVars[iMarker][iVertex][iVar];
            
            /*--- Exit once we find the point. ---*/
            
            return Area;
            
          }
        }
      }
    }
    
  }
  
  /*--- If we don't find a match, then the child point is not on the
   current inlet boundary marker. Return zero area so this point does
   not contribute to the restriction operator and continue. ---*/
  
  return Area;
  
}

void CScalarSolver::SetUniformInlet(const CConfig* config, unsigned short iMarker) {
  
  for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Inlet_ScalarVars[iMarker][iVertex][iVar] = Scalar_Inf[iVar];
  }

}


