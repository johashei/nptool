/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: V. Alcindor mail: valcindor@ikp.tu-darmstadt.de          *
 *                                                                           *
 * Creation Date  : July 2022                                          *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This files is for the generation of 15F Breit Wigner shape in the case   *
 *  of the sequential 2p emission 14O+p -> 15F -> 14O(1-) + p -> 13N + p     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_coulomb.h>

// dirty way of avoiding a header or input file...
namespace input_conditions {
// Incoming beam, here 14O
NPL::Particle Beam(8, 14);
// Target particle, here protons
NPL::Particle Target(1, 1);
// Intermediate state in the sequential emission (14O(1-), Ex=5.173 MeV)
NPL::Particle Intermediate(8, 14);
double        IntermediateState = 5.173;
// Final particle in the exit channel, here 13N
NPL::Particle Final(7, 13);
} // namespace input_conditions

using namespace input_conditions;

double CoulombFunction(std::string type, double E, double R, double mu,
                       Int_t Z1, Int_t Z2, double L, Int_t k = 0) {
  gsl_sf_result F_res;
  gsl_sf_result Fp_res;
  gsl_sf_result G_res;
  gsl_sf_result Gp_res;
  double        expF;
  double        expG;
  double        eta, rho;
  eta = 0.15748 * Z1 * Z2 * sqrt(mu / E);
  rho = 0.218735 * sqrt(mu * E) * R;
  gsl_sf_coulomb_wave_FG_e(eta, rho, L, k, &F_res, &Fp_res, &G_res, &Gp_res,
                           &expF, &expG);
  if (type == "F")
    return F_res.val * exp(expF);
  else if (type == "G")
    return G_res.val * exp(expG);
  else
    return 1;
}

double rho(double* E, double* par) {
  double R      = par[0];
  double mu_uma = par[1];

  return R * .218735 * sqrt(mu_uma * E[0]);
} // computation of the reduced variable for the penetrability

double P_l(double* E, double* par) {

  gsl_sf_result F_res;
  gsl_sf_result Fp_res;
  gsl_sf_result G_res;
  gsl_sf_result Gp_res;

  double R             = par[0];
  double mu_uma        = par[1];
  double Z1            = par[2];
  double Z2            = par[3];
  double L             = par[4];
  double par_P[2]      = {par[0], par[1]};
  double coulF         = CoulombFunction("F", E[0], R, mu_uma, Z1, Z2, L, 0);
  double coulG         = CoulombFunction("G", E[0], R, mu_uma, Z1, Z2, L, 0);
  double Penetrability = rho(E, par_P) / (pow(coulF, 2.) + pow(coulG, 2.));

  return Penetrability;
} // computation of the penetrability cf Abrahamovitz & Stegun Ch14 et C.

double f_theoretical_proton_width(double* x, double* par) {
  double Er = par[0];
  double l  = par[1];
  double z  = par[2];
  double a  = par[3];
  double Z  = par[4];
  double A  = par[5];
  double E  = x[0];

  double r0 = 1.2; // fm
  double R  = r0 * (pow(A, 1. / 3.) + pow(a, 1. / 3.)); // fm

  double mBeam_uma = Beam.GetA();
  +Beam.GetMassExcess() / 1000.;
  double mTarget_uma = Target.GetA();
  +Target.GetMassExcess() / 1000.;
  double mCN_uma = (mBeam_uma * mTarget_uma) / (mBeam_uma + mTarget_uma);

  double P_par[5] = {R, mCN_uma, z, Z, l};
  double P        = P_l(&E, P_par);

  double h_bar_c = 197.3; // MeV.fm
  double max_proton_width
      = 3. * pow(h_bar_c, 2.) / (((A * a) / (A + a) * 931.5) * pow(R, 2.));

  double proton_width = max_proton_width * P;

  if (E > 3.339 && E < 3.99) {
    double x1 = 3.339;
    double x2 = 3.99;
    P         = P_l(&x1, P_par);
    double y1 = max_proton_width * P;
    P         = P_l(&x2, P_par);
    double y2 = max_proton_width * P;

    return (y1 - y2) / (x1 - x2) * E + (y2 - x2 * (y1 - y2) / (x1 - x2));
  } // trying to correct GSLCoulomb problem in this region

  return proton_width;
} // Computation of proton width using COULOMB function from GSL

double f_proton_width(double* x, double* par) {
  double Er        = par[0];
  double res_width = par[1];
  double l         = par[2];
  double z         = par[3];
  double a         = par[4];
  double Z         = par[5];
  double A         = par[6];
  double E         = x[0];

  double r0 = 1.2; // fm
  double R  = r0 * (pow(A, 1. / 3.) + pow(a, 1. / 3.)); // fm

  double mBeam_uma   = Beam.GetA() + Beam.GetMassExcess() / 1000.;
  double mTarget_uma = Target.GetA() + Target.GetMassExcess() / 1000.;
  double mCN_uma     = (mBeam_uma * mTarget_uma) / (mBeam_uma + mTarget_uma);

  double P_par[5] = {R, mCN_uma, z, Z, l};
  double P        = P_l(&E, P_par);
  double Pr       = P_l(&Er, P_par);

  double proton_width = (res_width * P / Pr);

  if (E > 3.339 && E < 3.99) {
    double x1 = 3.339;
    double x2 = 3.99;
    P         = P_l(&x1, P_par);
    double y1 = (res_width * P / Pr);
    P         = P_l(&x2, P_par);
    double y2 = (res_width * P / Pr);

    return (y1 - y2) / (x1 - x2) * E + (y2 - x2 * (y1 - y2) / (x1 - x2));
  } // trying to correct GSLCoulomb problem in this region for 15F

  return proton_width;
} // Computation of proton width using COULOMB function from GSL

double f_double_sequential_proton_cross_section(double* x, double* par) {
  double E = x[0]; // MeV
  if (E < 0.) {
    std::cerr << "ENERGY IS NEGATIVE" << std::endl;
    return 1.;
  } else {

    double Er        = par[0]; // MeV
    double res_width = par[1]; // MeV

    double J1 = par[2];
    double J2 = par[3];
    double J  = J1 + J2;
    double w  = (2. * J + 1.) / ((2. * J1 + 1.) * (2. * J2 + 1.));

    double li = par[4];
    double lf = par[5];

    double p      = std::sqrt(2. * 14. / 15. * E); // MeV/c
    double h      = 4.14e-15 * 1.0e-6; // MeV.s
    double lambda = (h / p);

    double EFinal = Beam.GetBindingEnergy() - Final.GetBindingEnergy();

    // Width of the entrance channel
    double x_Gp1[1]   = {E};
    double par_Gp1[7] = {Er,
                         res_width,
                         li,
                         (double)Target.GetZ(),
                         (double)Target.GetA(),
                         (double)Beam.GetZ(),
                         (double)Beam.GetA()};
    double Gp1        = f_proton_width(x_Gp1, par_Gp1);

    // Width of the direct 2p channel
    double Gp2 = 0.;
    if (E - EFinal > 0.) {
      double x_Gp2[1] = {E - EFinal};
      // 2*Target because I considered that 2p = 1H + 1H, so Z=2, A=2
      double par_Gp2[6] = {Er - EFinal,
                           lf,
                           2 * (double)Target.GetZ(),
                           2 * (double)Target.GetA(),
                           (double)Final.GetZ(),
                           (double)Final.GetA()};
      Gp2               = f_theoretical_proton_width(x_Gp2, par_Gp2);
    }

    // Width of the sequential 2p channel
    double EIntermediate = Intermediate.GetExcitationEnergy();
    Intermediate.SetExcitationEnergy(IntermediateState);
    double Gp3 = 1.e-3;
    if (E - EIntermediate < 0.) {
      Gp3 = 0.;
    } else if (E - EIntermediate > 0.) {
      double x_Gp3[1]   = {E - EIntermediate};
      double par_Gp3[6] = {Er - EIntermediate,    lf,
                           (double)Target.GetZ(), (double)Target.GetA(),
                           (double)Beam.GetZ(),   (double)Beam.GetA()};
      Gp3               = f_theoretical_proton_width(x_Gp3, par_Gp3);
    }

    double Gt = Gp1 + Gp2 + Gp3;

    double mBeam_uma   = Beam.GetA() + Beam.GetMassExcess() / 1000.;
    double mTarget_uma = Target.GetA() + Target.GetMassExcess() / 1000.;
    double mCN_uma     = (mBeam_uma * mTarget_uma) / (mBeam_uma + mTarget_uma);

    double CN_formation = 0.657 / (mCN_uma * E) * w * (Gp1 * Gt)
                          / (pow((E - Er), 2.) + 1. / 4. * pow(Gt, 2.));
    double CN_decay         = Gp3 / Gt;
    double BW_cross_section = CN_formation * CN_decay;
    return BW_cross_section;
  }
} // Breit-Wigner cross section of the reaction 14O(p,p)14O(p)13N

// "Main of thsd
void ThreeBodyCS(double Er = 6.3, double res_width = 0.02, int Npx = 1000) {

  // Generating BW cross section for the 15F sequential emission, here 14O is
  // taken as the reference, E(14Og.s.) = 0 MeV, as a consequence, the ground
  // state of 13N is 4.628 MeV.

  double J1 = 0.;
  double J2 = 5. / 2.;
  double li = 3.;
  double lf = 2.;

  TF1* g_double_cross_section
      = new TF1("g_double_sequential_cross_section",
                f_double_sequential_proton_cross_section, 0.1, 10.1, 6, 1);

  g_double_cross_section->SetParameters(Er, res_width, J1, J2, li, lf);
  g_double_cross_section->SetNpx(Npx);
  g_double_cross_section->Draw();
  g_double_cross_section->SetName("BW");
  TFile* fout = new TFile("BW.root", "recreate");
  g_double_cross_section->Write();
}
