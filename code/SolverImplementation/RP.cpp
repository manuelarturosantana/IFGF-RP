
#include "../solver2.h"

void Solver::beta(const double r_0, const double r_1, const double r_2,
                  const long long q, const int flag_u_loc, const int flag_v_loc, 
                  const double u_a_loc, const double u_b_loc, const double v_a_loc, const double v_b_loc,
                  const double ubar_loc, const double vbar_loc,
                  std::vector<std::complex<double>>& prec)
{

    std::vector<double> ui(Nu_prec_), si(Nu_prec_), muwi(Nu_prec_);
    std::vector<double> vj(Nv_prec_), tj(Nv_prec_), muwj(Nv_prec_);

    for (int i = 0; i < Nu_prec_; i++) {

        ui[i] = xi(ubar_loc, fejer_nodes_u_prec_[i]);
        si[i] = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, eta(ui[i], flag_u_loc));
        muwi[i] = der_xi(ubar_loc, fejer_nodes_u_prec_[i]) * fejer_weights_u_prec_[i];
    
    }
    
    for (int j = 0; j < Nv_prec_; j++) {

        vj[j] = xi(vbar_loc, fejer_nodes_v_prec_[j]);
        tj[j] = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, eta(vj[j], flag_v_loc));
        muwj[j] = der_xi(vbar_loc, fejer_nodes_v_prec_[j]) * fejer_weights_v_prec_[j];
    
    }
    
    std::vector<std::complex<double>> H(Nu_prec_*Nv_prec_);
    
    if (GEOMETRY == 0) {

        for (int i = 0; i < Nu_prec_; i++) {
            for (int j = 0; j < Nv_prec_; j++) {                        
                
                const double s = si[i];
                const double t = tj[j];

                double px, py, pz;
                double nx, ny, nz;
                double dxdsx, dxdsy, dxdsz;
                double dxdtx, dxdty, dxdtz;
                double dxdsdsx, dxdsdsy, dxdsdsz;
                double dxdsdtx, dxdsdty, dxdsdtz;
                double dxdtdtx, dxdtdty, dxdtdtz;

                parametrization_q(s, t, q, px, py, pz);
                normal_q(s, t, q, nx, ny, nz);
                dxds_q(s, t, q, dxdsx, dxdsy, dxdsz);
                dxdt_q(s, t, q, dxdtx, dxdty, dxdtz);
                dxdsds_q(s, t, q, dxdsdsx, dxdsdsy, dxdsdsz);
                dxdsdt_q(s, t, q, dxdsdtx, dxdsdty, dxdsdtz);
                dxdtdt_q(s, t, q, dxdtdtx, dxdtdty, dxdtdtz);                        

                dxdsx *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);
                dxdsy *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);
                dxdsz *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);

                dxdtx *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdty *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdtz *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);

                dxdsdsx *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);
                dxdsdsy *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);
                dxdsdsz *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);

                dxdsdtx *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdsdty *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdsdtz *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);

                dxdtdtx *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdtdty *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdtdtz *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);

                HH(r_0, r_1, r_2,
                px, py, pz,
                nx, ny, nz,
                dxdsx, dxdsy, dxdsz,
                dxdtx, dxdty, dxdtz,
                dxdsdsx, dxdsdsy, dxdsdsz,
                dxdsdtx, dxdsdty, dxdsdtz,
                dxdtdtx, dxdtdty, dxdtdtz,
                coupling_parameter_,
                H[i*Nv_prec_+j]);

                
                H[i*Nv_prec_+j] *= muwi[i]*muwj[j];

            }
        }

    } else {

        InterpPatch elem = interp_surface_[q];

        std::vector<double> px(Nu_prec_*Nv_prec_), py(Nu_prec_*Nv_prec_), pz(Nu_prec_*Nv_prec_);
        std::vector<double> nx(Nu_prec_*Nv_prec_), ny(Nu_prec_*Nv_prec_), nz(Nu_prec_*Nv_prec_);
        std::vector<double> dxdsx(Nu_prec_*Nv_prec_), dxdsy(Nu_prec_*Nv_prec_), dxdsz(Nu_prec_*Nv_prec_);
        std::vector<double> dxdtx(Nu_prec_*Nv_prec_), dxdty(Nu_prec_*Nv_prec_), dxdtz(Nu_prec_*Nv_prec_);

        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.x, si, tj, &px[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.y, si, tj, &py[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.z, si, tj, &pz[0]);

        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.nuX, si, tj, &nx[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.nuY, si, tj, &ny[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.nuZ, si, tj, &nz[0]);

        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dxdu, si, tj, &dxdsx[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dydu, si, tj, &dxdsy[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dzdu, si, tj, &dxdsz[0]);

        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dxdv, si, tj, &dxdtx[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dydv, si, tj, &dxdty[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dzdv, si, tj, &dxdtz[0]);

        for (int i = 0; i < Nu_prec_; i++) {
            for (int j = 0; j < Nv_prec_; j++) {

                HH(r_0, r_1, r_2,
                px[i*Nv_prec_+j], py[i*Nv_prec_+j], pz[i*Nv_prec_+j],
                nx[i*Nv_prec_+j], ny[i*Nv_prec_+j], nz[i*Nv_prec_+j],
                dxdsx[i*Nv_prec_+j] * der_eta(si[i], flag_u_loc) * 0.5 * (u_b_loc - u_a_loc), dxdsy[i*Nv_prec_+j] * der_eta(si[i], flag_u_loc) * 0.5 * (u_b_loc - u_a_loc), dxdsz[i*Nv_prec_+j] * der_eta(si[i], flag_u_loc) * 0.5 * (u_b_loc - u_a_loc),
                dxdtx[i*Nv_prec_+j] * der_eta(tj[j], flag_v_loc) * 0.5 * (v_b_loc - v_a_loc), dxdty[i*Nv_prec_+j] * der_eta(tj[j], flag_v_loc) * 0.5 * (v_b_loc - v_a_loc), dxdtz[i*Nv_prec_+j] * der_eta(tj[j], flag_v_loc) * 0.5 * (v_b_loc - v_a_loc),
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                coupling_parameter_,
                H[i*Nv_prec_+j]);
                
                H[i*Nv_prec_+j] *= muwi[i]*muwj[j];

            }
        }

    }

    std::vector<std::complex<double>> Tn_mat(Nu_prec_*Nu_int_), Tm_mat(Nv_prec_*Nv_int_);
    
    cheb_evals<Nu_prec_, Nu_int_>(ui, Tn_mat);
    cheb_evals<Nv_prec_, Nv_int_>(vj, Tm_mat);

    std::vector<std::complex<double>> matprod1(Nv_prec_*Nu_int_, {0.0, 0.0});

    std::complex<double> one(1.0, 0.0);
    std::complex<double> zero(0.0, 0.0);

    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                Nu_int_, Nv_prec_, Nu_prec_, 
                &one, &Tn_mat[0], Nu_prec_, &H[0], Nv_prec_, &zero, &matprod1[0], Nv_prec_);

    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                Nu_int_, Nv_int_, Nv_prec_,
                &one, &matprod1[0], Nv_prec_, &Tm_mat[0], Nv_prec_, &zero, &prec[0], Nv_int_);
    
}

