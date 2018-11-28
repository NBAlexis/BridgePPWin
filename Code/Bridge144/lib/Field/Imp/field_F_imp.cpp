#include "BridgeLib_Private.h"
#if USE_IMP

/*!
@file    $Id:: field_F_imp.cpp #$

@brief

@author  Hideo Matsufuru (matsufuru)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 00:49:46 #$

@version $LastChangedRevision: 1561 $
*/

#include "Field/field_F.h"

using std::valarray;

//#include "ResourceManager/threadManager_OpenMP.h"

#if defined USE_GROUP_SU3
#include "field_F_imp_SU3.inc"
#elif defined USE_GROUP_SU2
#include "field_F_imp_SU2.inc"
#elif defined USE_GROUP_SU_N
#include "field_F_imp_SU_N.inc"
#endif

//====================================================================
void Field_F::check()
{
    //  assert(NC == CommonParameters::Nc());
    //  assert(ND == CommonParameters::Nd());

    check_Nc();
}


//====================================================================
void mult_Field_Gn(Field_F& y, const int ex,
    const Field_G& U, const int ex1,
    const Field_F& x, const int ex2)
{
    assert(ex < y.nex());
    assert(ex1 < U.nex());
    assert(ex2 < x.nex());
    assert(U.nvol() == y.nvol());
    assert(x.nvol() == y.nvol());

    const int Nd = x.nd();
    const int Nc = x.nc();
    const int Nc2 = x.nc2();
    const int Ncd2 = Nc2 * Nd;
    const int Ndf = Nc2 * Nc;

    //  int Nthread = ThreadManager_OpenMP::get_num_threads();
    //  int i_thread = ThreadManager_OpenMP::get_thread_id();
    int Nthread = 1;
    int i_thread = 0;

    int is = y.nvol() * i_thread / Nthread;
    int ns = y.nvol() * (i_thread + 1) / Nthread - is;

    double       *v = y.ptr(0, is, ex);
    const double *g = U.ptr(0, is, ex1);
    const double *w = x.ptr(0, is, ex2);

    for (int site = 0; site < ns; ++site) {
        int ig = Ndf * site;
        int iv = Ncd2 * site;
        for (int s = 0; s < Nd; ++s) {
            for (int ic = 0; ic < Nc; ++ic) {
                int ig2 = ic * Nc2 + ig;
                int iv2 = s * Nc2 + iv;
                v[2 * ic + iv2] = mult_Gn_r(&g[ig2], &w[iv2], Nc);
                v[2 * ic + 1 + iv2] = mult_Gn_i(&g[ig2], &w[iv2], Nc);
            }
        }
    }
}


//====================================================================
void mult_Field_Gd(Field_F& y, const int ex,
    const Field_G& U, const int ex1,
    const Field_F& x, const int ex2)
{
    assert(ex < y.nex());
    assert(ex1 < U.nex());
    assert(ex2 < x.nex());
    assert(U.nvol() == y.nvol());
    assert(x.nvol() == y.nvol());

    const int Nd = x.nd();
    const int Nc = x.nc();
    const int Nc2 = x.nc2();
    const int Ncd2 = Nc2 * Nd;
    const int Ndf = Nc2 * Nc;

    //  int Nthread = ThreadManager_OpenMP::get_num_threads();
    //  int i_thread = ThreadManager_OpenMP::get_thread_id();
    int Nthread = 1;
    int i_thread = 0;

    int is = y.nvol() * i_thread / Nthread;
    int ns = y.nvol() * (i_thread + 1) / Nthread - is;

    double       *v = y.ptr(0, is, ex);
    const double *g = U.ptr(0, is, ex1);
    const double *w = x.ptr(0, is, ex2);

    for (int site = 0; site < ns; ++site) {
        int ig = Ndf * site;
        int iv = Ncd2 * site;
        for (int s = 0; s < Nd; ++s) {
            for (int ic = 0; ic < Nc; ++ic) {
                int ig2 = ic * 2 + ig;
                int iv2 = s * Nc2 + iv;
                v[2 * ic + iv2] = mult_Gd_r(&g[ig2], &w[iv2], Nc);
                v[2 * ic + 1 + iv2] = mult_Gd_i(&g[ig2], &w[iv2], Nc);
            }
        }
    }
}


//====================================================================
void multadd_Field_Gn(Field_F& y, const int ex,
    const Field_G& U, const int ex1,
    const Field_F& x, const int ex2,
    const double a)
{
    assert(ex < y.nex());
    assert(ex1 < U.nex());
    assert(ex2 < x.nex());
    assert(U.nvol() == y.nvol());
    assert(x.nvol() == y.nvol());

    const int Nd = x.nd();
    const int Nc = x.nc();
    const int Nc2 = x.nc2();
    const int Ncd2 = Nc2 * Nd;
    const int Ndf = Nc2 * Nc;

    //  int Nthread = ThreadManager_OpenMP::get_num_threads();
    //  int i_thread = ThreadManager_OpenMP::get_thread_id();
    int Nthread = 1;
    int i_thread = 0;

    int is = y.nvol() * i_thread / Nthread;
    int ns = y.nvol() * (i_thread + 1) / Nthread - is;

    double       *v = y.ptr(0, is, ex);
    const double *g = U.ptr(0, is, ex1);
    const double *w = x.ptr(0, is, ex2);

    for (int site = 0; site < ns; ++site) {
        int ig = Ndf * site;
        int iv = Ncd2 * site;
        for (int s = 0; s < Nd; ++s) {
            for (int ic = 0; ic < Nc; ++ic) {
                int ig2 = ic * Nc2 + ig;
                int iv2 = s * Nc2 + iv;
                v[2 * ic + iv2] += a * mult_Gn_r(&g[ig2], &w[iv2], Nc);
                v[2 * ic + 1 + iv2] += a * mult_Gn_i(&g[ig2], &w[iv2], Nc);
            }
        }
    }
}


//====================================================================
void multadd_Field_Gd(Field_F& y, const int ex,
    const Field_G& U, const int ex1,
    const Field_F& x, const int ex2,
    const double a)
{
    assert(ex < y.nex());
    assert(ex1 < U.nex());
    assert(ex2 < x.nex());
    assert(U.nvol() == y.nvol());
    assert(x.nvol() == y.nvol());

    const int Nd = x.nd();
    const int Nc = x.nc();
    const int Nc2 = x.nc2();
    const int Ncd2 = Nc2 * Nd;
    const int Ndf = Nc2 * Nc;

    //  int Nthread = ThreadManager_OpenMP::get_num_threads();
    //  int i_thread = ThreadManager_OpenMP::get_thread_id();
    int Nthread = 1;
    int i_thread = 0;

    int is = y.nvol() * i_thread / Nthread;
    int ns = y.nvol() * (i_thread + 1) / Nthread - is;

    double       *v = y.ptr(0, is, ex);
    const double *g = U.ptr(0, is, ex1);
    const double *w = x.ptr(0, is, ex2);

    for (int site = 0; site < ns; ++site) {
        int ig = Ndf * site;
        int iv = Ncd2 * site;
        for (int s = 0; s < Nd; ++s) {
            for (int ic = 0; ic < Nc; ++ic) {
                int ig2 = ic * 2 + ig;
                int iv2 = s * Nc2 + iv;
                v[2 * ic + iv2] += a * mult_Gd_r(&g[ig2], &w[iv2], Nc);
                v[2 * ic + 1 + iv2] += a * mult_Gd_i(&g[ig2], &w[iv2], Nc);
            }
        }
    }
}


//====================================================================
void mult_GM(Field_F& y, const GammaMatrix& gm, const Field_F& x)
{
    assert(x.nex() == y.nex());
    assert(x.nvol() == y.nvol());

    int Nd = x.nd();
    int Nc = x.nc();
    int Nc2 = x.nc2();
    int Ncd2 = Nc2 * Nd;
    
    int *id = (int *)alloca(sizeof(int) * Nd);
    int *idc_r = (int *)alloca(sizeof(int) * Nd);
    int *idc_i = (int *)alloca(sizeof(int) * Nd);
    double *gv_r = (double *)alloca(sizeof(double) * Nd);
    double *gv_i = (double *)alloca(sizeof(double) * Nd);

    //int    id[Nd];
    //int    idc_r[Nd];
    //int    idc_i[Nd];
    //double gv_r[Nd];
    //double gv_i[Nd];

    
    
    for (int s = 0; s < Nd; ++s) 
    {
        id[s] = gm.index(s);
        gv_r[s] = gm.value_r(s);
        gv_i[s] = gm.value_i(s);
        idc_r[s] = gm.index_c(s);
        idc_i[s] = 1 - idc_r[s];
    }

    //  int Nthread = ThreadManager_OpenMP::get_num_threads();
    //  int i_thread = ThreadManager_OpenMP::get_thread_id();
    int Nthread = 1;
    int i_thread = 0;

    int is = y.nvol() * i_thread / Nthread;
    int ns = y.nvol() * (i_thread + 1) / Nthread - is;

    const int Nex = y.nex();

    for (int ex = 0; ex < Nex; ++ex) 
    {
        double       *v = y.ptr(0, is, ex);
        const double *w = x.ptr(0, is, ex);

        for (int site = 0; site < ns; ++site) 
        {
            int iv = Ncd2 * site;

            for (int s = 0; s < Nd; ++s) 
            {
                int iv2 = s * Nc2 + iv;
                int iw2 = id[s] * Nc2 + iv;

                for (int ic = 0; ic < Nc; ++ic) 
                {
                    v[2 * ic + iv2] = gv_r[s] * w[2 * ic + idc_r[s] + iw2];
                    v[2 * ic + 1 + iv2] = gv_i[s] * w[2 * ic + idc_i[s] + iw2];
                }
            }
        }
    }
}


//====================================================================
void mult_iGM(Field_F& y, const GammaMatrix& gm, const Field_F& x)
{
    assert(x.nex() == y.nex());
    assert(x.nvol() == y.nvol());

    int Nd = x.nd();
    int Nc = x.nc();
    int Nc2 = x.nc2();
    int Ncd2 = Nc2 * Nd;

    int *id = (int *)alloca(sizeof(int) * Nd);
    int *idc_r = (int *)alloca(sizeof(int) * Nd);
    int *idc_i = (int *)alloca(sizeof(int) * Nd);
    double *gv_r = (double *)alloca(sizeof(double) * Nd);
    double *gv_i = (double *)alloca(sizeof(double) * Nd);

    //int    id[Nd];
    //int    idc_r[Nd];
    //int    idc_i[Nd];
    //double gv_r[Nd];
    //double gv_i[Nd];

    for (int s = 0; s < Nd; ++s)
    {
        id[s] = gm.index(s);
        gv_r[s] = -gm.value_i(s);
        gv_i[s] = gm.value_r(s);
        idc_r[s] = 1 - gm.index_c(s);
        idc_i[s] = 1 - idc_r[s];
    }

    //  int Nthread = ThreadManager_OpenMP::get_num_threads();
    //  int i_thread = ThreadManager_OpenMP::get_thread_id();
    int Nthread = 1;
    int i_thread = 0;

    int is = y.nvol() * i_thread / Nthread;
    int ns = y.nvol() * (i_thread + 1) / Nthread - is;

    const int Nex = y.nex();

    for (int ex = 0; ex < Nex; ++ex) 
    {
        double       *v = y.ptr(0, is, ex);
        const double *w = x.ptr(0, is, ex);

        for (int site = 0; site < ns; ++site) 
        {
            int iv = Ncd2 * site;

            for (int s = 0; s < Nd; ++s) 
            {
                int iv2 = s * Nc2 + iv;
                int iw2 = id[s] * Nc2 + iv;
                for (int ic = 0; ic < Nc; ++ic) 
                {
                    v[2 * ic + iv2] = gv_r[s] * w[2 * ic + idc_r[s] + iw2];
                    v[2 * ic + 1 + iv2] = gv_i[s] * w[2 * ic + idc_i[s] + iw2];
                }
            }
        }
    }
}


//====================================================================
void mult_GMproj(Field_F& y,
    const int pm, const GammaMatrix& gm,
    const Field_F& w)
{
    assert(w.nvol() == y.nvol());
    assert(w.nex() == y.nex());

    mult_GM(y, gm, w);

    if (pm == 1) {
        axpy(y, 1.0, w);  // y += w;
        scal(y, 0.5);     // y *= 0.5;
    }
    else if (pm == -1) {
        axpy(y, -1.0, w); // y -= w     i.e. y = (gamma - 1) * w
        scal(y, -0.5);    // y *= -0.5  i.e. y = (1 - gamma)/2 * w
    }
    else {
        vout.crucial("Error at %s: wrong pm.\n", __func__);
        exit(EXIT_FAILURE);
    }
}


//====================================================================
void mult_GMproj2(Field_F& y,
    const int pm, const GammaMatrix& gm,
    const Field_F& w)
{
    assert(w.nvol() == y.nvol());
    assert(w.nex() == y.nex());

    mult_GM(y, gm, w);

    if (pm == 1) {
        axpy(y, 1.0, w);  // y += w;
    }
    else if (pm == -1) {
        axpy(y, -1.0, w); // y -= w;
        scal(y, -1.0);    // y *= -1.0;
    }
    else {
        vout.crucial("Error at %s: wrong pm.\n", __func__);
        exit(EXIT_FAILURE);
    }
}


//====================================================================
void mult_GMproj2(Field_F& y,
    const double nu_s,
    const int pm,
    const double r_s,
    const GammaMatrix& gm,
    const Field_F& w)
{
    assert(w.nvol() == y.nvol());
    assert(w.nex() == y.nex());

    mult_GM(y, gm, w);
    scal(y, nu_s);  // y *= nu_s;

    if (pm == 1) {
        axpy(y, r_s, w);  // y += r_s * w;
    }
    else if (pm == -1) {
        axpy(y, -r_s, w); // y -= r_s * w  i.e. y = (nu_s * gamma - r_s) * w
        scal(y, -1.0);    // y *= -1.0     i.e. y = (r_s - nu_s * gamma) * w
    }
    else {
        vout.crucial("Error at %s: wrong pm = %d\n", __func__, pm);
        exit(EXIT_FAILURE);
    }
}


//====================================================================
//============================================================END=====
#endif //COMMUNICATOR_USE_MPI
