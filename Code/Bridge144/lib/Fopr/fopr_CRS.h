/*!
@file    $Id:: fopr_CRS.h #$

@brief

@author  Hideo Matsufuru (matsufuru)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 18:35:38 #$

@version $LastChangedRevision: 1571 $
*/


#ifndef FOPR_CRS_INCLUDED
#define FOPR_CRS_INCLUDED

#include <string>

#include "fopr.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Fermion operator with CRS matrix format.

/*!
This fermion operator class is defined by CRS matrix format
which is widely used in studies of linear algorithms.
The matrix is defined by giving fermion operator, or reading
from a file (filename is provided at construction).
In present implementation, this class works only on single
node.  Parallel version should be implemented.
[07 Dec 2011  H.Matsufuru]
*/

class BAPI Fopr_CRS : public Fopr
{
public:
    static const std::string class_name;

private:
    int                 m_Nin, m_Nvol, m_Nex;
    int                 m_Nsize, m_Nnz;
    std::vector<int>    m_rowidx_nz;
    std::vector<int>    m_column_nz;
    std::vector<double> m_elem_nz;
    std::string         m_mode;
    Fopr                *m_fopr;

public:

    Fopr_CRS(Fopr *fopr)
        : Fopr(), m_fopr(fopr)
    {
        set_matrix();
    }

    Fopr_CRS(unique_ptr<Fopr>& fopr)
        : Fopr(), m_fopr(fopr.get())
    {
        set_matrix();
    }

    Fopr_CRS(std::string fname)
        : Fopr(), m_fopr(0)
    {
        set_matrix(fname);
    }

    void set_parameters(const Parameters&);

    void write_matrix(std::string);

    void set_config(Field *U)
    {
        if (m_fopr == 0) {
            vout.crucial(m_vl, "Error at %s: fopr is not set.\n", class_name.c_str());
            exit(EXIT_FAILURE);
        }
        else {
            m_fopr->set_config(U);
        }
    }

    void set_config(unique_ptr<Field_G>& U)
    {
        if (m_fopr == 0) {
            vout.crucial(m_vl, "Error at %s: fopr is not set.\n", class_name.c_str());
            exit(EXIT_FAILURE);
        }
        else {
            m_fopr->set_config(U.get());
        }
    }

    void set_mode(std::string mode)
    {
        m_mode = mode;
    }

    std::string get_mode() const
    {
        return m_mode;
    }

    void mult(Field& v, const Field& f)
    {
        if (m_mode == "D") {
            D(v, f);
        }
        else if (m_mode == "DdagD") {
            DdagD(v, f);
        }
        else if (m_mode == "Ddag") {
            Ddag(v, f);
        }
        else {
            vout.crucial(m_vl, "Error at %s: mode unknown: '%s'.\n", class_name.c_str(), m_mode.c_str());
            exit(EXIT_FAILURE);
        }
    }

    void mult_dag(Field& v, const Field& f)
    {
        if (m_mode == "D") {
            Ddag(v, f);
        }
        else if (m_mode == "DdagD") {
            DdagD(v, f);
        }
        else if (m_mode == "Ddag") {
            D(v, f);
        }
        else {
            vout.crucial(m_vl, "Error at %s: mode unknown: '%s'.\n", class_name.c_str(), m_mode.c_str());
            exit(EXIT_FAILURE);
        }
    }

    void DdagD(Field&, const Field&);
    void D(Field&, const Field&);
    void Ddag(Field&, const Field&);
    void H(Field&, const Field&);

    int field_nvol() { return m_Nvol; }
    int field_nin() { return m_Nin; }
    int field_nex() { return m_Nex; }

private:
    void set_matrix();
    void set_matrix(std::string);

    void set_matrix_1row(int&, std::vector<int>&,
        std::vector<double>&, Field&);
};
#endif
