/*!
        @file    $Id:: io_format_gauge.h #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef IO_FORMAT_GAUGE_INCLUDED
#define IO_FORMAT_GAUGE_INCLUDED

#include "io_format.h"
#include "Parameters/commonParameters.h"

namespace IO_Format {
  namespace Gauge {
/**
  ILDG_Format -- 4-dimensional SU(3) gauge field configuration

  layout: (u_\mu(n))_ab
    complex x b x a x mu x n
*/

    class ILDG_Format : public Format {
     public:
      ILDG_Format()
      {
//    int Nc = CommonParameters::Nc();
//    int Ndim = CommonParameters::Ndim();

        int Nc   = 3;
        int Ndim = 4;

        m_nin    = 2 * Nc * Nc * Ndim;
        m_nex    = 1;
        m_matrix = 2 * Nc * Nc;
      }

      int nin() const { return m_nin; }
      int nex() const { return m_nex; }

      void file_to_field(int& s, int& t, const int i, const int j) const
      {
        s = i % m_matrix;
        t = i / m_matrix;
      }

     private:
      int m_nin;
      int m_nex;
      int m_matrix;
    };

/**
  JLQCD_Format -- JLQCD gauge configuration

  layout: (u_\mu(n))_ab
    complex x b x a x site x mu

  (same as present Field_G format. thus adopt Trivial_Format class.)
 */

    class JLQCD_Format : public Trivial_Format {};

//----------------------------------------------------------------
// predefined formats

    extern const Format *ILDG;
    extern const Format *JLQCD;
  } // namespace Gauge
}   // namespace IO_Format
#endif /* IO_FORMAT_GAUGE_INCLUDED */
