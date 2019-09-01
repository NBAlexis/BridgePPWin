/*!
        @file    io_format.h

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef IO_FORMAT_INCLUDED
#define IO_FORMAT_INCLUDED

//! IO_Format for data layout conversion from/to file.

/**
   Format class BAPI provides data layout conversion between file and field data.

   Field data is structured so that it has inner index, space-time index,
   and outer index. Format class BAPI provides method that maps inner and outer
   indexes between them.

   Format class BAPI is an abstract base class BAPI that defines interfaces.

   A subclass, Trivial Format, provides trivial conversion (the layouts of
   file and field are the same).
   Another subcalss for gauge configuration would be provided separately.

   N.B. A Set of Format classes should be provided for a specific field layout.
 */

namespace IO_Format {
  class BAPI Format {
   public:
    virtual ~Format() {}

    virtual int nin() const = 0;
    virtual int nex() const = 0;

    virtual void file_to_field(int& s, int& t, const int i, const int j) const = 0;
  };

  class BAPI Trivial_Format : public Format {
   public:
    virtual ~Trivial_Format() {}

    virtual int nin() const { return 0; }
    virtual int nex() const { return 0; }

    virtual void file_to_field(int& s, int& t, const int i, const int j) const
    {
      s = i;
      t = j;
    }
  };

//----------------------------------------------------------------
// predefined formats

  BAPI extern const Format *Trivial;

//----------------------------------------------------------------
}
#endif /* IO_FORMAT_INCLUDED */
