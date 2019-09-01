/*!
        @file    unique_pointer.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/


#ifndef UNIQUE_POINTER_INCLUDED
#define UNIQUE_POINTER_INCLUDED

#include <cstddef>
#include <cassert>
#define _assert(expr)    assert(expr)

namespace Bridge {
  template<typename T>
  class BAPI unique_ptr  // noncopyable
  {
   public:
    typedef T * pointer;
    typedef T   element_type;

   private:
    pointer _ptr;

    unique_ptr(unique_ptr const&);
    unique_ptr& operator=(unique_ptr const&);

    typedef unique_ptr<T> this_type;

    // safe-bool idiom to avoid unintentional conversion
    typedef element_type * (this_type::*unspecified_bool_type)() const;

   public:

    explicit
    unique_ptr(pointer p = pointer()) : _ptr(p) {}

    ~unique_ptr()
    {
      if (_ptr != pointer()) delete _ptr;
      _ptr = pointer();
    }

    void reset(pointer p = pointer())
    {
      _assert(p == pointer() || p != _ptr);
      this_type(p).swap(*this);
    }

    element_type& operator*() const
    {
      _assert(_ptr != pointer());
      return *_ptr;
    }

    pointer operator->() const
    {
      _assert(_ptr != pointer());
      return _ptr;
    }

    pointer get() const
    {
      return _ptr;
    }

    pointer release()
    {
      pointer _p = _ptr;

      _ptr = pointer();
      return _p;
    }

    void swap(unique_ptr& _u)
    {
      pointer tmp = _u._ptr;

      _u._ptr = _ptr;
      _ptr    = tmp;
    }

    // operator bool() const {
    //   return _ptr == pointer() ? false : true;
    // }

    operator unspecified_bool_type() const {
      return _ptr == 0 ? 0 : &this_type::get;
    }
  };

  template<typename T>
  class BAPI unique_ptr<T[]>  // noncopyable
  {
   public:
    typedef T *pointer;
    typedef T  element_type;

   private:
    pointer _ptr;

    unique_ptr(unique_ptr const&);
    unique_ptr& operator=(unique_ptr const&);

    typedef unique_ptr<T[]> this_type;

   public:

    explicit
    unique_ptr(pointer p = pointer()) : _ptr(p) {}

    ~unique_ptr()
    {
      if (_ptr != pointer()) delete [] _ptr;
      _ptr = pointer();
    }

    void reset(pointer p = pointer())
    {
      _assert(p == pointer() || p != _ptr);
      this_type(p).swap(*this);
    }

    element_type& operator[](std::ptrdiff_t i) const
    {
      _assert(_ptr != pointer());
      _assert(i >= 0);
      return _ptr[i];
    }

    pointer get() const
    {
      return _ptr;
    }

    operator bool() const {
      return _ptr == pointer() ? false : true;
    }

    void swap(unique_ptr& _u)
    {
      pointer tmp = _u._ptr;

      _u._ptr = _ptr;
      _ptr    = tmp;
    }
  };
} // end namespace Bridge

#undef _assert
#endif /* UNIQUE_POINTER_INCLUDED */
