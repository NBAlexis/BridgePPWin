/*!
        @file    dataIO_Text_impl.h

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2014-10-09 17:20:36 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef DATAIO_TEXT_IMPL_INCLUDED
#define DATAIO_TEXT_IMPL_INCLUDED

template<typename T>
void DataIO_Text::read_file_base(T *v, const size_t n, const string& filename)
{
  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "reading data to %s.\n", filename.c_str());

    std::fstream fs(filename.c_str(), std::ios::in);
    if (!fs.is_open()) {
      vout.crucial(m_vl, "Error at dataIO_Text_impl: read file failed.\n");
      exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < n; ++i) {
      fetch_data<T>(v[i], fs);
    }

    fs.close();
  }

  vout.general(m_vl, "read successful.\n");
}


template<typename T>
void DataIO_Text::write_file_base(const T *v, const size_t n, const string& filename, const bool append)
{
  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "writing data to %s.\n", filename.c_str());

    std::fstream fs(filename.c_str(), std::ios::out | (append ? std::ios::app : std::ios::trunc));
    if (!fs.is_open()) {
      vout.crucial(m_vl, "Error at dataIO_Text_impl: write file failed.\n");
      exit(EXIT_FAILURE);
    }

    fs.setf(std::ios_base::scientific, std::ios_base::floatfield);
    fs.precision(m_format_precision);

    for (size_t i = 0; i < n; ++i) {
      store_data<T>(v[i], fs);
    }

    fs.close();
  }

  vout.general(m_vl, "write successful.\n");
}


template<>
inline void
DataIO_Text::store_data<double> (const double& v, fstream& fs)
{
  fs << v << std::endl;
}


template<>
inline void
DataIO_Text::store_data<dcomplex> (const dcomplex& v, fstream& fs)
{
  fs << real(v) << "\t" << imag(v) << std::endl;
}


template<>
inline void
DataIO_Text::fetch_data<double> (double& v, fstream& fs)
{
  double val;

  fs >> val;
  v = val;
}


template<>
inline void
DataIO_Text::fetch_data<dcomplex> (dcomplex& v, fstream& fs)
{
  double val_r, val_i;

  fs >> val_r >> val_i;
  v = cmplx(val_r, val_i);
}


#endif /* DATAIO_TEXT_IMPL_INCLUDED */
