/*!
@file    $Id:: randomNumbers.h #$

@brief

@author  Hideo Matsufuru (matsufuru)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 18:35:38 #$

@version $LastChangedRevision: 1571 $
*/

#ifndef RANDOMNUMBERS_INCLUDED
#define RANDOMNUMBERS_INCLUDED

#include "Field/field.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Base class of random number generators.

/*!
This class defines the interface of random number
generator, and implements common methods.
Practical methods to generate random numbers are
defined in subclasses.
This class also implements Gaussian random number and
method to set a global field of Gaussian random numbers
and cut it out to the local field for the own node
(gauss_lex_global()) which is useful in HMC etc.
[25 Dec 2011 H.Matsufuru]
*/

class BAPI RandomNumbers
{
public:
    static const std::string class_name;

protected:
    Bridge::VerboseLevel m_vl;

public:

    RandomNumbers()
        : m_vl(CommonParameters::Vlevel()) {}

    virtual ~RandomNumbers() {}

private:
    // non-copyable
    RandomNumbers(const RandomNumbers&);
    RandomNumbers& operator=(const RandomNumbers&);

public:

    void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

    virtual double get() = 0;

    void gauss(double& rand1, double& rand2);


    //! gaussian random number defined on global lattice.
    virtual void gauss_lex_global(Field&);

    //! gaussian noise for even-odd perconditioned field (S.UEDA)
    virtual void gauss_eo_global(Field&);

    //! uniform random number defined on global lattice.
    virtual void uniform_lex_global(Field&);

    //! save and load random number status.
    virtual void read_file(const std::string&) = 0;
    virtual void write_file(const std::string&) = 0;

    //! reset state with new seed.
    virtual void reset(unsigned long seed) = 0;


protected:

    class rand_gauss_even
    {
    public:
        rand_gauss_even(Field& f, RandomNumbers *rand_gauss)
            : m_rand_gauss(rand_gauss), m_ptr(f.ptr(0)), m_block(f.nin()) {}
        void operator()(const bool do_fill);
        size_t block_size() const;

    private:
        RandomNumbers * m_rand_gauss;
        double        *m_ptr;
        size_t        m_block;
    };

    class rand_gauss_odd
    {
    public:
        rand_gauss_odd(Field& f, RandomNumbers *rand_gauss)
            : m_rand_gauss(rand_gauss), m_ptr(f.ptr(0)), m_block(f.nin()) {}
        void operator()(const bool do_fill);
        size_t block_size() const;

    private:
        RandomNumbers * m_rand_gauss;
        double        *m_ptr;
        size_t        m_block;
    };

    class rand_uniform
    {
    public:
        rand_uniform(Field& f, RandomNumbers *rand_gauss)
            : m_rand_gauss(rand_gauss), m_ptr(f.ptr(0)), m_block(f.nin()) {}
        void operator()(const bool do_fill);
        size_t block_size() const;

    private:
        RandomNumbers * m_rand_gauss;
        double        *m_ptr;
        size_t        m_block;
    };

private:
    template<typename InnerGenerator>
    void generate_global(Field& f);
};
#endif
