#include "BridgeLib_Private.h"

/*!
@file    $Id: fieldIO_Binary.cpp #$

@brief

@author  Hideo Matsufuru (matsufuru)
$LastChangedBy: aoym $

@date    $LastChangedDate: 2014-06-06 18:09:02 #$

@version $LastChangedRevision: 1561 $
*/

#include "fieldIO_Binary.h"

#include <fstream>
#if BRIDGE_WIN
#include <string.h>
#else
#include <strings.h>
#endif

#include "bridgeIO.h"
using Bridge::vout;

//typedef unsigned short int  n_uint16_t;
//typedef unsigned int        n_uint32_t;
//typedef unsigned long int   n_uint64_t;

const std::string FieldIO_Binary::class_name = "FieldIO_Binary";

//====================================================================
void FieldIO_Binary::read_file(Field *v, std::string filename)
{
    int nin_field = v->nin();
    int nex_field = v->nex();

    int nin_file = m_format->nin();
    int nex_file = m_format->nex();

    if ((nin_file == 0) || (nex_file == 0)) {
        nin_file = nin_field;
        nex_file = nex_field;
    }

    // temporary buffer: allocated only at I/O node.
    Field vtmp;

    if (Communicator::is_primary()) {
        vout.detailed(m_vl, "reading field data from %s\n", filename.c_str());

        int Lvol = CommonParameters::Lvol();
        vtmp.reset(nin_field, Lvol, nex_field);

        bool do_swap = (is_bigendian() == false);
        if (do_swap) {
            vout.detailed(m_vl, "host endian is little: byte swap performed.\n");
        }

        const int block_size = nin_file;

        //char      buf[sizeof(double) * block_size];
        char *buf = (char *)alloca(sizeof(double) * block_size);

        std::ifstream infile(filename.c_str(), std::ios::in | std::ios::binary);
        if (!infile) {
            vout.crucial(m_vl, "Error at %s: file open failed, %s may not exist.\n", class_name.c_str(), filename.c_str());
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < nex_file; ++j) {
            for (int isite = 0; isite < Lvol; ++isite) {
                // read 1 block
                infile.read(buf, sizeof(double) * block_size);

                if (!infile) {
                    if (infile.eof()) {  // short file
                        vout.crucial(m_vl, "Error at %s: file size of %s is too small.\n", class_name.c_str(), __func__);
                    }
                    else {
                        vout.crucial(m_vl, "Error at %s: io error of %s.\n", class_name.c_str(), __func__);
                    }

                    exit(EXIT_FAILURE);
                }

                if (do_swap) {
                    byte_swap(buf, sizeof(double), block_size);
                }

                double *ptr = (double *)buf;

                for (int i = 0; i < nin_file; ++i) {
                    int s, t;
                    m_format->file_to_field(s, t, i, j);

                    vtmp.set(s, isite, t, ptr[i]);
                }
            }
        }

        infile.close();
    }

    FieldIO::deliver(v, &vtmp);

    vout.detailed(m_vl, "read successful\n");
}


//====================================================================
void FieldIO_Binary::write_file(Field *v, std::string filename)
{
    int nin_field = v->nin();
    int nex_field = v->nex();

    int nin_file = m_format->nin();
    int nex_file = m_format->nex();

    if ((nin_file == 0) || (nex_file == 0)) {
        nin_file = nin_field;
        nex_file = nex_field;
    }

    int Lvol = CommonParameters::Lvol();

    Field vtmp;
    if (Communicator::is_primary()) {
        vtmp.reset(nin_field, Lvol, nex_field);
    }

    FieldIO::gather(&vtmp, v);

    if (Communicator::is_primary()) {
        vout.detailed(m_vl, "writing field data to %s\n", filename.c_str());

        bool do_swap = (is_bigendian() == false);
        if (do_swap) {
            vout.detailed(m_vl, "host endian is little: byte swap performed.\n");
        }

        const int block_size = nin_file;
        //char      buf[sizeof(double) * block_size];
        char *buf = (char *)alloca(sizeof(double) * block_size);

        std::ofstream outfile(filename.c_str(), std::ios::out | std::ios::binary);
        if (!outfile) {
            vout.crucial(m_vl, "Error at %s: file open of %s failed\n", class_name.c_str(), filename.c_str());
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < nex_file; ++j) {
            for (int isite = 0; isite < Lvol; ++isite) {
                double *ptr = (double *)buf;

                for (int i = 0; i < nin_file; ++i) {
                    int s, t;
                    m_format->file_to_field(s, t, i, j);

                    ptr[i] = vtmp.cmp(s, isite, t);
                }

                if (do_swap) {
                    byte_swap(buf, sizeof(double), block_size);
                }

                outfile.write(buf, sizeof(double) * block_size);

                if (!outfile) {  // error
                    vout.crucial(m_vl, "Error at %s: io error of %s.\n", class_name.c_str(), __func__);
                    exit(EXIT_FAILURE);
                }
            }
        }

        outfile.close();
    }

    vout.detailed(m_vl, "write succeeded.\n");
}


//====================================================================
//============================================================END=====
