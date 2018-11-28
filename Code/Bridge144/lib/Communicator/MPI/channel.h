/*!
        @file    $Id: channel.h #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2012-11-15 16:42:21 #$

        @version $LastChangedRevision: 1678 $
*/

#ifndef CHANNEL_INCLUDED
#define CHANNEL_INCLUDED

#include <mpi.h>
#include <cassert>
#include <vector>
#include <valarray>

#include "Communicator/communicator.h"

static const int max_dimension = 8;

// forward declaration
class ChannelSet;

//! Channel class for asynchronous communication

/**
   Channel class defines communication channel between sender and receiver
   for asynchronous data transfer. It actually provides a half-channel of
   sender part or receiver part.
   This class wraps MPI persistent communication, as well as send/receive
   buffer.

   An instance of channel class is generated through Communicator send_int()
   and recv_init() methods.
   To access the buffer, operator[] is provided as if it is an ordinary
   container or array.

   start() method to start async data transfer, and wait() method to wait
   for the completion of operation. the channel instance will be re-used
   repeatedly for the same communication channel once it is created.
*/

class Channel {
 public:
  typedef char                          element_type;    //!< data transfer is byte-wise.
  typedef std::valarray<element_type>   container_type;

  Channel();                //!< constructor.
  Channel(const int count); //!< constructor with buffer size (count bytes)
  virtual ~Channel();       //!< destructor

  int start();              //!< start asynchronous communication
  int wait();               //!< wait for completion

  //! accessor to buffer
  inline element_type& operator[](unsigned int idx) { return m_buf[idx]; }
  //! accessor to buffer
  inline element_type operator[](unsigned int idx) const { return m_buf[idx]; }
  //! accessor to buffer; returns pointer to the first element.
  inline const element_type *ptr() { return &m_buf[0]; }

 private:
  container_type m_buf;   //!< buffer

  MPI_Request m_request;  //!< handler to MPI persistent communication
  MPI_Status  m_status;   //!< handler to MPI status information

  friend class Communicator_impl;
  friend class ChannelSet;
};

//! ChannelSet class for a collection of channels

/**
   ChannelSet defines a collection of channel class instances
   to invoke start and wait methods collectively for a set of channels.
*/
class ChannelSet {
 public:
  ChannelSet(int nchannel = 8); //!< constructor. default number of channels is 8 for upward and downward in 4 dimensions.

  int append(Channel *const p); //!< append channel to the set. there is no way to remove a channel.

  int start();                  //!< collective start
  int wait();                   //!< collective wait

 private:
  std::vector<MPI_Request> m_array;  //!< a collection of MPI request held in channels.
  std::vector<MPI_Status> m_status;  //!< a collection of MPI status.
  unsigned int             m_nreq;   //!< number of channels to hold.
};
#endif /* _CHANNEL_H_ */
