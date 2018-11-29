/*!
        @file    $Id: channel.h #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2012-11-15 16:42:21 #$

        @version $LastChangedRevision: 1679 $
*/

#ifndef CHANNEL_INCLUDED
#define CHANNEL_INCLUDED

#include <cassert>

#include "Communicator/communicator.h"

static const int max_dimension = 8;

// forward declaration
class BAPI ChannelSet;

//! Channel class BAPI for asynchronous communication

/**
   Channel class BAPI for single process.
   This implementation acts similarly to the Channel class BAPI with MPI,
   while practically copys data through the buffer.

   Channel class BAPI defines communication channel between sender and receiver
   for asynchronous data transfer. It actually provides a half-channel of
   sender part or receiver part.
   Original version of Channel class BAPI wraps MPI persistent communication,
   as well as send/receive buffer.

   An instance of channel class BAPI is generated through Communicator send_int()
   and recv_init() methods.
   To access the buffer, operator[] is provided as if it is an ordinary
   container or array.

   In single process version, start() method copys data to buffer,
   and wait() method does nothing.
                                                [01 Sep 2017 H.Matsufuru]
*/

class BAPI Channel {
 public:
  typedef char element_type;    //!< data transfer is byte-wise.

  Channel();                //!< constructor.
  Channel(const int count); //!< constructor with buffer size (count bytes)
  Channel(void *buf); //!< constructor with buffer
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
  element_type *m_buf;  //!< buffer
  bool m_my_buf;  //!< whether buffer is owned by this instance.

  friend class BAPI Communicator;
  friend class BAPI ChannelSet;
};

//! ChannelSet class BAPI for a collection of channels

/**
   ChannelSet defines a collection of channel class BAPI instances
   to invoke start and wait methods collectively for a set of channels.
*/
class BAPI ChannelSet {
 public:
  ChannelSet(int nchannel = 8); //!< constructor. default number of channels is 8 for upward and downward in 4 dimensions.

  int append(Channel *const p); //!< append channel to the set. there is no way to remove a channel.

  int start();                  //!< collective start
  int wait();                   //!< collective wait
};
#endif /* _CHANNEL_H_ */
