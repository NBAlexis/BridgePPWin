/*!
        @file    $Id: channel.h #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2012-11-15 16:42:21 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef CHANNEL_INCLUDED
#define CHANNEL_INCLUDED

#include <vector>
#include <cassert>

#include "Communicator/communicator.h"

static const int max_dimension = 4;

class Channel {
 public:
  typedef char                        element_type;
  typedef std::vector<element_type>   container_type;
  enum channel_mode { SEND, RECV };
 private:

  /*
  element_type* m_buf;
  //  container_type m_buf;
  static std::vector<container_type> m_buf_send;
  static std::vector<container_type> m_buf_recv;
  */
  static std::vector<container_type> m_buf;
  int m_count;
  int m_ibuf;

  struct bgnet_IDs
  {
    int fifoID;
    int sendBufID;
    int sendOffset;
    int size;
    int destRank;
    int groupID;
    int recvBufID;
    int recvOffset;
    int counterID;
  };

  bgnet_IDs              m_bgnet_ids;
  std::vector<bgnet_IDs> m_bgnet_ids_thread;

  static int m_set_id;  //!< identifies set of channels.

 public:

  Channel();
  Channel(const int count);
  Channel(const int count, const int tag, const channel_mode mode);
  virtual ~Channel();

  int start();
  int wait();

  //! To be called before defining new set of channels for a operator.
  static void increment_channel_set_id();

  void set_thread(int, const std::vector<int>&,
                  const std::vector<int>&, const std::vector<int>&);
  int start_thread(int);
  int wait_thread(int);

  inline element_type& operator[](unsigned int idx)
  {
    return m_buf[m_ibuf][idx];
  }

  //  inline element_type operator[] (unsigned int idx) const { return m_buf[idx]; }

  inline const element_type *ptr()
  {
    return &m_buf[m_ibuf][0];
  }

  inline const element_type *ptr(int offset)
  {
    return &m_buf[m_ibuf][offset];
  }

  friend class Communicator_impl;
};
#endif /* CHANNEL_INCLUDED */
