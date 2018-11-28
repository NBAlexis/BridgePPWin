/*!
        @file    $Id: counter.h #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef COUNTER_INCLUDED
#define COUNTER_INCLUDED

//! Counter of time and GFlops for Blue Gene/Q at KEK.

/*!
   This class wraps libkek.h counter probably available only
   at KEK.
                                   [7 Feb 2014 H.Matsufuru]
 */
class Counter {
 private:
  int        m_id;       //!< counter id.
  static int cid_static; //!< static id to asign unique id to the instance.

 public:
  Counter();
  void start();
  void finish();
  void finish(double& time, double& gflops);
};
#endif /* #ifndef COUNTER_INCLUDED */
