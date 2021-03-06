/*!
        @file    action_list.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef ACTION_LIST_INCLUDED
#define ACTION_LIST_INCLUDED

#include "Action/action.h"

#include "IO/bridgeIO.h"

/*! \class BAPI ActionList
    \brief lists of actions at respective integrator levels.

    ActionList represents an array of lists of actions, each list holds
    a set of actions for the corresponding integrator level.

    (the pointer to) an action is stored at a specified level by calling
    append() method.
    the element_type of integrator for individual level or all levels is set
    by set_integrator_type() method.

    ActionList is to be passed to Builder_Integrator as an argument.

    unique_ptr is introduced to avoid memory leaks
                                         [21 Mar 2015 Y.Namekawa]
 */

typedef std::vector<Action *> ActionSet;

class BAPI ActionList
{
 public:
  static const std::string class_name;

 public:
  ActionList(const int nlevel);

 public:
  bool append(const int level, Action *action);
  bool append(const int level, unique_ptr<Action>& action);
  bool append(const int level, const std::vector<Action *>& actions);

  bool set_integrator_type(const int level, const std::string& type);
  bool set_integrator_type(const std::string& type);  // for all levels

  ActionSet get_actions() const;
  ActionSet get_actions(const int level) const;
  std::string get_integrator_type(const int level) const;
  int get_levels() const;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  int m_Nlevels;
  std::vector<ActionSet> m_actions;
  std::vector<std::string> m_integrator_types;
};
#endif
