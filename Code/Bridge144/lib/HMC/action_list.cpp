/*!
        @file    $Id:: action_list.cpp #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "action_list.h"

const std::string ActionList::class_name = "ActionList";

//====================================================================
ActionList::ActionList(const int nlevel)
  : m_vl(CommonParameters::Vlevel()),
    m_num_level(nlevel), m_actions(nlevel), m_integrator_types(nlevel)
{
}


//====================================================================
bool ActionList::append(const int level, Action *action)
{
  if ((level < 0) || (level >= m_num_level)) return false;

  m_actions[level].push_back(action);

  return true;
}


//====================================================================
bool ActionList::append(const int level, unique_ptr<Action>& action)
{
  return append(level, action.get());
}


//====================================================================
bool ActionList::append(const int level, const std::vector<Action *>& actions)
{
  if ((level < 0) || (level >= m_num_level)) return false;

  for (size_t i = 0, n = actions.size(); i < n; ++i) {
    m_actions[level].push_back(actions[i]);
  }

  return true;
}


//====================================================================
bool ActionList::set_integrator_type(const int level, const std::string& type)
{
  if ((level < 0) || (level >= m_num_level)) return false;

  m_integrator_types[level] = type;

  return true;
}


//====================================================================
bool ActionList::set_integrator_type(const std::string& type)
{
  for (size_t i = 0; i < m_num_level; ++i) {
    m_integrator_types[i] = type;
  }

  return true;
}


//====================================================================
ActionSet ActionList::get_actions() const
{
  ActionSet v;

  for (size_t i = 0; i < m_num_level; ++i) {
    v.insert(v.end(), m_actions[i].begin(), m_actions[i].end());
  }

  return v;
}


//====================================================================
ActionSet ActionList::get_actions(const int level) const
{
  if ((level < 0) || (level >= m_num_level)) return ActionSet();

  return m_actions[level];
}


//====================================================================
std::string ActionList::get_integrator_type(const int level) const
{
  if ((level < 0) || (level >= m_num_level)) return std::string();

  return m_integrator_types[level];
}


//====================================================================
int ActionList::get_levels() const
{
  return m_num_level;
}


//====================================================================
//============================================================END=====
